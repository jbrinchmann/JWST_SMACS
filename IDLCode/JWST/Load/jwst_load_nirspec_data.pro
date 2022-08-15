

   

function jwst_load_nirspec_data, DIR, root, border=border, $
                                 blank_negative=blank_negative, $
                                 keep_bad=keep_bad
   ;;
   ;; Given a root like jw02736-o008_s10635
   ;; this finds the NIRSPEC data for that object. 
   ;;
   ;; By default this loads the 1D and 2D spectra 
   ;;

   if (n_elements(border) eq 0) then $
    border = 15

   
   oneD_files = findfile(DIR+root+'*x1d.fits', count=N_1D)
   twoD_files = findfile(DIR+root+'*s2d.fits', count=N_2D)
   cal_files = findfile(DIR+root+'*cal.fits')

   ;; Add reading of optional bad pixel files
   bad_files = findfile(DIR+'bad-*'+root+'*.h5', count=N_bad)
   if (N_1d ne N_2d) then stop
   if (N_bad gt 0 and N_bad ne N_1d) then stop

   
   ;; I sort the filenames here - they will then be sequential in
   ;; filters. This is insufficient when there is a mixture of high
   ;; and low resolution data which overlap (but then this automatic
   ;; handling is not desirable anyhow
   si = sort(oneD_files)
   oneD_files = oneD_files[si]
   twoD_files = twoD_files[si]
   cal_files = cal_files[si]
   if (N_bad gt 0) then $
    bad_files = bad_files[si]
   
   
   filters = strarr(N_1D)
   gratings = strarr(N_1D)
   Rmin = fltarr(N_1d)
   Rmax= fltarr(N_1d)
   L_min = fltarr(N_1d)
   L_max = fltarr(N_1d)
   
   all_obs = Dictionary()

   for i=0L, n_elements(oneD_files)-1 do begin

      info = jwst_parse_filename(oneD_files[i])
      if size(info, /tname) ne 'STRUCT' then begin
         ;; Name did not match expectations so no info
         have_info = 0
      endif else begin
         filters[i] = info.filter
         gratings[i] = info.grating
         info_grating = jwst_get_grating_data(info=info)
         have_info = 1
      endelse

;      print, 'Loading '+oneD_files[i]
      sp_1d = jwst_load_file(oneD_files[i])
      if not have_info then begin
         h = sp_1d.primary_header ; just to save on typing
         info = {progid: long(sxpar(h, 'PROGRAM')), $
                 progIDstr: sxpar(h, 'PROGRAM'), $
                 grating: sxpar(h, 'GRATING'), $
                 filter: sxpar(h, 'FILTER')}
      endif
      
      if (n_bad gt 0) and not keyword_set(keep_bad) then begin
         ;; SHOULD REALLY HAVE A CHECK HERE FOR FILTER ETC!!!
;         print, 'Loading '+bad_files[i]
         s_bad = hdf5_to_struct(bad_files[i])
;         stop
         if (n_elements(s_bad.bad) eq 0) then stop

         if (n_elements(s_bad.bad) gt 0) then begin
            sp_1d.flux[s_bad.bad] = !values.f_nan
            sp_1d.flux_error[s_bad.bad] = !values.f_nan
         endif
            
      endif
;      stop

      
      ;; The units by default are Janskys but since platefit wants
      ;; something per Angstrom as input we want to change here. I
      ;; will transform to erg/s/cm/AA therefore.
      ;;
      conv = 2.99792458E-05/(1e4*sp_1d.wavelength)^2
      conv = conv/1e-10         ; there seems to be a problem with the units
      sp_1d.flux = sp_1d.flux*conv
      sp_1d.flux_error = sp_1d.flux_error*conv
      

      ;; Finally, I would like to estimate an empirical uncertainty
      emp_sig = empirical_sigma(sp_1d.wavelength, sp_1d.flux)
      sp_1d = struct_add_key(sp_1d, 'dflux_emp', emp_sig)
      
      ;; Then add to this. We want the resolution spectrum
      ;; attached to the spectrum
      R = interpol(info_grating.curve.r, info_grating.curve.wavelength, $
                   sp_1d.wavelength)
      sp_1d = struct_add_key(sp_1d, 'R', R)

      ;; Then we also want FWHM and delta log lambda - in AA! which is what
      ;; platefit wants
      fwhm = sp_1d.wavelength/R
      sigma = fwhm/sqrt(8.0*alog(2.0))
      delta_logl = sigma/(alog(10.0)*(sp_1d.wavelength))
      sp_1d = struct_add_key(sp_1d, 'FWHM', fwhm)
      sp_1d = struct_add_key(sp_1d, 'DELTA_LOGL', delta_logl)
;      stop

      ;; And then let's get the 2D data as well.
      sp_2d = jwst_load_file(twoD_files[i])
      sp_2d = create_struct(sp_2d, 'wavelength', sp_1d.wavelength)

      ;;
      ;; We also want to convert the 2D array into cgs/arcsec^2
      ;;
      area = !radeg^2*(60.0*60.0)^2 ; arcsec per steradian
      dims = size(sp_2d.flux, /dimen)
      for i2d=0L, dims[1]-1 do begin
         sp_2d.flux[*, i2d] = sp_2d.flux[*, i2d]*conv/area
         sp_2d.dflux[*, i2d] = sp_2d.dflux[*, i2d]*conv/area
      endfor
;      stop

      if keyword_set(blank_negative) then begin
         bad = where(sp_1d.flux lt -1e4, n_bad_blank)
         if (n_bad_blank gt 0) then begin
            sp_1d.flux[bad] = !values.f_nan
            sp_1d.flux_error[bad] = 0.0
            sp_2d.flux[bad, *] = !values.f_nan
         endif
      endif

      ii = where(sp_1d.flux ne 0.0 and finite(sp_1d.flux)) ; 0.0 is missing data. 
      l_min[i] = min(sp_1d.wavelength[ii])
      l_max[i] = max(sp_1d.wavelength[ii])


      all_obs[filters[i]] = {oned: sp_1d, twod: sp_2d}

      
   endfor

   

   
   ;; If there are multiple gratings, we need to glue the spectra
   ;; together.  There are some overlap regions. I will assume they
   ;; are pairwise and sequential. Given the sorting we made earlier
   ;; the filters are sorted

   Nf = n_elements(filters)

;   stop
   if Nf eq 1 then begin
      sp = all_obs[filters[0]].oned

      l_final = sp.wavelength
      f_final = sp.flux
      df_final = sp.flux_error
      df_emp_final = sp.flux_error
      R_final = sp.R
      delta_logl_final = sp.delta_logl

   endif else begin
   
      ratio = fltarr(Nf-1)
      for i=0L, Nf-1 do begin
         ;; Get this and the subsequent spectrum. We look for overlaps at
         ;; the red edge of the blue (spB) and blue end of red (spR)
         spB = all_obs[filters[i]]
         x1 = spB.oned.wavelength
         y1 = spB.oned.flux
         ;; I never trust the edges of spectra so I have a border.
         high_blue = where(abs(x1-l_max[i]) lt 1e-3)
         high_blue = high_blue[0]-border ; These are in pixels
         l_high_blue = x1[high_blue]

         
         ;; Then the same with the spectrum to the right unless we are
         ;; at the end.
         
         if i < NF-1 then begin
            spR = all_obs[filters[i+1]]
            x2 = spR.oned.wavelength
            y2 = spR.oned.flux
            
            low_red = where(abs(x2 -l_min[i+1]) lt 1e-3)
            low_red = low_red[0]+border
            l_low_red = x2[low_red]
            overlapR = where(x2 lt l_high_blue and x2 gt l_low_red, n_overlapR)

            overlapB = where(x1 gt l_low_red and x1 lt l_high_blue, n_overlapB)

            ;;
            ;; Here a decision needs to be taken. Do we want to average
            ;; together? The dispersion solutions may vary significantly so I
            ;; will instead keep the blue as long as possible and then
            ;; insert the red. However I would like to make a check on the
            ;; flux levels
            ;;
            if n_overlapB gt 0 then $
             ratio[i] = median(y1[overlapB]/y2[overlapR]) $
            else ratio[i] = -9999.
         endif else begin
            overlapB = where(x1 lt l_high_blue, n_overlapB)

         endelse



         
         ;;
         ;; Now we deal with assembling the full wavelength axis but this
         ;; needs different handling if we are the first round through.
         ;;
         if i eq 0 then begin
            keep = where(x1 le l_high_blue)
            
            l_final = x1[keep]
            f_final = y1[keep]
            df_final = spB.oned.flux_error[keep]
            df_emp_final = spB.oned.dflux_emp[keep]
            R_final = spB.oned.R[keep]
            delta_logl_final = spB.oned.delta_logl[keep]

            twod_f = spB.twod.flux[keep, *]
            twod_df = spB.twod.dflux[keep, *]

;            stop
         endif else begin
            ;; We need to cut on the left and the right
            keep = where((x1 gt max(l_finaL)) and (x1 le l_high_blue))
            
            l_final = [l_final, x1[keep]]
            f_final = [f_final, y1[keep]]
            df_final = [df_final, spB.oned.flux_error[keep]]
            df_emp_final = [df_emp_final, spB.oned.dflux_emp[keep]]
            R_final = [R_final, spB.oned.R[keep]]
            delta_logl_final = [delta_logl_final, spB.oned.delta_logl[keep]]


            dims_old = size(twod_f, /dimen)
            newim = fltarr(dims_old[0]+n_elements(keep), dims_old[1])
            newdim = fltarr(dims_old[0]+n_elements(keep), dims_old[1])

            newim[0:dims_old[0]-1, *] = twod_f
            newdim[0:dims_old[0]-1, *] = twod_df

            newim[dims_old[0]:*, *] = spB.twod.flux[keep, *]
            newdim[dims_old[0]:*, *] = spB.twod.dflux[keep, *]
;            stop
            
            twod_f = newim
            twod_df = newdim

            
         endelse
         
      endfor

   endelse


   ;; Mimick the format of the oned/twod division we have for the
   ;; individual filters.
   oned = {wavelength: l_final, flux: f_final, flux_error: df_final, $
           dflux_emp: df_emp_final, $
           dflux: df_final, R: R_final, delta_logl: delta_logl_final, $
           inst_res: delta_logl_final}
   twod = {wavelength: l_final, flux: twod_f, dflux: twod_df, $
           flux_error: twod_df}
   sp_final = {oned: oned, twod: twod}

   all_obs['combined'] = sp_final

   return, all_obs
   
end
