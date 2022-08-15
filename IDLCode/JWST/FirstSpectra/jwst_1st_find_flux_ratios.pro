pro jwst_1st_find_flux_ratios, plot=plot, ps=ps, $
                               use_idlphot=use_idlphot, $
                               use_sextractor=use_sextractor
   ;;
   ;; Find flux ratios relative to NIRCam and rectify spectra.
   ;;
   ;;

   
   ROOT = '/data2/jarle/JWST/ERO/SMACS/'
   SPECDIR = ROOT+'Spectra/'

   t = mrdfits(ROOT+'match-ids.fits', 1)
   zcat = mrdfits(ROOT+'redshift-for-sources.fits', 1)
   spm = hdf5_to_struct('/data2/jarle/JWST/ERO/SMACS/raw-specmags.h5')


   ;;; Need to add to this because the photometry is missing for the
   ;;; following objects:
   ;;
   ;;  SPECID          I_ACS      I_WF3        ZZ
   ;; 2736_2038         5313        1214     -9.99900
   ;; 2736_5735         2197         196     0.222480
   ;; 2736_7677           -1        1218      4.24909
   ;; 2736_8277         4644        1066     -9.99900
   ;; 2736_8717         4103         828     0.786000
   ;; 2736_8730         4142          -1     -9.99900
   ;; 2736_8883           -1          -1     -9.99900
   ;; 2736_8886         3703         788     -9.99900
   ;; 2736_9483         3225         586      1.16162
   ;; 2736_9721           -1         468      2.11809
   ;; 2736_10380          -1          -1     -9.99900
   ;; 2736_10380          -1          -1     -9.99900
   ;; 2736_10444        1940         145     -9.99900
   ;;
   ;; 8140 also clearly has detections in other bands,

   ;; The original catalogues
   nc_mag = hdf5_to_struct(ROOT+'NIRCam-photometry-NIRSpec.h5')

   ;; My IDL aper photometry
   aper = mrdfits(ROOT+'Stamps/AperturePhotometry/photometry_idl.fits', 1)
   aper = arr_of_struct_to_struct_of_arr(aper)

   ;; The SExtractor photometry - this is in separate files
   sex_f200w = mrdfits(ROOT+'Stamps/SExtractor/phot-f200w-aligned.fits', 1)
   sex_f277w = mrdfits(ROOT+'Stamps/SExtractor/phot-f277w-aligned.fits', 1)
   sex_f356w = mrdfits(ROOT+'Stamps/SExtractor/phot-f356w-aligned.fits', 1)
   sex_f444w = mrdfits(ROOT+'Stamps/SExtractor/phot-f444w-aligned.fits', 1)
   
   
   dims = size(nc_mag.aper50_flux_f200w, /dimen)
   nobj = dims[0]
   
   alpha = dblarr(4, nobj)
   alpha_p16 = alpha
   alpha_p84 = alpha

   diff = dblarr(4, nobj)
   diff_p16 = alpha
   diff_p84 = alpha

   ;; Pivot wavelengths from https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-filters
   lpivot = [1.989, 2.762, 3.568, 4.408]

   if keyword_set(plot) and keyword_set(ps) then begin
      OUTDIR = '/data2/jarle/JWST/ERO/SMACS/Figures/'
      if keyword_set(use_idlphot) then suffix = '-idlphot' $
      else  if keyword_set(use_sextractor) then suffix = '-sextractor' $
      else suffix = ''
      ps_on, OUTDIR+'illustration-plots-for-flux-normalisation'+suffix+'.ps', /times, /color, /landscape
      !p.font = 0
      !x.thick = 2
      !y.thick = 2
   endif

   
   for i=0L, nobj-1 do begin


      
      
;      if i le 33 then continue

      
      dum = strsplit(t[i].specid, '_', /extract)
      propid = long(dum[0])
      objid = long(dum[1])

      if objid ne 4580 then continue

      
;      if objid ne 4590 then continue

      if t[i].nircam_f200w_id lt 0 then begin
         ;; Skip if no match exists
         alpha[*, i] = !values.f_nan
         alpha_p16[*, i] = !values.f_nan
         alpha_p84[*, i] = !values.f_nan


         diff[*, i] = !values.f_nan
         diff_p16[*, i] = !values.f_nan
         diff_p84[*, i] = !values.f_nan


         if not keyword_set(use_idlphot) and not keyword_set(use_sextractor) then $
          continue
      endif 

      ;; Get fluxes from NIRCAM

      if keyword_set(use_idlphot) then begin
         ;; Get the fluxes and scale to the provided catalog fluxes.
         i_obj = where(aper.source eq t[i].specid, n_match)
         if (n_match eq 0) then stop
         
         f = [aper.f200w[i_obj], aper.f277w[i_obj], aper.f356w[i_obj], aper.f444w[i_obj]]*8.3359649e-9 
         df = [aper.f200w[i_obj], aper.f277w[i_obj], aper.f356w[i_obj], aper.f444w[i_obj]]*8.3359649e-9/10. ; Arbitrary uncertainty!
      endif else if keyword_set(use_sextractor) then begin
         i_obj = where(sex_f200w.specid eq objid)
         f = [sex_f200w[i_obj].flux_aper, sex_f277w[i_obj].flux_aper, $
              sex_f356w[i_obj].flux_aper, sex_f444w[i_obj].flux_aper]*8.3359649e-9 
         df = [sex_f200w[i_obj].fluxerr_aper, sex_f277w[i_obj].fluxerr_aper, $
               sex_f356w[i_obj].fluxerr_aper, sex_f444w[i_obj].fluxerr_aper]*8.3359649e-9/10.0
      endif else begin
         f = [nc_mag.aper50_flux_f200W[i], nc_mag.aper50_flux_f277w[i], $
              nc_mag.aper50_flux_f356w[i], nc_mag.aper50_flux_f444w[i]]
         
         df = [nc_mag.aper50_flux_err_f200W[i], nc_mag.aper50_flux_err_f277w[i], $
               nc_mag.aper50_flux_err_f356w[i], nc_mag.aper50_flux_err_f444w[i]]
         
         mags_nc = [nc_mag.aper50_abmag_f200W[i], nc_mag.aper50_abmag_f277w[i], $
                    nc_mag.aper50_abmag_f356w[i], nc_mag.aper50_abmag_f444w[i]]
         dmags_nc = [nc_mag.aper50_abmag_err_f200W[i], nc_mag.aper50_abmag_err_f277w[i], $
                     nc_mag.aper50_abmag_err_f356w[i], nc_mag.aper50_abmag_err_f444w[i]]
      endelse
      
      ok_nc = where(f gt 0, n_ok_nc)
      if n_ok_nc eq 0 then stop
      
      alpha[ok_nc, i] = f/reform(spm.flux[i, ok_nc])-1.0        
      diff[ok_nc, i] = f-spm.flux[i, ok_nc]

      ;;
      ;; Unfortunately the various galaxies require quite different
      ;; handling... So here we go on a case-by-case basis. 
      ;;

      ;;
      ;; Load the spectrum and figure out the blue and red grisms
      ;;
      d = jwst_1st_get_coadded(objid, blank=0)
      diffR=d.oned.R-shift(d.oned.R, 1)
      diffR[0] = 0
      big_jump = where(abs(diffR) gt 100, n_jump)
      if (n_jump eq 0) or (n_jump gt 1) then stop
      l_jump = d.oned.wavelength[big_jump[0]]
      blue = where(d.oned.wavelength lt l_jump, complement=red)
      
      ;; We need fNu for this normalisation - since I will renormalise
      ;; below I am leaving out the normalisation factor (units &
      ;; speed of light)
      fnu_orig = d.oned.flux*d.oned.wavelength^2
      dfnu_orig = d.oned.flux_error*d.oned.wavelength^2


      
      ;; We will also normalise to the F277W pivot wavelength as long
      ;; as that is detected
      if f[1] gt 0 then i_pivot = 1 else i_pivot = min(where(f gt 0))
      

      ;;
      ;; Case-by-case settings
      ;;
      check = 0
      norm_range = 150.0/1e4    ; micron
      case objid of
         1370: begin
            scale_type = 'single'
            i_pivot = 0
            yrange = [-5e-8, 5e-8]
         end
         1679: begin
            scale_type = 'spline'
            i_pivot = 1
            yrange = [-5e-8, 1e-7]
         end
         1917: begin
            scale_type = 'spline'
            i_pivot = 1
            yrange = [-1e-7, 3e-7]
         end
         2572: begin
            scale_type = 'shift_linear'
            i_pivot = 1
            yrange = [-1e-6, 2e-6]
         end
         2653: begin
            scale_type = 'linear'
            i_pivot = 1
            yrange = [0, 1e-6]
            check = 1
         end
         3042: begin
            scale_type = 'linear'
            i_pivot = 1
            yrange = [0, 2e-6]
         end
         3772: begin
            scale_type = 'linear'
            i_pivot = 1
            yrange = [-3e-7, 6e-7]
            check = 1
         end
         4580: begin
            scale_type = 'single'
            i_pivot = 0
            yrange = [-2e-8, 2e-8]
         end
         4590: begin
            scale_type = 'linear'
            i_pivot = 2
            yrange = [-3e-7, 1e-6]
            check = 1
         end
         4798: begin
            scale_type = 'linear'
            i_pivot = 1
            yrange = [-5e-8, 2e-7]
            check = 1
         end
         5144: begin
            if keyword_set(use_idlphot) then $
             scale_type = 'linear' $
            else $
             scale_type = 'none'
            i_pivot = 2
            yrange = [-2e-7, 1e-6]
            check = 1
         end
         5735: begin
            scale_type = 'linear'
            i_pivot = 1
         end
         6113: begin
            scale_type = 'single'
            i_pivot = 2
            yrange = [-0.5e-8, 2e-8]
            check = 1
         end
         6355: begin
            scale_type = 'linear'
            i_pivot = 2
            yrange = [-2e-7, 5e-6]
         end
         7570: begin
            ;; This is actually quite bad but has no redshift
            scale_type = 'spline'
            i_pivot = 2
            yrange = [-1e-6, 3e-6]
            check = 1
         end
         8140: begin
            scale_type = 'single'
            i_pivot = 0
            yrange = [-1e-7, 6e-7]
         end
         8311: begin
            scale_type = 'none'
            i_pivot = 2
            yrange = [-1e-7, 6e-7]
         end
         8498: begin
            scale_type = 'single'
            i_pivot = 0
            yrange = [-1e-7, 2e-7]
         end
         8506: begin
            scale_type = 'linear'
            i_pivot = 2
            yrange = [-1e-6, 2e-6]
         end
         9239: begin
            scale_type = 'single'
            if keyword_set(use_idlphot) then $
             i_pivot = 1 $
            else $
             i_pivot = 0
            yrange = [0, 2e-7]
         end
         9922: begin
            scale_type = 'single'
            i_pivot = 0
            yrange = [0, 1e-6]
         end
         10511: begin
            ;; Gives more or less the same result whether linear or
            ;; spline, in either case rather bad
            scale_type = 'linear'
            
            i_pivot = 1
            yrange = [-1e-7, 4e-7]
            check = 1
         end
         10612: begin
            scale_type = 'linear'
            i_pivot = 2
            yrange = [-1e-6, 2e-6]
         end
         
      endcase

      
      

      
      
      dum = min(abs(d.oned.wavelength-lpivot[i_pivot]), i_norm)
      l_norm = d.oned.wavelength[i_norm]
      
      ii_norm = where(abs(d.oned.wavelength-l_norm) lt norm_range and $
                      finite(d.oned.flux) eq 1 and d.oned.flux ne 0.0, n_norm)
      fnu_orig_norm = median(fnu_orig[ii_norm])
      orig_scale = f[i_pivot]/fnu_orig_norm
      fnu_orig_scaled = fnu_orig*orig_scale
      dfnu_orig_scaled = dfnu_orig*orig_scale
      
      stop
      
      alpha_pred = d.oned.wavelength*0.0
      shift = 0
      case scale_type of
         'single': begin
            this_alpha = mean(alpha[ok_nc, i])
            alpha_pred[*] = this_alpha
         end
         'linear': begin
            sixlin, lpivot[0:1], alpha[0:1,i], a_b, siga_b, b_b, sigb_b
            sixlin, lpivot[2:3], alpha[2:3,i], a_r, siga_r, b_r, sigb_r

            alpha_pred[blue]=a_b[0]+b_b[0]*d.oned.wavelength[blue]
            alpha_pred[red]=a_r[0]+b_r[0]*d.oned.wavelength[red]
         end
         'quadratic': begin
            ;; Over the whole range
            pf = poly_fit(lpivot,  alpha[*, i], 2)
            alpha_pred = poly(d.oned.wavelength, pf)
         end
         'spline': begin
            spl = spl_init(lpivot, alpha[*, i])
            alpha_pred = spl_interp(lpivot, alpha[*, i], spl, d.oned.wavelength)
         end
         'shift_linear': begin
            ;; These are tricky cases. Here we use the normalised
            ;; spectrum. Then get the fluxes at the lpivot and
            ;; calculate the differences relative to f and then use
            ;; this to shift the flux
            diff_tmp = fltarr(4)
            for jj=0L, 4-1 do begin
               ii_norm_this = where(abs(d.oned.wavelength-lpivot[jj]) lt norm_range and $
                                    finite(d.oned.flux) eq 1 and d.oned.flux ne 0.0, n_norm_this)
               if (n_norm_this eq 0) then stop
               diff_tmp[jj] = median(f[jj]-fnu_orig_scaled[ii_norm_this])
            endfor
               
            sixlin, lpivot[0:1], diff_tmp[0:1], a_b, siga_b, b_b, sigb_b
            sixlin, lpivot[2:3], diff_tmp[2:3], a_r, siga_r, b_r, sigb_r
            delta_pred = alpha_pred*0.0
            delta_pred[blue] = a_b[0]+b_b[0]*d.oned.wavelength[blue]
            delta_pred[red] = a_r[0]+b_r[0]*d.oned.wavelength[red]
            shift = 1
         end
         'none': begin
            ;; Do not have an idea what to do 
            alpha_pred[*] = 0.0
         end
         else: begin
            print, 'Unknown scaling type!'
            stop
         end
      endcase

      dum = min(abs(d.oned.wavelength-lpivot[i_pivot]), i_norm)
      l_norm = d.oned.wavelength[i_norm]

      ii_norm = where(abs(d.oned.wavelength-l_norm) lt norm_range and $
                      finite(d.oned.flux) eq 1 and d.oned.flux ne 0.0, n_norm)
      if (n_norm eq 0) then stop

      if shift then begin
         fnu_total_pred = fnu_orig_scaled+delta_pred
         dfnu_total_pred = dfnu_orig_scaled
      endif else begin
         fnu_total_pred = (1+alpha_pred)*fnu_orig
         dfnu_total_pred = (1+alpha_pred)*dfnu_orig
      endelse

;      stop
      ;; Normalise to one filter.
      fnu_total_pred = fnu_total_pred*f[i_pivot]/median(fnu_total_pred[ii_norm])
      dfnu_total_pred = dfnu_total_pred*f[i_pivot]/median(fnu_total_pred[ii_norm])
      
      stop

;      if objid eq 4590 then stop
;      if check then stop
      
      ;; Save this.
      if keyword_set(use_idlphot) then infix = '-idl-' $
      else if keyword_set(use_sextractor) then infix = '-sex-' $
      else infix = ''

      outfile = ROOT+'Spectra/renorm-'+infix+string(format='(I5.5,"_",I5.5)', propid, objid)+'.fits'
      tmp = {wavelength: 0.0d0, flux_nu: 0.0d0, dflux_nu: 0.0d0, $
             flux_nu_error: 0.0d0, inst_res: 0.0d0, flux: 0.0d0, dflux: 0.0}
      sp = replicate(tmp, n_elements(d.oned.wavelength))
      sp.wavelength = d.oned.wavelength
      sp.flux_nu = fnu_total_pred
      sp.dflux_nu = dfnu_total_pred
      sp.flux_nu_error = dfnu_total_pred

      ;; Then convert to Flambda
      Flambda = fnu_total_pred*2.99792456e-5/(sp.wavelength*1e4)^2
      dFlambda = dfnu_total_pred*2.99792456e-5/(sp.wavelength*1e4)^2
      sp.flux = Flambda
      sp.dflux = dFlambda
      sp.inst_res = d.oned.inst_res
      mwrfits, sp, outfile, /create

;      if objid eq 1917 then stop
      

      old = 0

      if keyword_set(use_idlphot) or keyword_set(use_sextractor) then $
       yrange = minmax([fnu_total_pred, fnu_orig_scaled])
      
      if keyword_set(plot) and not old then begin
         symbols, 30, 1

         pos = get_position_arr(0, nx=1, ny=2, ygap=0, tickf=tf)
         plot, d.oned.wavelength, alpha_pred, position=pos, xtickformat=tf[0], $
               ytickformat=tf[1], ytitle='Alpha', yrange=yrange_alpha, $
               title=t[i].specid
         oplot, lpivot, alpha[*, i], psym=8, color=jb_colour('red')

         pos = get_position_arr(1, tickf=tf)
         plot, d.oned.wavelength, fnu_orig_scaled, position=pos, xtickformat=tf[0], $
               ytickformat=tf[1], ytitle='Flux', /nodata, /noerase, $
               yrange=yrange, /ys
         oplot, d.oned.wavelength, fnu_orig_scaled, color=jb_colour('gray50')
         oplot, d.oned.wavelength, fnu_total_pred
         
         oplot, lpivot, f, psym=8, color=jb_colour('red')
         myoploterr2, lpivot, f, lpivot, lpivot, f-df, f+df, psym=8, color=jb_colour('red')
;         stop
      endif

      
      
      ;;
      ;; Before doing that, however, I want to do a automatic plot
      ;; first to check
      if keyword_set(plot) and old then begin

         if (n_ok_nc gt 1) then $
          yrange_alpha = minmax(alpha[ok_nc, i]) $
         else $
          yrange_alpha = [alpha[ok_nc, i]/2, alpha[ok_nc, i]*2]
          
         if n_ok_nc eq 1 then begin
            ;;
            ;; In this case it is easy - we just use the scale
            ;; calculated above.
            ;;
            fnu_total_pred = fnu_orig*f[i_pivot]/fnu_orig_norm
            alpha_pred[*] = (alpha[ok_nc, i])[0]

         endif else if n_ok_nc eq 4 then begin
            ;;
            ;; This is also straightforward
            ;;
            sixlin, lpivot[0:1], alpha[0:1,i], a_b, siga_b, b_b, sigb_b
            sixlin, lpivot[2:3], alpha[2:3,i], a_r, siga_r, b_r, sigb_r

            alpha_pred[blue]=a_b[0]+b_b[0]*d.oned.wavelength[blue]
            alpha_pred[red]=a_r[0]+b_r[0]*d.oned.wavelength[red]
            
            fnu_total_pred = (1+alpha_pred)*fnu_orig

            ;; Normalise to one filter.
            fnu_total_pred = fnu_total_pred*f[i_pivot]/median(fnu_total_pred[ii_norm])
            
         endif else stop

         symbols, 30, 1

         pos = get_position_arr(0, nx=1, ny=2, ygap=0, tickf=tf)
         plot, d.oned.wavelength, alpha_pred, position=pos, xtickformat=tf[0], $
               ytickformat=tf[1], ytitle='Alpha', yrange=yrange_alpha, $
               title=t[i].specid
         oplot, lpivot, alpha[*, i], psym=8, color=jb_colour('red')

         pos = get_position_arr(1, tickf=tf)
         plot, d.oned.wavelength, fnu_orig_scaled, position=pos, xtickformat=tf[0], $
               ytickformat=tf[1], ytitle='Flux', /nodata, /noerase
         oplot, d.oned.wavelength, fnu_orig_scaled, color=jb_colour('gray50')
         oplot, d.oned.wavelength, fnu_total_pred
         
         oplot, lpivot, f, psym=8, color=jb_colour('red')
         myoploterr2, lpivot, f, lpivot, lpivot, f-df, f+df, psym=8, color=jb_colour('red')

;         if i eq 9 then stop
      endif

      
      

      
      
      
      

      
   endfor

   if keyword_set(ps) then begin
      ps_off
      my_cleanplot, /silent
   endif
   
;   stop
   

end
