
pro find_useful_lines, z, objname, lambda_low=lambda_low, lambda_high=lambda_high, pf=pf

   if (n_elements(lambda_low) eq 0) then lambda_low = 1.8e4   ; a bit conservative
   if (n_elements(lambda_high) eq 0) then lambda_high = 5.1e4 ; a bit conservative
   

   readcol, getenv('HOME')+'/IDL/platefit/etc/linelist-jwst-full.txt', lname, $
            lc, ll, ul, type, format='A,F,F,F,A', /silent
   for i=0L, n_elements(z)-1 do begin
      l_obs = lc*(1+z[i])

      keep = where(l_obs ge lambda_low and l_obs le lambda_high, n_keep)

      print, ''
      print, ''
      print, ''
      print, format='(A,"   Redshift=",F6.4,"    N(lines)=",I0)', objname[i], z[i], n_keep
      print, '-------------------------------------------------------'

      for i_k=0L, n_keep-1 do begin

         j = keep[i_k]
         if n_elements(pf) gt 0 then begin
            if tag_exist(pf, lname[j]+'_FLUX') then begin
               x = struct_var(pf, lname[j]+'_FLUX')
               dx = struct_var(pf, lname[j]+'_FLUX_ERR')
               sn = x[i]/dx[i]
               print, format='(A10,2X,F7.1,2X,F6.1)', lname[j], lc[j], sn
            endif               ; Otherwise silently skip
         endif else $
          print, format='(A10,2X,F7.1)', lname[j], lc[j]
      endfor
      

   endfor

end

function _subset, pf, name
   ;;
   ;; Subset a table for that source
   ;;

   if size(pf, /tname) ne 'STRUCT' then begin
      return, -1
   endif else begin
      i_match = where(strtrim(pf.specid) eq name, n_match)
      if (n_match eq 0) then begin
         ;; Check for renorm match
         i_match = where(strtrim(pf.specid) eq name+' (renorm)', n_match)
         if n_match eq 0 then  $
          return, -1 
      endif
   endelse

   return, pf[i_match[0]]
end

function _get_val, pf, name
   ;;
   ;; Get a value and return meaningful results also if the table does
   ;; not exist.
   ;;
   if size(pf, /tname) ne 'STRUCT' then $
    x = -9999.9 $
   else begin
      x = struct_var(pf, name, missing=-9999.0, /silent)
   endelse
   return, x
end

pro _overview_for_all
   ROOT = '/data2/jarle/JWST/ERO/SMACS/'
   t = mrdfits(ROOT+'redshift-for-sources.fits', 1)

   for i=0L, n_elements(t.object)-1 do begin
      if t[i].confidence lt 2 then continue
      tmp = strsplit(t[i].object, '_', /extract)
      srcid = tmp[1]
;      stop
      jwst_1st_overview_for_source, srcid, /luminosity
      jwst_1st_overview_for_source, srcid, /norm
      jwst_1st_overview_for_source, srcid
   endfor

end


pro jwst_1st_overview_for_source, sourceid, norm=norm, luminosity=luminosity
   ;;
   ;; For the given source, summarise the fluxes in different lines
   ;; for the different platefit runs I made.
   ;;

   ROOT = '/data2/jarle/JWST/ERO/SMACS/'
   PFDIR = ROOT+'Spectra/Platefit/'
   OUTDIR = PFDIR+'FluxOverviews/'
   if not file_test(OUTDIR) then $
    mkdir, OUTDIR

   
   ;;
   ;; The different pfit outputs
   ;;
   suffix = ''
   if size(sourceid, /tname) eq 'STRING' then begin
      ;; Check if the last is a 'b'
      if strmid(sourceid, 0, /reverse) eq 'b' then $
       suffix = '-second'
   endif
   bit = string(format='(I5.5)', sourceid)


;   stop
   individual_fixed = PFDIR+'fit-results-02736_'+bit+'-fixedlinelist'+suffix+'-line.fits'
   individual = PFDIR+'fit-results-02736_'+bit+suffix+'.fits'

   renorm_individual_fixed = PFDIR+'renorm-fit-results-02736_'+bit+'-fixedlinelist'+suffix+'-line.fits'
   renorm_individual = PFDIR+'renorm-fit-results-02736_'+bit+suffix+'.fits'

   full = PFDIR+'platefit-results-lines-fixedlinelist.fits'
   renorm_full = PFDIR+'platefit-results-lines-fixedlinelist-renorm.fits'
   renorm_idl_full = PFDIR+'platefit-results-lines-fixedlinelist-idlrenorm.fits'
   renorm_sex_full = PFDIR+'platefit-results-lines-fixedlinelist-sexrenorm.fits'


   ;;
   ;; Read in the different fits - subsetting when the file has
   ;; multiple objects.
   ;;

   pfs = LIST()
   files = [individual, individual_fixed, renorm_individual, $
            renorm_individual_fixed, full, renorm_full, renorm_idl_full, $
            renorm_sex_full]
;   stop
   for i=0L, n_elements(files)-1 do begin
      if not file_test(files[i]) then begin
         pf = -1
      endif else begin
         pf = mrdfits(files[i], 1, /silent)
;         stop
         if (n_elements(pf) gt 1) then $
          pf = _subset(pf, sourceid)
      endelse
      pfs.Add, pf
   endfor
   

   labels = ['Individual', 'Individual fixed', 'Renorm individua', $
             'Renorm individual fixed', $
             'Full', 'Full renorm', 'Full IDL renorm', 'Full SExtractor renorm']

   labels_short = ['ind', 'indfix', 'reind', 'reindfix', 'full', 'refull', 'refull_idl', $
                   'refull_sex']
   
;   stop

   
   z = pfs[1].z
   
   ;;
   ;; From the redshift, decide on the linelist to use.
   ;;
   
   if (z lt 3.0) then begin
      readcol, ROOT+'list-of-lines-zlt3.txt', lname, lc, sn_once, format='A,F,F', /silent
      norm_line = 'SIII_9533'
   end else begin
      readcol, ROOT+'list-of-lines-zgt3.txt', lname, lc, sn_once, format='A,F,F', /silent
      norm_line = 'H_BETA'
   endelse
   i_norm_line = where(strtrim(lname) eq norm_line)


   n_lines = n_elements(lname)
   print, '#   Line                Individual       |         Individual fixed    |       Individual renorm     |      Individ. fixed renorm  |                 Full        |          Full renorm        |        Full renorm IDL      |    Full renorm SExtractor '


   N_runs = n_elements(pfs)
   
   values = fltarr(n_lines, N_runs)-9999.
   dvalues = fltarr(n_lines, N_runs)-9999.


   norm_flux = dblarr(N_runs)
   missing = -9999.             ; The default
   if keyword_set(luminosity) then begin
      ;;
      ;; We calculate a conversion factor here for a default cosmology.
      ;;
      conv = l_from_flux(1.0d0, z, omega=0.3, lambda=0.7, h=0.7)

      norm_flux[*] = 1d17/conv
      missing = -1d30
      
   endif else if keyword_set(norm) then begin
      ;;
      ;; First, let us get the normalising factors
      ;;
      for j=0L, N_runs-1 do begin
         norm_flux[j] = _get_val(pfs[j], norm_line+'_FLUX')
      endfor      
      
   endif else norm_flux[*] = 1.0

   for i=0L, n_lines-1 do begin

      ;; Assemble the line bit by bit
      line = string(format='(A10)', lname[i])

      for j=0L, N_runs-1 do begin
         x = _get_val(pfs[j], lname[i]+'_FLUX')
         dx = _get_val(pfs[j], lname[i]+'_FLUX_ERR')

         if x eq 0.0 or dx lt 0 then begin
            values[i, j] = missing
            dvalues[i, j] = missing
         endif else begin
            values[i, j] = x/norm_flux[j]
            dvalues[i, j] = dx/norm_flux[j]
         endelse
         
         ;;        line = line+string(format='(" |  ",F10.2,2X,F10.3)', x, dx)
         line = line+string(format='(" |  ",E12.4,2X,E12.4)', values[i, j], dvalues[i, j])
      endfor
      print, line
      
   endfor
   
   ;; For saving we want a transposed format.
   res = {labels: labels_short, full_labels: labels}
   for i=0L, n_lines-1 do begin
      res = create_struct(res, lname[i], reform(values[i, *]))
      res = create_struct(res, 'd'+lname[i], reform(dvalues[i, *]))
   endfor
   res = arr_of_struct_to_struct_of_arr(res)
   if keyword_set(luminosity) then $
    outfile = OUTDIR+'luminosities-'+bit $
   else $
    outfile = OUTDIR+'fluxes-'+bit
   if keyword_set(norm) then outfile = outfile+'-norm'
   outfile = outfile+suffix+'.fits'
   mwrfits, res, outfile, /create
   
   
   

   ;; Next I want to show ratios. 
   for i=0L, N_runs-1 do begin
      x = double(values[*, i])
      x_norm = double(values[*, 5])

      ratio = x/x_norm
      ok = where(finite(ratio) eq 1 and x gt -9000 and x_norm ne 0.0 and x_norm gt -9000, n_ok)
      if (n_ok eq 0) then begin
         print, 'BAD BAD!'
         continue
      endif
      scale = median(ratio[ok])

      missing = where(dx lt 0, n_missing)
      values[*, i] = values[*, i]/scale
      dvalues[*, i] = dvalues[*, i]/scale

      values[missing, i] = missing
      dvalues[missing, i] = missing

      
   endfor


;   plot_norm_diff, lname, values, dvalues
   ;; Do .r jwst_1st_ratio_illustrations
   jwst_1st_show_line_ratio_variations, lname, values, dvalues, i_norm_line, yrange=[-4, 1], $
                                        /ys, labels=labels
   

;   stop
   
end

