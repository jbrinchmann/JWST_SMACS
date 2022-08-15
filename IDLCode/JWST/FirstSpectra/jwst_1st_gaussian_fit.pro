pro jwst_1st_gaussian_fit, objid, pipeline_error=pipeline_error, d_in=d_in, idlrenorm=idlrenorm, $
                           sexrenorm=sexrenorm
   ;;
   ;; Given an object, allow the user to fit Gaussians interactively
   ;; and output results.
   ;;

   if size(objid, /tname) ne 'STRING' then $
    objid = string(format='(I0)', objid)

   target = '2736_'+objid
   propid = 2736
   specid = long(objid)
   bit = string(format='(I5.5,"_",I5.5)', propid, specid)
   suffix_pf = ''
   if strmid(objid, 0, /reverse) eq 'b' then begin
      bit_pf = bit
      bit = bit+'b'
      suffix_pf = '-second'
   endif else bit_pf = bit
   
   jwst_1st_dirs, dd
   OUTDIR = dd.SPECDIR+'GaussianFits/'
   
   t = mrdfits(dd.ROOT+'redshift-for-sources.fits', 1)
   

   i_match = where(t.object eq bit, n_match)

   if (n_match eq 0) then begin
      print, 'No match for '+objid
      stop
      return
   endif
   z = t[i_match[0]].redshift

;   if specid eq 5144 then stop
   ;z = 6.37884





   
   d = jwst_1st_get_coadded(specid)
   ;; We get the de-reddened and continuum fit data from the platefit
   ;; run.

   if n_elements(d_in) eq 0 then begin
      if keyword_set(sexrenorm) then prefix = 'sexrenorm-' $
      else if keyword_set(idlrenorm) then prefix = 'idlrenorm-' $
      else prefix = ''
      fname = prefix+"fit-spec-"+bit_pf+ $
              '-fixedlinelist.fits'
      fsp = mrdfits(dd.SPECDIR+'/Platefit/'+fname, 1)
      
      sigma = empirical_sigma(fsp.wavelength, fsp.flux)
      
      wave = fsp.wavelength
      flux = fsp.flux-(fsp.continuum+fsp.resid_cont)

      if not keyword_set(pipeline_error) then $
       dflux = sigma $
      else $
       dflux =  fsp.dflux
      
      R = interpol(d.oned.R, d.oned.wavelength, wave/1e4)
   endif else begin
      wave = d_in.oned.wavelength
      flux = d_in.oned.flux
      dflux = d_in.oned.dflux
      R = d_in.oned.R
   endelse

;;   stop
   
   jwst_1st_fit_gaussian_interactive, wave, flux, dflux, z, result=result, LSF_R=R
   
   ;;
   ;; Store results as a text file.
   ;;

   if keyword_set(pipeline_error) then suffix = '-pipe_err' else suffix = ''
   outfile = OUTDIR+'gaussfit-'+bit+suffix+'.txt'

   if file_test(outfile) then begin
      print, 'The output file '+outfile
      print, 'Do you want to move this out of the way?'
      stop
   endif
   openw, lun, outfile, /get_lun
   printf, lun, '# Line Flux dFlux Lcenter dLcenter Sigma dSigma Level dLevel'

   f = '(A12,2X,8(E13.4,2X))'
   tn = tag_names(result)
   for i=0L, n_elements(tn)-1 do begin
      key = strupcase(tn[i])
      if key eq 'LINES' then continue ; skip this
      
      x = struct_var(result, key)

      printf, lun, format=f, key, x.flux, x.dflux, x.lc, x.dlc, x.sigma, x.dsigma, $
              x.level, x.dlevel
   endfor
   free_lun, lun


   
end
