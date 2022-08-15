pro _run_all
   
   todo = ['A', 'B', 'C', 'D']

   ROOT = getenv('HOME')+'/Work/JWST/FirstSpectra/Platefit/'

   
   
   for i=0L, n_elements(todo)-1 do begin
      parfile = ROOT+'pfit_init-'+todo[i]+'.par'
      platefit_init_jwst, parfile, /mpa, /no_named_fiber, fiber=fiber, info=info


      sp = jwst_1st_load(todo[i])
      res = jwst_1st_platefit(sp, fiber, fit_spec=fit_spec)
      
      outfile = ROOT+'fit-results-'+todo[i]+'.fits'
      fitspec = ROOT+'fit-spec-'+todo[i]+'.h5'
      
      mwrfits, res, outfile, /create
      struct_to_hdf5, fit_spec, fitspec
      
   endfor

   

end


function jwst_1st_platefit, sp, fiber, fit_spec=fit_spec
   ;;
   ;; Run platefit on the JWST first release spectra. 
   ;;

   common burstinfo, burst_lib, burst_wl, wavelength_range_for_fit, zval, $
    ztags, ssp_ages, burst_model_file, ssp_norms, model_dispersion
   
   common lfit_settings, max_sigma, max_shift, max_iter, lines_in_vacuum, use_c_code


   if (n_elements(fiber) eq 0) then begin
      print, 'For JWST is usually best to create an unnamed fibre and provide this here.'
      stop
      fiber = {fiberspec_jwst}
      if (n_elements(z) eq 0) then begin
         print, 'You need to provide the redshift if you do not give a fibre!'
         return, -1
      endif
      fiber.z = z
;      fiber.plateid = i
   endif

   fiber.z = sp.z
   
   z = fiber.z
   mnmx = minmax(sp.wave)
   wavelength_range_for_fit = [mnmx[0]/(1+z[0]), mnmx[1]/(1+z[0])]



   ;;------------------------------------------------------------------
   ;; Interpolate onto uniform log lambda scale
   ;;
   ;; Note that this assumes that the input cube is in vacuum wavelengths
   ;;------------------------------------------------------------------
   vacwl_orig = sp.wave
   
   tmp_logwl = alog10(vacwl_orig)
   dw=min(abs(tmp_logwl-shift(tmp_logwl,1)), mi)
   if (dw eq 0.0) then begin
      print, 'The wavelengths are identical somewhere!'
      stop
   endif
   if dw lt 1e-4 then dw = 1e-4                      ; Seems a logical minimum 
   xnew = mkarr(min(tmp_logwl), max(tmp_logwl), dw ) ; vacuum wavelengths 

   ;; Get the fluxes
   flux_orig = sp.flux
   err_orig = sp.dflux


   ;;
   ;; Insert S/N information
   ;;
   fiber.sn_median = median(flux_orig/err_orig)
   
   use = where(err_orig gt 0 and err_orig lt 1e5 and finite(flux_orig) eq 1, $
               n_use)
   if (n_use eq 0) then return, fiber
   ynew = interpol(flux_orig[use], vacwl_orig[use], 10.0^xnew)
   errnew = interpol(err_orig[use], vacwl_orig[use], 10.0^xnew)

   logwl = xnew                 ; vacuum
   flux = ynew
   err = errnew
   

   ;;
   ;; The stellar spectra are in air so we also need an air wavelength axis.
   ;;
   airwl = 10.0^xnew            ; still vacuum
   vactoair, airwl              ; After this it is in air.

   restwl_air = airwl/(1+fiber.z)
   restwl = 10.0^logwl/(1+fiber.z)

   ;; Added a check for finite flux.
   ok = where(err gt 0 and finite(err^2) eq 1 $
              and err lt 1e10 and $
              finite(flux)  eq 1)

   ;;
   ;; De-redden and adjust for 1+z scaling.
   ;;
   ;; Use O'donnell extinciton curve to get e^(Tau) -- this is what SFD used!
   A_v = 3.1 * fiber.e_bv_sfd
   a_odonnell = ext_odonnell(airwl, 3.1)

   e_tau = exp(a_odonnell * A_v / 1.086)
   
   ;; Scale by galactic reddening.
   flux = flux*e_tau
   err = err*e_tau

   ;; Scale by (1+z)
   flux = flux*(1+fiber.z)
   err = err*(1+fiber.z)

   ;;------------------------------------------------------------------
   ;; LSF correction
   ;; 
   ;; By default we will apply a correction for the LSF to deconvolve
   ;; from the emission line widths. This is done using the inst_res
   ;; option to platefit. This needs a bit of care because platefit
   ;; assumes this follows the SDSS convention which has inst_res in
   ;; log lambda units (actually in pixels but this conversion is done
   ;; by platefit (in plate_lfit.pro) when reading an SDSS type
   ;; spectrum file.
   ;;
   ;; I will need to get the FWHM in AA and then convert this into
   ;; sigmas so divide this by sqrt(8*ln(2)). And then I need to
   ;; convert this to delta log lambda.
   ;;
   ;;    Delta log l = log l_1/logl_0
   ;;
   ;; Thus from Delta l = l_1-l_0 = l_0 (l1/l0-1), I get
   ;;
   ;;    Delta log l = Delta l/(ln(10) * l_0)
   ;;
   ;; to first order
   ;;
   ;; However at the moment I do NOT have this function so use 0.
   ;;
   ;;------------------------------------------------------------------

   print, "Current no NIRSPEC LSF implemented"
   inst_res = flux*0.0

      ;;-----------------------------------------------------------
   ;; Fit continuum using NNLS
   ;; FIBER will contain the parameters required to reconstruct the
   ;; continuum, but the continuum is also returned in the
   ;; continuum variable.
   ;;-----------------------------------------------------------

   if keyword_set(subtracted) then begin
      ;; Ok, this spectrum has been continuum subtracted already.
      continuum = flux[ok]*0+median(flux) ; smooth(flux[ok], 100, /nan)
      fiber.best_model_z = 0.0
   endif else begin
      extreme = 0
      fiber=fiber_continfit(fiber, logwl[ok], restwl[ok], $
                            flux[ok], err[ok], $
                            cont=continuum, title=title, $
                            /adjust_plotrange , /noplot, extreme=extreme, $
                           status=status)
      

      if (fiber.best_model_z eq 0.0 or status lt 0) then begin
;         stop
         print, 'Continuum fitting failed so I will use a simple continuum'
         continuum = flux[ok]*0+median(flux) ; smooth(flux[ok], 100, /nan)
      endif
   endelse
   

      
   ;;------------------------------------------------------------
   ;; Calculate indices from data
   ;; Since at this point we might have emission lines, this has
   ;; substantial uncertainties.
   ;;------------------------------------------------------------
   fiber_index, restwl[ok], flux[ok], err[ok], fiber = fiber

   ;;---------------------------------------------------------------
   ;; Calculate indices from model fit - this is only necessary if
   ;; it is needed
   ;;---------------------------------------------------------------
   fiber_index, restwl[ok], continuum, $
                fltarr(n_elements(restwl[ok])) + 1, $
                fiber = fiber, auxtag = '_MODEL'
         
   ;;------------------------------------------
   ;; EMISSION LINE FITTING
   ;;------------------------------------------

   ;;----------------------------------------
   ;; Here we need the inverse variance
   ;;----------------------------------------
   invvar = 1.0/err^2
   cont1 = continuum
   
   ;;---------------------------------------------------------------
   ;; Fit emission lines.
   ;; This routine also calculates the equivalent width of lines,
   ;; plus continuum indices on the emission line subtracted
   ;; spectra.
   ;;
   ;; The nonsdss keyword means that the FIBER_DERED routine is not
   ;; called - in other words, the spectra are expected to be
   ;; corrected for galactic reddening etc.
   ;;---------------------------------------------------------------
   
   ;; The tying of lines together is poorly implemented and needs some
   ;; thought for the re-implementation.

   if keyword_set(split_voff) then begin
      ;; sv decides whether we tie velocity offsets together or not
      sv = 0
   endif else begin
      if (fiber.z lt 0.9) then sv = 'Balmer' $
      else sv = 'Forbidden'     ; Note that He & Ly are not Balmer!
   endelse

;   stop

;  fiber.z = fiber.z+5e-4
   fiber = fiber_lfit(fiber, logwl[ok], flux[ok], invvar[ok], $
                      continuum=continuum, /nonsdss, $
                      title=title, nebular=nebular, $
                      specfit=specfit, noplot=noplot, $
                      stellar=stellar, /adjust_plotrange, $
                      resid_cont=resid_cont, /use_input_cont, $
                      same_voff_sigma=sv, inst_res=inst_res[ok], $
                      do_not_smooth=do_not_smooth, doublepeak=doublepeak, $
                      doubleinfo=doubleinfo)

   if (keyword_set(doublepeak)) then begin
      nebular1 = doubleinfo.neb1
      nebular2 = doubleinfo.neb2
   endif

   ;; Copy back the fitted continuum.
   continuum = cont1 
;   stop

   ;; Finally, if the FIT_SPEC argument is present, restore the
   ;; fitted spectrum on the wavelength range of the actual spectrum. 
   if arg_present(fit_spec) then begin

      ;; This rest wavelength is needed for input to
      ;; restore_spectrum_fit, but it is not the "correct" restwl

      this_restwl =  vacwl_orig/(1+fiber.z) ; sp.wavelength/(1+fiber.z)
      v = restore_spectrum_fit(fiber, this_restwl, status=status, /wavelength_in_vac)

      ;;
      ;; Now incorporate our velocity offset from above to get the
      ;; best-estimate of a rest-wavelength - in vacuum
      ;;
      z_corr= (1+fiber.z)*10.0^(fiber.v_off_forbidden/(2.99792e5*alog(10.0)))-1
      restwl_vac = vacwl_orig/(1+z_corr)
      
      if (status lt 0) then $
       v.continuum[*] = median(flux[ok])
      ;; We also want the Galactic foreground correction and the
      ;; residual continuum.
      this_e_tau = exp(ext_odonnell(sp.wave, 3.1)*A_V/1.086)
      this_rc = interpol(resid_cont, 10.0^logwl[ok], sp.wave)

      if keyword_set(keep_dereddened) then begin
         ;; In this case we will return a fit_spec structure that
         ;; corresponds to the dereddened spectrum. 
         corr_spec = this_e_tau*(1+fiber.z)
         corr = 1.0
      endif else begin
         ;; Then adjust the continuum & emission line spectra to match
         ;; the observed spectrum.
         corr = 1.0/(this_e_tau*(1+fiber.z))
         corr_spec = 1.0
      endelse
         
      v.continuum = v.continuum*corr
      v.emlines.nebular_clip = v.emlines.nebular_clip*corr
      v.emlines.nebular = v.emlines.nebular*corr
      ;; AND THE RESIDUAL CONTINUUM!
      this_rc = this_rc*corr
      

      ;; This is the spectrum fit.
      this_specfit = v.continuum + v.emlines.nebular_clip + this_rc

      ;; And this is the "stellar" part
      this_stellar = v.continuum + this_rc
      
      ;; Create the output structure format.
      tmp = {continuum: 0.0, nebular_clip: 0.0, $
             nebular_all: 0.0, specfit: 0.0, stellar: 0.0, $
             e_tau: 0.0, resid_cont: 0.0, restwl: 0.0, wavelength: 0.0, $
             flux: 0.0, dflux: 0.0, air_wavelength: 0.0, restwl_noshift: 0.0}

      if keyword_set(complex) then begin
         tmp = create_struct(tmp, 'nebular_complex_only', 0.0, $
                             'nebular_complex_clip', 0.0, $
                             'nebular_complex_all', 0.0, $
                             'specfit_complex', 0.0)

      endif


      
;      stopb
      
      fit_spec = replicate(tmp, n_elements(sp.wave))
      fit_spec.continuum = v.continuum
      fit_spec.nebular_clip = v.emlines.nebular_clip
      fit_spec.nebular_all = v.emlines.nebular
      fit_spec.specfit = this_specfit
      fit_spec.stellar = this_stellar
      fit_spec.e_tau = this_e_tau
      fit_spec.resid_cont = this_rc
      fit_spec.restwl = restwl_vac ; In vacuum, corrected for velocity shift!
      fit_spec.restwl_noshift = this_restwl ; In vacuum, not corrected for velocity shift!
      fit_spec.air_wavelength = sp.wave   ; in air
      fit_spec.wavelength = vacwl_orig

      ;; Interpolate onto the vacuum wavelengths so everything refers
      ;; to vacuum. 
      fit_spec.flux = interpol(flux_orig, sp.wave, vacwl_orig)*corr_spec
      fit_spec.dflux = interpol(err_orig, sp.wave, vacwl_orig)*corr_spec; err_orig

   endif

   return, fiber
end
