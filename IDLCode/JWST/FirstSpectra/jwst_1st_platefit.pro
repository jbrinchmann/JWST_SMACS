pro _split_array, fibers, f_tmp, indices
   
   f_tmp = arr_of_struct_to_struct_of_arr(fibers)
   
   tags = tag_names(f_tmp)
   
   first = 1
   for i=0L, n_elements(tags)-1 do begin
      txt = tags[i]
      if stregex(txt, '^LICK.*', /bool) or $
       stregex(txt, '^DTT.*', /bool) or $
       stregex(txt, '^BH', /bool) then begin 
         
         ;; Ok ,this is to move so first copy
         x = struct_var(f_tmp, txt)
         if first then begin
            indices = create_struct(txt, x)
            first = 0
         endif else $
          indices = create_struct(indices, txt, x)
         
         ;; And then remove
         f_tmp = struct_remove_keys(f_tmp, keys=txt)
         
      endif
      
   endfor

;   stop
   if (n_elements(f_tmp.z) gt 1) then begin
      f_tmp = struct_of_arr_to_arr_of_struct(f_tmp)
      indices = struct_of_arr_to_arr_of_struct(indices)
   endif
   
end


pro _generate_par_files, fixed_linelist=fixed_linelist
   ROOT = '/data2/jarle/JWST/ERO/SMACS/'

   SPECDIR = ROOT+'Spectra/'
   PFDIR = SPECDIR+'Platefit/'

   ;; Get the table with redshifts
   t = mrdfits(ROOT+'redshift-for-sources.fits', 1)
   todo = where(t.confidence gt 1, n_todo)

   if keyword_set(fixed_linelist) then begin
      linelist_file = 'linelist-fixed.txt'
      suffix = '-fixedlinelist'
   endif else suffix = ''
   

   for itmp=0L, n_todo-1 do begin

      i = todo[itmp]
      specid = strmid(t[i].object, 5)

      openw, lun, PFDIR+'pfit_init-'+specid+suffix+'.par', /get_lun
      printf, lun, '# Parameter file for platefit_init'
      printf, lun, '# Automatically generated by _generate_par_files in jwst_1st_platefit.pro'
      printf, lun, '#'
      printf, lun, "pipelinedir = getenv('HOME')+'/IDL/platefit/'"
      printf, lun, "twoddir = '/data1/jarle/SDSS_Spectra/DR7/'"
      printf, lun, "oneddir = '/data1/jarle/SDSS_Spectra/DR7/'"
      printf, lun, "dustdir = getenv('HOME')+'/IDL/SFD_Dust_Maps/maps/'"
      printf, lun, "platefitdir = '"+PFDIR+"'"
      printf, lun, "# Use a model file with good resolution way into the UV and IR coverage"
      printf, lun, "burst_model_file  = 'bc_models_subset_cb08_milesx__bursts_extlam_jwst.fit'"
      printf, lun, '# The wavelength range to set is actually set by jwst_1st_platefit_one so no need to set here.'
      printf, lun, "wavelength_range_for_fit = [3600, 9300]"
      printf, lun, "use_z = [0.004, 0.008, 0.017, 0.04]"
      printf, lun, "fgdust_cor = 1"
      printf, lun, "use_fortran_indices = 0"
      printf, lun, ""
      printf, lun, "# Use this linelist."
      if keyword_set(fixed_linelist) then $
       printf, lun, "linelist = platefitdir+'"+linelist_file+"'" $
      else $
       printf, lun, "linelist = platefitdir+'linelist-"+specid+".txt'"
      printf, lun, "# lines_in_vacuum = 1"
      printf, lun, "indexlist = pipelinedir+'etc/indexlist.txt'"
      free_lun, lun
   endfor
   
end


pro _run_individual_observations, ps=ps, pipeline_error=pipeline_error
   ;;
   ;; This runs platefit for each separate observation instead of
   ;; adding them together. 
   ;;

   ROOT = '/data2/jarle/JWST/ERO/SMACS/'

   SPECDIR = ROOT+'Spectra/'
   PFDIR = SPECDIR+'Platefit/'

   t = mrdfits(ROOT+'redshift-for-sources.fits', 1)

   todo = where(t.confidence gt 1, n_todo)

   print, 'I will do a total of ', n_todo, ' galaxies.'

   ;; For this I will NOT do any renormalisation - I want to look at
   ;; the difference between re-observations.


   ;; I also use a fixed line list but a smaller one than 
   parfile = PFDIR+'pfit_init-fixed-strong.par'
   platefit_init_jwst, parfile, /mpa, /no_named_fiber, fiber=fiber, info=info
   
   fibers07 = replicate(fiber, n_elements(todo))
   fibers08 = replicate(fiber, n_elements(todo))
   

   jwst_1st_dirs, dd
   
   if keyword_set(pipeline_error) then suffix = '-pipeline-err' else suffix = ''
   if keyword_set(ps) then begin
      Outfile = dd.figdir+'platefit_log-individual-obs'+suffix+'.ps'
      ps_on, outfile, /landscape, /color, /times
      !p.font = 0
      !x.thick = 2
      !y.thick = 2
   endif
      
   for itmp=0L, n_elements(todo)-1 do begin
      ;;
      ;; Get name
      ;;
      i = todo[itmp]


      tmp = strsplit(t[i].object, '_', /extract)
      propid = long(tmp[0])
      objid = long(tmp[1])
      specid = strmid(t[i].object, 5)


      r = jwst_1st_get_each_obs(objid)

      for iobs=0, 1 do begin
         if iobs eq 0 then begin
            d = r.obs07['combined']
            fiber = fibers07[itmp]
         endif else begin
            d = r.obs08['combined']
            fiber = fibers08[itmp]
         endelse

         ;; Insert informational data
         fiber.specid = specid
         fiber.z = t[i].redshift

         if not keyword_set(pipeline_error) then begin
            ;; By default replace the pipeline uncertainty with an
            ;; empirical one.
            df = empirical_sigma(d.oned.wavelength, d.oned.flux)
            d.oned.dflux = df
         endif

         scale = 1e-8
         sp = d.oned
         sp.wavelength = sp.wavelength*1e4 ; Convert to AA
         sp.flux = sp.flux*scale
         sp.dflux = sp.dflux*scale
      
         res = jwst_1st_platefit(sp, fiber, /show_fit)

         if iobs eq 0 then $
          fibers07[itmp] = res $
         else $
          fibers08[itmp] = res
      endfor
      
   endfor


   mwrfits, fibers07, PFDIR+'individual-obs07'+suffix+'.fit', /create
   mwrfits, fibers08, PFDIR+'individual-obs08'+suffix+'.fit', /create


   if keyword_set(ps) then begin
      ps_off
      my_cleanplot, /silent
   endif
   
end


pro _run_all, ps=ps, fixed_linelist=fixed_linelist, $
                 renormidl=renormidl, renormsex=renormsex


   ROOT = '/data2/jarle/JWST/ERO/SMACS/'

   SPECDIR = ROOT+'Spectra/'
   PFDIR = SPECDIR+'Platefit/'

   t = mrdfits(ROOT+'redshift-for-sources.fits', 1)

   todo = where(t.confidence gt 1, n_todo)

   print, 'I will do a total of ', n_todo, ' galaxies.'

;   stop

   if keyword_set(renormidl) then $
    renorminfix = '-idl-' $
   else if keyword_set(renormsex) then $
    renorminfix = '-sex-' $
   else renorminfix = ''

   n_renorm = 0
   for itmp=0L, n_todo-1 do begin
      i = todo[itmp]
      tmp = strsplit(t[i].object, '_', /extract)

      if strmid(tmp[1], 0, /reverse) eq 'b' then begin
         tmp[1] = strmid(tmp[1], 0, strlen(tmp[1])-1)
;         stop
      endif
      propid = long(tmp[0])
      objid = long(tmp[1])

      outfile = ROOT+'Spectra/renorm-'+renorminfix+string(format='(I5.5,"_",I5.5)', propid, objid)+'.fits'
      if file_test(outfile) then n_renorm = n_renorm+1
   endfor

   ;; Get the table with redshifts

   
   if keyword_set(fixed_linelist) then begin
      suffix = '-fixedlinelist'
   endif else suffix = ''
   

   OUTDIR = ROOT+'Figures/'


   ;; The behaviour is rather different if we have a fixed_linelist
   ;; because then we create a single fibre array and save this at the
   ;; end. 

   if keyword_set(fixed_linelist) then begin
      parfile = PFDIR+'pfit_init'+suffix+'.par'
      platefit_init_jwst, parfile, /mpa, /no_named_fiber, fiber=fiber, info=info

      fibers = replicate(fiber, n_elements(todo))
      fibers_renorm = replicate(fiber, n_renorm)

      if keyword_set(ps) then begin
         Outfile = OUTDIR+'platefit_log'+suffix+'.ps'
         ps_on, outfile, /landscape, /color, /times
         !p.font = 0
         !x.thick = 2
         !y.thick = 2
      endif
      
      
   endif

   

   i_renorm = 0
   old_suffix = suffix
   for itmp=0L, n_elements(todo)-1 do begin
      ;;
      ;; Get name
      ;;
      i = todo[itmp]

      suffix = old_suffix
      
      object = t[i].object
      original_object = object  ; Object might be modfied below
      tmp = strsplit(object, '_', /extract)
      
      secondary = 0
      if strmid(tmp[1], 0, /reverse) eq 'b' then begin
         tmp[1] = strmid(tmp[1], 0, strlen(tmp[1])-1)
         secondary = 1
         object = string(format='("2736_",I0)', long(tmp[1]))
         old_suffix = suffix
         suffix = suffix+'-second'
;         stop
      endif else begin
         object = string(format='("2736_",I0)', long(tmp[1]))
      endelse

      propid = long(tmp[0])
      objid = long(tmp[1])


      
      specid_full = strmid(original_object, 6)
      specid = strmid(object, 5)
      

;      if objid ne 5144 then continue


      if keyword_set(fixed_linelist) then begin
         fiber = fibers[itmp]
      endif else begin
         parfile = PFDIR+'pfit_init-'+specid_full+suffix+'.par'
         platefit_init_jwst, parfile, /mpa, /no_named_fiber, fiber=fiber, info=info
      endelse

      ;; Insert informational data
      fiber.specid = specid_full
      fiber.z = t[i].redshift

;      stop
      d = jwst_1st_get_coadded(long(specid))

      if t[i].reextract gt 0 then begin
         ;; See first if there is a bad file for this re-extraction
         badfile = SPECDIR+'bad-reextract-'+object+'.h5'
         if (not file_test(badfile)) then begin
            stop_for_bad = 1
         endif else begin
            stop_for_bad = 0
            bad = hdf5_to_struct(badfile)
         endelse
         
         d = jwst_nirspec_reextract(d, t[i].reextract, $
                                    narrow_sum=t[i].width, /ap_background, $
                                    bad=bad)

         ;; This is in jwst_load_nirspec_data (could be factored out)
         df = empirical_sigma(d.oned.wavelength, d.oned.flux)
         d.oned.dflux = df

         if stop_for_bad then stop
         scale = 1e-8
      endif else scale = 1e-4

      infix = string(format='(I5.5,"_",I5.5)', propid, objid)
      if keyword_set(ps) and not keyword_set(fixed_linelist) then $
       psfile = OUTDIR+'platefit_log-'+infix+suffix+'.ps'
      
      sp = d.oned
      sp.wavelength = sp.wavelength*1e4 ; Convert to AA
      sp.flux = sp.flux*scale
      sp.dflux = sp.dflux*scale

      
      if t[i].reextract le 0 and objid eq 4590 then begin
         bad_x = [1785, 1786, 1787, 1978, 1979, 1980, 1981, 1982, $
                  1983, 1984, 1985]
         sp.flux[bad_x] = !values.d_nan
      endif
      
      res = jwst_1st_platefit(sp, fiber, fit_spec=fit_spec, /show_fit, $
                              psfile=psfile)

      outfile = PFDIR+'fit-results-'+infix+suffix+'.fits'
      fitspec_file = PFDIR+'fit-spec-'+infix+suffix+'.fits'


      if keyword_set(fixed_linelist) then begin
         ;;
         ;; In this case it needs to be split
         ;;
         outfile_line = PFDIR+'fit-results-'+infix+suffix+'-line.fits'
         outfile_index = PFDIR+'fit-results-'+infix+suffix+'-idx.fits'
         fitspec_file = PFDIR+'fit-spec-'+infix+suffix+'.fits'


         _split_array, res, f_tmp, indices

         mwrfits, f_tmp, outfile_line, /create
         mwrfits, indices, outfile_index, /create
         mwrfits, fit_spec, fitspec_file, /create
         
      endif else begin
         outfile = PFDIR+'fit-results-'+infix+suffix+'.fits'
         fitspec_file = PFDIR+'fit-spec-'+infix+suffix+'.fits'
         mwrfits, res, outfile, /create
         mwrfits, fit_spec, fitspec_file, /create
      endelse
;;      struct_to_hdf5, fit_spec, fitspec_file

      if keyword_set(fixed_linelist) then $
       fibers[itmp] = fiber

      ;; Then check whether there is a renormalised version of this
      ;; spectrum.
      outfile = ROOT+'Spectra/renorm-'+renorminfix+string(format='(I5.5,"_",I5.5)', propid, objid)+'.fits'
      if file_test(outfile) then begin
         sp = mrdfits(outfile, 1)
         print, 'Reading renorm file'+outfile
         sp.wavelength = sp.wavelength*1e4 ; Convert to AA

         ;; The flux here is in Flambda units but we want to scale a
         ;; bit
         sp.flux = sp.flux*1e17
         ;; I also prefer my empirical uncertainty.
         df = empirical_sigma(sp.wavelength, sp.flux)
         sp.dflux = df

         
         fiber.specid = specid_full+' (renorm)'

         if keyword_set(renormidl) then $
          txtbit = 'idlrenorm' $
         else if keyword_set(renormsex) then $
          txtbit = 'sexrenorm' $
         else $
          txtbit = 'renorm'
         if keyword_set(ps) and not keyword_set(fixed_linelist) then $
          psfile = OUTDIR+'platefit_log-'+infix+'-'+txtbit+'.ps'
         res = jwst_1st_platefit(sp, fiber, fit_spec=fit_spec, /show_fit, $
                                psfile=psfile)
         
         if keyword_set(fixed_linelist) then begin
            ;;
            ;; In this case it needs to be split
            ;;
            outfile_line = PFDIR+txtbit+'-fit-results-'+infix+suffix+'-line.fits'
            outfile_index = PFDIR+txtbit+'-fit-results-'+infix+suffix+'-idx.fits'
            fitspec_file = PFDIR+txtbit+'-fit-spec-'+infix+suffix+'.fits'
            

            _split_array, res, f_tmp, indices

            mwrfits, f_tmp, outfile_line, /create
            mwrfits, indices, outfile_index, /create
            mwrfits, fit_spec, fitspec_file, /create

            fibers_renorm[i_renorm] = res
            i_renorm = i_renorm+1
            
         endif else begin
            outfile = PFDIR+txtbit+'-fit-results-'+infix+suffix+'.fits'
            fitspec_file = PFDIR+txtbit+'-fit-spec-'+infix+suffix+'.h5'
            
            mwrfits, res, outfile, /create
            mwrfits, fit_spec, fitspec_file, /create
         endelse
;         struct_to_hdf5, fit_spec, fitspec
      endif
         


      
   endfor

   if keyword_set(fixed_linelist) then begin
      ;;
      ;; Move some keys into a different structure because there are
      ;; too many keys in it
      ;;
      if keyword_set(ps) then begin
         ps_off
         my_cleanplot, /silent
      endif

      _split_array, fibers, f_tmp, indices
      mwrfits, f_tmp, PFDIR+'platefit-results-lines-fixedlinelist.fits', /create
      mwrfits, indices, PFDIR+'platefit-results-indices-fixedlinelist.fits', /create


      _split_array, fibers_renorm, f_tmp, indices
      mwrfits, f_tmp, PFDIR+'platefit-results-lines-fixedlinelist-'+txtbit+'.fits', /create
      mwrfits, indices, PFDIR+'platefit-results-indices-fixedlinelist-'+txtbit+'.fits', /create

      
;      stop
   endif

   
end


function jwst_1st_platefit, sp, fiber, fit_spec=fit_spec, show_fit=show_fit, $
                            psfile=psfile
   ;;
   ;; Run platefit on the JWST first release spectra. 
   ;;

   common burstinfo, burst_lib, burst_wl, wavelength_range_for_fit, zval, $
    ztags, ssp_ages, burst_model_file, ssp_norms, model_dispersion
   
   common lfit_settings, max_sigma, max_shift, max_iter, lines_in_vacuum, use_c_code
   common fitinfo, linepars, indexpars 


   noplot = 1                   ; Turn off the default plotting
   
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

   z = fiber.z
   mnmx = minmax(sp.wavelength)
   wavelength_range_for_fit = [mnmx[0]/(1+z[0]), mnmx[1]/(1+z[0])]



   ;;------------------------------------------------------------------
   ;; Interpolate onto uniform log lambda scale
   ;;
   ;; Note that this assumes that the input cube is in vacuum wavelengths
   ;;------------------------------------------------------------------
   vacwl_orig = sp.wavelength
   
   tmp_logwl = alog10(vacwl_orig)
   dw=min(abs(tmp_logwl-shift(tmp_logwl,1)), mi)
   if (dw eq 0.0) then begin
      print, 'The wavelengths are identical somewhere!'
      stop
   endif
   if dw lt 1e-4 then dw = 1e-4 ; Seems a logical minimum

   xnew = mkarr(min(tmp_logwl), max(tmp_logwl), dw ) ; vacuum wavelengths 

   ;; Get the fluxes
   flux_orig = sp.flux
   err_orig = sp.dflux


   ;;
   ;; Insert S/N information
   ;;
   fiber.sn_median = median(flux_orig/err_orig)
   
   use = where(err_orig gt 0 and err_orig lt 1e10 and finite(flux_orig) eq 1, $
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

   ;; Use nominal R values given on JWST page
   inst_res = sp.inst_res

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
;      stop
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

   if keyword_set(show_fit) then begin

      if n_elements(psfile) gt 0 then begin
         ps_on, psfile, /landscape, /color, /times
         !p.font = 0
         !x.thick = 2
         !y.thick = 2
      endif

      
      ;;
      ;; Show the fit result
      ;;
      xmin = 0.1
      xmax = 0.9
      ymin = 0.1
      ymax = 0.9
      dy = 0.05


      mnmx = minmax(restwl)
      keep_lines = where(linepars.centr gt mnmx[0] and $
                         linepars.centr lt mnmx[1], n_lines)
      if (n_lines eq 0) then stop

      this_linepars = linepars[keep_lines]

      if (n_lines gt 10) then begin
         yspec = 0.4
         pos = get_position_arr(0, nx=floor(n_lines/3.+0.5), ny=3, $
                                xmin=xmin, xmax=xmax, ymin=yspec+dy, $
                                ymax=ymax, ygap=0.03, xgap=0.03)
      endif  else if (n_lines gt 5) then begin
         yspec = 0.5 
         pos = get_position_arr(0, nx=floor(n_lines/2.+0.5), ny=2, $
                                xmin=xmin, xmax=xmax, ymin=yspec+dy, $
                                ymax=ymax, ygap=0.03, xgap=0.03) 
      endif else begin
         yspec = 0.6
         pos = get_position_arr(0, nx=n_lines, xmin=xmin, xmax=xmax, $
                                ymin=yspec+dy, ymax=ymax, ygap=0.03)
      endelse
      
      pos_spec = [xmin, ymin, xmax, yspec]
      
      plot, restwl, flux, xtitle='Wavelength', ytitle=TeXtoIDL('F_\lambda'), $
            position=pos_spec, charsize=1.5
      oplot, restwl[ok], specfit, color=jb_colour('CornFlowerBlue')
      oplot, restwl[ok], continuum+resid_cont, color=jb_colour('orange')
      jb_text, -0.05, -0.1, String(format='("Z=",F0.4)', fiber.z), $
               align=1, charsize=1.6, /relative, /fraction, $
               color=jb_colour('red')


      for i_line=0L, n_lines-1 do begin
         pos = get_position_arr(i_line)
         info = this_linepars[i_line]
         xrange = info.centr+[-40, 40]
         inside = where(restwl ge xrange[0] and restwl le xrange[1])
         mnmx = minmax(flux[inside])
         plot, restwl, flux, xrange=xrange, $
               xtickformat='noticks', ytickformat='noticks', position=pos, $
               /noerase, yrange=[mnmx[0], 1.15*mnmx[1]], /ys
         oplot, restwl[ok], specfit, color=jb_colour('orange')
         jb_text, 0.5, -0.1, this_linepars[i_line].line, color=jb_colour('CornFlowerBlue'), $
                  /fraction, /relative, align=0.5, charsize=1.6

      endfor

      xyouts, 0.5*(xmax+xmin), ymax+0.04, "SpecID="+fiber.specid, charsize=1.5, /normal

;      stop
      if n_elements(psfile) gt 0 then begin
;         stop
         ps_off
         my_cleanplot, /silent
      endif 
      
   endif


   ;; Copy back the fitted continuum.
   continuum = cont1 

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
      this_e_tau = exp(ext_odonnell(sp.wavelength, 3.1)*A_V/1.086)
      this_rc = interpol(resid_cont, 10.0^logwl[ok], sp.wavelength)

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
      
      fit_spec = replicate(tmp, n_elements(sp.wavelength))
      fit_spec.continuum = v.continuum
      fit_spec.nebular_clip = v.emlines.nebular_clip
      fit_spec.nebular_all = v.emlines.nebular
      fit_spec.specfit = this_specfit
      fit_spec.stellar = this_stellar
      fit_spec.e_tau = this_e_tau
      fit_spec.resid_cont = this_rc
      fit_spec.restwl = restwl_vac ; In vacuum, corrected for velocity shift!
      fit_spec.restwl_noshift = this_restwl ; In vacuum, not corrected for velocity shift!
      fit_spec.air_wavelength = sp.wavelength   ; in air
      fit_spec.wavelength = vacwl_orig

      ;; Interpolate onto the vacuum wavelengths so everything refers
      ;; to vacuum. 
      fit_spec.flux = interpol(flux_orig, sp.wavelength, vacwl_orig)*corr_spec
      fit_spec.dflux = interpol(err_orig, sp.wavelength, vacwl_orig)*corr_spec; err_orig

   endif





   
   return, fiber
end