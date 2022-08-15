;; Useful commands:
;;
;; my_lineid_plot, l/(1+z_guess), spec, l_line, lname, /extend,
;;elinestyle=1, lcolor='00ffff'x, dflux=dspec, xrange=[4500, 5100]
;;
;; muse_show_smooth_nb, l_view, root=root, name=name, data={cube: cs}, /sn, /remove, smooth=1.5

pro jwst_view_spec_findz, d, result=result, z=z, $
                          wavelet=wavelet, smooth=smooth, $
                          info=info, reread=reread

   ;;
   ;; A program to quickly inspect JWST NIRSpec spectra and get their
   ;; redshifts. 
   ;;

   
   common linelist, lname, l_line, type

   if (n_elements(l_line2) eq 0 or n_elements(type) eq 0) or keyword_set(reread) then begin
      readcol, getenv('HOME')+'/IDL/platefit/etc/linelist-jwst-full.txt', lname1, $
               l_line1, llimit, ulimit, type1, format='(A,F,F,F,A)', $
               /silent
      
      readcol, getenv('HOME')+'/IDL/platefit/etc/uv-linelist-strong.txt', lname2, $
               l_line2, lower, upper, type2, ref, format='(A,F,F,F,A,I)', $

               /silent
      type2[*] = 'is'           ; NOT CORRECT!!!
      type = [type1, type2]
      l_line = [l_line1, l_line2]
      lname = [lname1, lname2]
   endif


   sp = d.oned
   twod = d.twod

   dims = size(twod.flux, /dimen)
   nx = dims[0]
   ny = dims[1]
   n_lambda = nx
   
   spec = sp.flux
   dspec = sp.flux_error
   l = sp.wavelength*1e4

   ;; JWST spectra are always in vacuum
   vacuum = 1

   ;; Not sure I ever need the spectrum in air, but here it is at least.
   lair = l
   vactoair, lair
   
   if (keyword_set(wavelet)) then begin
      w = jb_atrous(spec)
      spec = spec-total(w.w[0:wavelet-1, *], 1) ; Remove high frequency noise
;      stop
   endif else if (n_elements(smooth) gt 0) then begin
      spec = wr_smooth(l, spec, smooth=smooth)
      dspec = wr_smooth(l, dspec, smooth=smooth)
   endif

   spec_orig = spec

   ;;
   ;; Redshift?
   ;;
   if (n_elements(z) gt 0) then $
    l = l/(1+z)
   

   ;;
   ;; Loop over plots - this is the part I really should have in
   ;; cgzplot, but no time right now. 
   ;;

   tlb = Widget_Base(Title='My Program')
   jb_cgzplot, l, spec, xsize=1300, ysize=700, parent=tlb, dy=dspec
   Widget_Control, tlb, /Realize
   xmanager, 'myprogram', tlb

   ans = 'y'
   z_guess = -9999.
   plot, l, spec, xrange=!x.crange
   oplot, l, dspec, color=jb_colour('orange')
   oplot, l, spec
   print, 'Lines: '
   print, '1: Ly-a, 2: C III] 3: [O II], 4: Hb, 5: [O III]5007, 6: Ha - otherwise wavelength'

   ;; 28/5/2016
   ;; Changed line wavelengths to be in vacuum.

   lambda_line = [1215.67, 1909, 3727, 4862.68, 5008.24 , 6564.61]

   l_oiia = 3727.09  ; 3726.032 
   l_oiib = 3729.88  ; 3728.815 
   l_ciiia = 1906.68
   l_ciiib = 1908.73
   sep_oii = l_oiib-l_oiia
   sep_ciii = l_ciiib-l_ciiia

   pos1 = [0.1, 0.12, 0.75, 0.95]
   pos2 = [0.755, 0.22, 0.95, 0.85]
   first = 1
   
   while (ans ne 'q') do begin
      
      print, 'i: Identify line, s: smooth, o: output, q: Quit, h: halt, f: check fit, p: plot all, a: run autoz'
      read, ans

      ansfull = strlowcase(ans)
      ans = strlowcase(strmid(ans, 0, 1))

      case ans of
         'a': begin
            
            ;; Run auto-z on the spectrum.

            if (first) then begin
               set_up_log_lambda, gap, log_lambda_rebin, velocity_gap=velocity_gap, $
                                  /highz
;               tdata_rebin = rebin_template_data(log_lambda_rebin, file=getenv('HOME') + $
;                                                 '/IDL/AutoZ/data/filtered-templates.fits')
               templates_to_use = [2+indgen(13), 16+indgen(7), 28, 29];, 33]
               tdata_rebin = rebin_template_data(log_lambda_rebin, templates_to_use, $
                                                 file=getenv('HOME') + $
                                                 '/IDL/AutoZ/data/templates-incl-tc-2016-03-22.fits')
               first = 0
            endif


            a_z = autoz_singlespec(spec, dspec, l, $
                                   num_peaks=num_peaks, $
                                   log_lambda_rebin=log_lambda_rebin, $
                                   gap=gap, tdata_rebin=tdata_rebin, $
                                    helio_velocity=0.0, $
                                   inv_correction=inv_correction, /high)
            autoz = convert_z_format(a_z)  
            print, format='("AUTOZ best redshift=",F0.5,"   p=",F0.4)', autoz.z, autoz.prob
            stop
            
            l_rest = l/(1+autoz.z)
            plot, l_rest, spec, /nodata,  $
                  xtitle='Wavelength', ytitle='Flux'
            for iline=0L, n_elements(lname)-1 do begin
               oplot, l_line[iline]+[0, 0], !y.crange, color=jb_colour('gray55')
            endfor
            oplot, l_rest, spec 
            
;            stop
         end
         'i': begin
            print, 'Click on the line!'
            cursor, xx, yy
            l_obs = xx

            print, format='($,"line: ")'
            read, line

            if (line lt 10) then $
             line_rest = lambda_line[line-1] $
            else $
             line_rest = line
            
            ;; We first fit in observed frame and then refine

            ;; First guess z - might be off a bit.
            z_guess = l_obs/line_rest-1.0

            ;; Refine the position by fitting a Gaussian.
            use = where(abs(l - l_obs) lt 15 and finite(spec) eq 1)
            xx = l[use]
            yy = spec[use]

            if (line_rest eq 3727) then begin
               ;; Special treatment for [O II]
               g = fit_two_gaussians_mp(xx, yy, a, $
                                        separation=sep_oii*(1+z_guess))
            endif else if (line_rest eq 1909) then begin
               g = fit_two_gaussians_mp(xx, yy, a, $
                                        separation=sep_ciii*(1+z_guess))
            endif else begin
               g = fit_one_gaussian_mp(xx, yy, a)
            endelse
            
            ;; Calculate redshift.
            if (line_rest eq 3727) then $
             z_guess = a[1]/l_oiib-1.0 $
            else if (line_rest eq 1909) then $
             z_guess = a[1]/l_ciiib-1.0 $
            else $
             z_guess = a[1]/line_rest

            l_fit = a[1]
;            stop
            ;; Redo fit.
            if (line_rest eq 3727) then begin
               ;; Special treatment for [O II]
               g = fit_two_gaussians_mp(xx, yy, a, $
                                        separation=sep_oii*(1+z_guess))
            endif else if (line_rest eq 1909) then begin
               g = fit_two_gaussians_mp(xx, yy, a, $
                                        separation=sep_ciii*(1+z_guess))
            endif else begin
               g = fit_one_gaussian_mp(xx, yy, a)
            endelse
            if (line_rest eq 3727) then $
             z_guess = a[1]/l_oiib-1.0 $
            else if (line_rest eq 1909) then $
             z_guess = a[1]/l_ciiib-1.0 $
            else $
             z_guess = a[1]/line_rest-1.0

            print, format='("  Z=",F8.5, "  for line at l=",F7.2," with wavelength ",F7.2," and EW=",F7.2,"A")', z_guess, $
                   a[1], line_rest, a[0]/a[3]

            l_rest = l/(1+z_guess)
            plot, l_rest, spec, /nodata, position=pos1, $
                  xtitle='Wavelength', ytitle='Flux'
            for iline=0L, n_elements(lname)-1 do begin
               oplot, l_line[iline]+[0, 0], !y.crange, color=jb_colour('gray55')
            endfor
            
            oplot, l_rest, spec

         end
         'f': begin
            if (n_elements(g) eq 0) then begin
               print, 'You need to identify a line first!'
            endif else begin
               use = where(abs(l - l_obs) lt 80)
               plot, l[use], spec[use]
               oplot, xx, g, color=jb_colour('red')
            endelse
         end
         'o': begin
            file = ''
            print, format='($,"Filename: ")'
            read, file
            
            openw, lun, file, /get_lun
            if (n_elements(l_rest) eq 0) then $
             lr = l $
            else lr = l_rest

            w = jb_atrous(spec)


            if (n_elements(z_guess) gt 0) then $
             printf, lun, format='("#  Z = ",F8.5)', z_guess
            printf, lun, '#'
            printf, lun, '#   (1): Observed wavelength'
            printf, lun, '#   (2): Rest-frame wavelength'
            printf, lun, '#   (3): Flux'
            printf, lun, '#   (4): Flux error'
            printf, lun, '#   (5): Wavelength reconstructed flux leaving out the highest freq.'
            printf, lun, '#   (6): Ditto, but leaving out the two highest freq.'
            printf, lun, '#'
            printf, lun, '#   (1)      (2)         (3)          (4)          (5)          (6)'
            
            f = '(2X,F7.2,2x,F7.2,2X,4(E11.4,2X))'
            for j=0L, n_elements(l)-1 do begin
               printf, lun, format=f, l[j], lr[j], spec[j], dspec[j], $
                       spec[j]-w.w[0, j], $
                       spec[j]-w.w[0, j]-w.w[1, j]
            endfor
            free_lun, lun
         end

         'h': stop
         's': begin
            print, format='($,"smooth: ")'
            read, smooth
            if (smooth le 0) then spec = spec_orig $
            else spec = wr_smooth(l, spec, smooth=smooth)
         end
         'x': begin
            ;; Show a JWST image - not yet implemented!
            stop
            if n_elements(info) gt 0 then begin

               if tag_exist(info, 'HST') then begin
                  if ansfull eq 'xv' then begin
                     pos = get_position_arr(0, nx=2, ny=1, xgap=0.1)
                     if (n_elements(l_fit) eq 0) then begin
                        print, format='($,"Wavelength to create NB image around: ")'
                        read, l_view
                     endif else l_view = l_fit
                     muse_show_smooth_nb, l_view, root=root, name=name, /remove,/sn, $
                                          data={cube: cs}, position=pos , dl=nb_width
                     overlay_cross, color='wheat'
                    noerase = 1
                     pos = get_position_arr(1)
                     tvim_true, asinh(info.hst_im), noerase=noerase, position=pos
                     overlay_cross, color='wheat'

                  endif else begin
                     noerase = 0
                     tvim_true, info.hst_im

                  endelse
                  
                  plots, info.x_hst, info.y_hst, psym=6
                  xyouts, info.x_hst, info.y_hst+1.5, string(format='(I0)', info.ids), $
                          align=0.5, color=jb_colour('ForestGreen')
                  
                  if ansfull eq 'xv' then $
                   !p.multi = 0
                  
               endif else print, 'No HST image found in INFO structure!'
            endif
         end
         'q': begin
            result = {z: z_guess}
            if (n_elements(a) gt 0) then $
             result = create_Struct(result, 'linefit', a)
         end
         else: begin
            print, 'Option not recognised!'
         end
      endcase

   endwhile
   

end
