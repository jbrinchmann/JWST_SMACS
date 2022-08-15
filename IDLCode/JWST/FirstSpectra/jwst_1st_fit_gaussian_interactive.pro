pro jwst_1st_fit_gaussian_interactive, wave, flux, dflux, z, result=result, $
                                       LSF_R=LSF_R

   ;;
   ;; Interactive fitting. Let the user step through as they see fit.
   ;; This will work better if the spectrum is continuum subtracted. 
   ;;

   common linelist, lname, l_line, type

   if (n_elements(l_line) eq 0 or n_elements(type) eq 0) or keyword_set(reread) then begin
      jwst_initialize_labels
   endif

   
   l_oiia = 3727.09  ; 3726.032 
   l_oiib = 3729.88  ; 3728.815 
   l_ciiia = 1906.68
   l_ciiib = 1908.73
   l_siia = 6718.29
   l_siib = 6732.67
   l_neiva = 2421.83
   l_neivb = 2424.42
   
   l_ha = 6564.61
   l_n2a = 6549.85
   l_n2b = 6585.28
   
   
   sep_oii = l_oiib-l_oiia
   sep_sii = l_siib-l_siia
   sep_ciii = l_ciiib-l_ciiia
   sep_neiv = l_neivb-l_neiva

   ;; Wavelengths to add to get to the other line
   sep_ha_a = l_n2a-l_ha
   sep_ha_b = l_n2b-l_ha
   
   
   ans = ''

   specfit = wave*0.0           ; We add the gaussians to this.

   ;;
   ;; Create an output structure
   ;;
   result = {lines: lname}


   xrange = minmax(wave/1e4)
   xrange_orig = xrange
   
   my_lineid_plot, wave/1e4, flux, l_line*(1+z)/1e4, lname, /extend, $
                   elinestyle=1, lcolor=jb_colour('yellow'), $
                   xrange=xrange, title=title, xtitle='Wavelength [micron]'

;   stop
   
   while (ans ne 'q') do begin
      print, 'f: fit a line, q: quit, h: halt, x: xrange, xr: reset x-range'
      read, ans

      case ans of
         'f': begin
            print, 'Click on the line!'
            cursor, xx_obs, yy

            xx_obs = xx_obs*1e4 ; Convert to Angstrom
            xx_rest = xx_obs/(1+z)
            
            ;; Find the closest line to this point. 
            dum = min(abs(l_line-xx_rest), i_close)
            current_line = strupcase(lname[i_close])
            print, "Fitting "+current_line

            ;;;
            ;; Doublets need some special care. The user will click
            ;; closest to one of them and we need to handle this
            ;; gracefully
            ;;;
            is_ha = 0
            is_o2 = 0
            is_c3 = 0
            is_s2 = 0
            is_ne4 = 0
            doublet_lines = []
            if (current_line eq 'OII_3726' or $
                current_line eq 'OII_3729') then begin
               is_o2 = 1
               sep = sep_oii
               doublet_lines = ['OII_3726', 'OII_3729']
            endif else if (current_line eq 'SII_6717' or $
                           current_line eq 'SII_6731') then begin
               is_s2 = 1
               sep = sep_sii
               doublet_lines = ['SII_6717', 'SII_6731']
            endif else if (current_line eq 'CIII_1907' or $
                           current_line eq 'CIII_1909') then begin
               is_c3 = 1
               sep = sep_ciii
               doublet_lines = ['CIII_1907', 'CIII_1909']
            endif else if (current_line eq 'NEIV_2422' or $
                           current_line eq 'NEIV_2424') then begin
               is_ne4 = 1
               sep = sep_neiv
               doublet_lines = ['NEIV_2422', 'NEIV_2424']
            endif else if (current_line eq 'H_ALPHA') then begin
               is_ha = 1
               sep_a = sep_ha_a
               sep_b = sep_ha_b
               triplet_lines = ['NII_6548', 'H_ALPHA', 'NII_6584']
            endif
            is_doublet = (is_o2 or is_s2 or is_c3 or is_ne4)
            is_triplet = (is_ha)
            if is_doublet then begin
               is_multiplet = 1
               multiplet_lines = doublet_lines
               width = sep*(1+z)*1.5
            endif else if is_triplet then begin
               is_multiplet = 1
               multiplet_lines = triplet_lines
               width = (sep_b-sep_a)*(1+z)*1.5
            endif else begin
               is_multiplet = 0
               width = 100
            endelse
            ;; We first fit in observed frame and then refine

            ;; Refine the position by fitting a Gaussian.
            use = where(abs(wave - xx_obs) lt width and finite(flux) eq 1)
            xx = wave[use]
            yy = flux[use]
            dyy = dflux[use]

            if n_elements(LSF_R) gt 0 then begin
               fwhm = mean(xx/LSF_R[use])
               intrinsic_sigma = fwhm/2.355 ; 2 sqrt(2 ln 2)
            endif else intrinsic_sigma = 0.0

            
            if is_doublet then begin
               g = fit_two_gaussians_mp_alt(xx, yy, a, $
                                            separation=sep*(1+z), $
                                            error=dyy, perr=da, $
                                            intrinsic_sigma=intrinsic_sigma, $
                                            /return_flux)
               
            endif else if is_triplet then begin
               g = fit_ha_complex_mp(xx, yy, dyy, a, $
                                     sep_a=sep_a*(1+z), $
                                     sep_b=sep_b*(1+z), $
                                      perr=da, $
                                     intrinsic_sigma=intrinsic_sigma, $
                                     /return_flux)

;               stop
            endif else begin
               g = fit_one_gaussian_mp(xx, yy, a, /return_flux, $
                                       err=dyy, perror=da, $
                                       intrinsic_sigma=intrinsic_sigma)
            endelse


            
            ;;
            ;; Store results.
            ;;
            ;; We store flux, central wavelength, width, and
            ;; background and uncertainties on these
            ;;
            sig = sqrt(intrinsic_sigma^2+a[2]^2)
            tmp = {flux: a[0], dflux:  da[0], lc: a[1], dlc: da[1], $
                   sigma: sig, dsigma: da[2], level: a[3], dlevel: da[3]}

            
            if is_doublet then begin
               ;; Doublets must be treated differently because we must
               ;; initialise both lines.
               
               for i_d=0L, n_elements(doublet_lines)-1 do begin
                  if not tag_exist(result, doublet_lines[i_d]) then $
                   result = create_struct(result, doublet_lines[i_d], tmp)

                  ;; The tmp struct must be recreated for doublets
                  ;;
                  tmp.flux = a[0+i_d*3]
                  tmp.dflux = da[0+i_d*3]
                  tmp.lc = a[1]+i_d*sep
                  tmp.level = a[n_elements(a)-1]
                  tmp.dlevel = da[n_elements(a)-1]
                  struct_set_var, result, doublet_lines[i_d], tmp
                  
               endfor
            
            endif else if is_ha then begin
               ;; Special handling.
               tmp.level = a[n_elements(a)-1]
               tmp.dlevel = da[n_elements(a)-1]
               for i_d=0L, n_elements(triplet_lines)-1 do begin
                  
                  if not tag_exist(result, triplet_lines[i_d]) then $
                   result = create_struct(result, triplet_lines[i_d], tmp)
               endfor
               result.h_alpha.flux = a[0]
               result.h_alpha.dflux = da[0]
               result.h_alpha.lc = a[1]

               result.nii_6548.flux = a[3]/2.95732
               result.nii_6548.dflux = da[3]
               result.nii_6548.lc = a[1]+sep_a

               result.nii_6584.flux = a[3]
               result.nii_6584.dflux = da[3]
               result.nii_6584.lc = a[1]+sep_b
;               stop
               
            endif else if is_triplet then stop $
            else begin            
               if not tag_exist(result, current_line) then $
                result = create_struct(result, current_line, tmp)
            
               struct_set_var, result, current_line, tmp
            endelse


            ;; Update the spectrum fit - note that this does not
            ;; handle overlapping lines!
            specfit[use] = g
            oplot, wave/1e4, specfit, color=jb_colour('orange')
         end
         'h': stop
         'xr': begin
            xrange = xrange_orig
         end
         'x': begin
            print, 'Xrange to show: '
            read, x1, x2
            xrange = [x1, x2]
            my_lineid_plot, wave/1e4, flux, l_line*(1+z)/1e4, lname, /extend, $
                            elinestyle=1,lcolor=jb_colour('yellow'), $
                            xrange=xrange, title=title, xtitle='Wavelength [micron]'
            oplot, wave/1e4, specfit, color=jb_colour('orange')
         end
         'r': begin
            my_lineid_plot, wave/1e4, flux, l_line*(1+z)/1e4, lname, /extend, $
                            elinestyle=1,lcolor=jb_colour('yellow'), $
                            xrange=xrange, title=title, xtitle='Wavelength [micron]'
            oplot, wave/1e4, specfit, color=jb_colour('orange')
         end
         'o': begin
            ;; Show overview of what lines have been fit.
            tn = tag_names(result)
            if n_elements(tn) eq 1 then begin
               print, 'No lines yet fit'
            endif else begin
               for i=0L, n_elements(tn)-1 do begin
                  key = strupcase(tn[i])
                  if key eq 'LINES' then continue ; skip this
                  
                  x = struct_var(result, key)
                  print, format='(A10," S/N=",F0.2)', key, x.flux/x.dflux
               endfor
            endelse
         end
         'q': begin
            if (n_elements(bad_x) gt 0) then begin
               si = sort(bad_x)
               ui = uniq(bad_x[si])
               bad_x = bad_x[si[ui]]
               print, 'You marked ', n_elements(bad_x), ' pixels as bad.'
            endif else bad_x = -1
         end
      endcase
   endwhile

end
