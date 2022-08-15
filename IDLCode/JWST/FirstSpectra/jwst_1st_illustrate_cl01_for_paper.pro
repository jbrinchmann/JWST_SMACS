pro jwst_1st_illustrate_cl01_fits_for_paper, ps=ps, grid_g16=grid_g16
   ;;
   ;; Create an illustration of the results of CL01 fitting for
   ;; the five high-z galaxies. I divide this into two: one for O/H
   ;; and another for other quantities. 
   ;;

   todo = ['04590', '05144', '06355', '08140', '10612']

   ;; Oxygen abundance from the direct method
   ab_oh = [[1.4560E-05,   2.0400E-05,   3.7370E-05], $
            [5.7469E-05,   8.7408E-05,   1.5198E-04], $
            [7.2083E-05,   9.7838E-05,   1.4168E-04], $
            [-9999, -9999, -9999.], $
            [4.9722E-05,   1.0905E-04 ,  2.8219E-04]]

   te_direct = [25209.3, 14455, 13922.5, -9999, 15370.4]/1e4
   

   trussler_sfr = [2.18, 0.32, 4.94, -9999., 0.94]

   ;; I correct these below
   tacchella_sfr = [4, -1, 38, -1, 7.0]

   
   ;; From Carnall et al - SHOULD UPDATE!
   lens_correction = [10.01, 2.89, 2.69, 1.67, 1.58]

   ;; From Trussler et al, except for 8140 where I use Carnall etal 
   lens_correction = [5.8, 3.2, 1.6, 1.67, 1.7]

   magn_tacchella = [3.74, -1, 1.23, -1, 1.34]
   for i=0L, 5-1 do $
    if (magn_tacchella[i] gt 0) then $
     tacchella_sfr[i] = tacchella_sfr[i]*magn_tacchella[i]/lens_correction[i]

;   stop
   
   keys = 'id'+todo
   n_obj = n_elements(todo)
   
   jwst_1st_dirs, dd
   outdir = dd.root+'CL01Fit/'
   restore, outdir+'combined-results.sav' ;; Gets d & s

   dim = hdf5_to_struct(outdir+'dim.h5')

   print, 'Depletion times: '
   for i=0L, n_elements(keys)-1 do begin
      ;; Convert to per kpc and assume emission within 1x1 kpc^2
      print, todo[i], 'Log years=', s[keys[i]].mugas.median+6-s[keys[i]].logsfr.median
   endfor
;   stop
   
   
   ;;
   ;; Show O/H, log U log SFR, \mu_gas, 
   ;; 

   to_show = ['OH', 'U', 'LOGSFR', 'MUGAS']
   labels = ['12 + Log O/H', 'Log U', 'Log SFR [M_{sun}/yr]', $
             'log \mu_{gas} [M_{sun}/pc^2]', '\tau_V']
   
   n_var = n_elements(to_show)

   an_col = '#2aba45'
   
   if keyword_set(ps) then begin
      ps_on, dd.figdir+'cl01-OH-comparison.ps', aspect=0.25, xlen=12, /times,/color
      !p.font = 0
      !x.thick = 2
      !y.thick = 2
      !p.thick = 4
      
   endif
   
   pos = get_position_arr(0, nx=n_obj, ny=1, xgap=0, ygap=0.07, $
                         ymax=0.9, ymin=0.15, xmin=0.05)
   count = 0L
   for iobj=0L, n_obj-1 do begin
      title = todo[iobj]

      tmp = d[keys[iobj]]
      ss = s[keys[iobj]]
      
      ;; Also load the corresponding G16 likelihood if relevant.
      p_G16 = jwst_1st_load_pifit_pdf(todo[iobj], 'OH', gridded=1)

      
      ;; Load the analogue file
      analogue_file = dd.root+'CL01Fit/cl01-analogues-'+string(format='(I5.5)', todo[iobj])+'.h5'
      if file_test(analogue_file) then begin
         h_an = hdf5_to_struct(analogue_file)
         have_an = 1
      endif else have_an = 0
;      stop

      pdf = struct_var(tmp, 'BIN_OH')
      pdf = pdf/max(pdf)
      x = struct_var(tmp, 'XOH')
      
      pos = get_position_arr(iobj, tickf=tf)
      plot, x, pdf, psym=10, position=pos, $
            ytickformat=tf[1], xtitle=TeXtoIDL('12+Log O/H'), $
            noerase=(iobj ne 0), $
            title=title, yrange=[0,1.15], $
            xrange=[7, 9.49], /xs

      oplot, p_G16.x, p_G16.pdf/max(p_G16.pdf), color=jb_colour('orange'), psym=10

      ;; Finally show the direct oxygen estimate.
      if (ab_oh[0, iobj] gt 0) then begin
         vv = 12+alog10(ab_oh[*, iobj])
         v_low = vv[0]
         v = vv[1]
         v_high = vv[1]
         ypos = 1.075
         symbols, 2 ,2
         myoploterr2, v, ypos, v_low, v_high, ypos, ypos, psym=8, $
                      color=jb_colour('CornFlowerBlue')
      endif

      if have_an then begin
         oplot, h_an.x_oh, h_an.pdf_oh/max(h_an.pdf_oh), color=jb_colour(an_col)
      endif
   endfor
   
   if keyword_set(ps) then begin
      ps_off
      my_cleanplot, /silent
   endif else stop

   ;;;
   ;; Next, we want to show several additional parameters
   ;;;

   grid_g16 = 1
   ;;
   ;; Show O/H, log U log SFR, \mu_gas, 
   ;; 

   to_show = ['U', 'LOGSFR', 'MUGAS', 'TE']
   labels = ['Log U', 'Log SFR [M_{sun}/yr]', $
             'log \mu_{gas} [M_{sun}/kpc^2]', 'T_e [10^4 K]']
   xranges = [[-3.9, -1], [-1, 1.85], [-1.3, 3.5], [8e3/1e4, 2.6e4/1e4]]
   
   n_var = n_elements(to_show)

   
   if keyword_set(ps) then begin
      ps_on, dd.figdir+'cl01-fit-results-multiple.ps', aspect=0.5, xlen=12, /times,/color
      !p.font = 0
      !x.thick = 2
      !y.thick = 2
      !p.thick = 4
      
   endif
   
   pos = get_position_arr(0, nx=n_obj, ny=n_var, xgap=0, ygap=0.07)
   count = 0L
   for iv=0L, n_var-1 do begin

      
      for iobj=0L, n_obj-1 do begin
         if iv eq 0 then $
          title = todo[iobj] $
         else title = ''
         
         tmp = d[keys[iobj]]
         ss = s[keys[iobj]]

         ;; Also load the corresponding G16 likelihood if relevant.
         show_pifit = 0
         if to_show[iv] eq 'OH' then begin
            p_G16 = jwst_1st_load_pifit_pdf(todo[iobj], 'OH', gridded=grid_g16)
            show_pifit = 1
         endif
         if to_show[iv] eq 'U' then begin
            p_G16 = jwst_1st_load_pifit_pdf(todo[iobj], 'LOGU', gridded=grid_g16)
            show_pifit = 1
         endif
         if to_show[iv] eq 'LOGSFR' then begin
            p_G16 = jwst_1st_load_pifit_pdf(todo[iobj], 'Log_SFR', gridded=grid_g16)
            show_pifit = 0      ; Not tested this enough to be sure
         endif
         if to_show[iv] eq 'MUGAS' and grid_g16 then begin
            p_G16 = jwst_1st_load_pifit_pdf(todo[iobj], 'Log_Sigma_gas', gridded=1)
            show_pifit = 1

         endif

         ;; Load the analogue file
         analogue_file = dd.root+'CL01Fit/cl01-analogues-'+string(format='(I5.5)', todo[iobj])+'.h5'
         if file_test(analogue_file) then begin
            h_an = hdf5_to_struct(analogue_file)
            have_an = 1
         endif else have_an = 0

         
         case to_show[iv] of
            'U': x = dim.logu
            'TAUV': x = dim.tauv
            else: begin
               x = struct_var(tmp, 'X'+to_show[iv])
            end
         endcase

         ;; Correct for lensing.
         if to_show[iv] eq 'LOGSFR' then $
          x = x-alog10(lens_correction[iobj])
         
         pdf = struct_var(tmp, 'BIN_'+to_show[iv])
         pdf = pdf/max(pdf)

         
         if to_show[iv] eq 'TE' then begin
            x = x/1e4
            ;; It was also calculated on too fine grid so we can add
            ;; some bins.
            x2 = []
            y2 = []
            for ii=0L, n_elements(x)-2, 2 do begin
               x2 = [x2, 0.5*(x[ii]+x[ii+1])]
               y2 = [y2, pdf[ii]+pdf[ii+1]]
            endfor
            x = x2
            pdf = y2/max(y2)
         endif



         
         pos = get_position_arr(count, tickf=tf)
         plot, x, pdf, psym=10, position=pos, $
               ytickformat=tf[1], xtitle=TeXtoIDL(labels[iv]), $
               noerase=(count ne 0), $
               title=title, yrange=[0,1.15], $
               xrange=xranges[*, iv], /xs

         if show_pifit then $
          oplot, p_G16.x, p_G16.pdf/max(p_G16.pdf), color=jb_colour('orange')

          if to_show[iv] eq 'LOGSFR' then begin
             if trussler_sfr[iobj] gt 0 then $
              oplot, alog10(trussler_sfr[iobj])+[0, 0], !y.crange, color=jb_colour('CornFlowerBlue'), $
                     linestyle=2

             if tacchella_sfr[iobj] gt 0 then $
              oplot, alog10(tacchella_sfr[iobj])+[0, 0], !y.crange, color=jb_colour('orange'), $
                     linestyle=2

             if have_an then $
              oplot, h_an.x_logsfr, h_an.pdf_logsfr/max(h_an.pdf_logsfr), color=jb_colour(an_col)
          endif
         ;; Finally show the direct oxygen estimate.
         if (to_show[iv] eq 'OH') then begin
            if (ab_oh[0, iobj] gt 0) then begin
               vv = 12+alog10(ab_oh[*, iobj])
               v_low = vv[0]
               v = vv[1]
               v_high = vv[1]
               ypos = 1.075
               symbols, 2 ,2
               myoploterr2, v, ypos, v_low, v_high, ypos, ypos, psym=8, $
                            color=jb_colour('CornFlowerBlue')
            endif

            if have_an then begin
               oplot, h_an.x_oh, h_an.pdf_oh/max(h_an.pdf_oh), color=jb_colour(an_col)
            endif
         endif

         if to_show[iv] eq 'TE' then begin
            if te_direct[iobj] gt 0 then begin
               oplot, te_direct[iobj]+[0, 0], !y.crange, color=jb_colour('purple'), thick=3
            endif

            ss = summarise_distribution(x, pdf)
            print, todo[iobj], te_direct[iobj], ss.median 

            
         endif
         
         if to_show[iv] eq 'MUGAS' then begin
            ;;
            ;; Also show the analogues
            ;;
            if have_an then begin
               oplot, h_an.x_mugas, h_an.pdf_mugas/max(h_an.pdf_mugas), color=jb_colour(an_col)
            endif

         endif
         if to_show[iv] eq 'TE' then begin
            ;;
            ;; Also show the analogues
            ;;
            if have_an then begin
               oplot, h_an.x_te/1e4, h_an.pdf_te/max(h_an.pdf_te), color=jb_colour(an_col)
            endif
;            stop
         endif
         if to_show[iv] eq 'U' then begin
            ;;
            ;; Also show the analogues
            ;;
            if have_an then begin
;               oplot, dim.logu, h_an.pdf_logu/max(h_an.pdf_logu), color=jb_colour(an_col)
            endif
;            stop
         endif

         count = count+1
         
      endfor
      

   endfor
   
   if keyword_set(ps) then begin
      ps_off
      my_cleanplot, /silent
   endif



   
end
