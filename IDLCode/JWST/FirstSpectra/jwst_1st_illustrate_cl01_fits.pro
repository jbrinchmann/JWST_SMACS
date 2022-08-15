pro jwst_1st_illustrate_cl01_fits , ps=ps, grid_g16=grid_g16
   ;;
   ;; Create a simple illustration of the results of CL01 fitting for
   ;; the five high-z galaxies.
   ;;

   todo = ['04590', '05144', '06355', '08140', '10612']

   ;; Oxygen abundance from the direct method


   ab_oh = [[1.4560E-05,   2.0400E-05,   3.7370E-05], $
            [5.7469E-05,   8.7408E-05,   1.5198E-04], $
            [7.2083E-05,   9.7838E-05,   1.4168E-04], $
            [-9999, -9999, -9999.], $
            [4.9722E-05,   1.0905E-04 ,  2.8219E-04]]
   

   trussler_sfr = [2.18, 0.32, 4.94, -9999., 0.94]
   
   ;; From Carnall et al - SHOULD UPDATE!
   lens_correction = [10.01, 2.89, 2.69, 1.67, 1.58]

   ;; From Trussler et al, except for 8140 where I use Carnall etal 
   lens_correction = [5.8, 3.2, 1.6, 1.67, 1.7]

   keys = 'id'+todo
   n_obj = n_elements(todo)
   
   jwst_1st_dirs, dd
   outdir = dd.root+'CL01Fit/'
   restore, outdir+'combined-results.sav' ;; Gets d & s

   dim = hdf5_to_struct(outdir+'dim.h5')
;   stop
   
   ;;
   ;; Show O/H, log U log SFR, \mu_gas, 
   ;; 

   to_show = ['OH', 'U', 'LOGSFR', 'MUGAS']
   labels = ['12 + Log O/H', 'Log U', 'Log SFR [M_{sun}/yr]', $
             'log \mu_{gas} [M_{sun}/kpc^2]', '\tau_V']
   xranges = [[7.0, 9.5], [-3.9, -1], [-1, 1.49], [-1.3, 3.5], [0, 0.99]]
   
   n_var = n_elements(to_show)

   an_col = 'HotPink'
   
   if keyword_set(ps) then begin
      ps_on, dd.figdir+'cl01-fit-results.ps', aspect=0.5, xlen=12, /times,/color
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

         if to_show[iv] eq 'MUGAS' then begin
            ;;
            ;; Also show the analogues
            ;;
            if have_an then begin
               oplot, h_an.x_mugas, h_an.pdf_mugas/max(h_an.pdf_mugas), color=jb_colour(an_col)
            endif

         endif

         count = count+1
         
      endfor
      

   endfor
   
   if keyword_set(ps) then begin
      ps_off
      my_cleanplot, /silent
   endif

end
