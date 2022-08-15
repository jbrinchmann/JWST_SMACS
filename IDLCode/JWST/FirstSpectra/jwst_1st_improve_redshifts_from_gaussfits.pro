pro jwst_1st_improve_redshifts_from_gaussfits, ps=ps
   ;;
   ;; Use the individual Gaussian fits to improve the redshifts of the
   ;; galaxies.
   ;;

   n_mc = 5000
   
   jwst_1st_dirs, dd
   common linelist, lname, l_line, type
   if (n_elements(l_line) eq 0 or n_elements(type) eq 0) or keyword_set(reread) then begin
      jwst_initialize_labels
   endif

   outdir = dd.specdir+"Redshifts/"
   t_z = mrdfits(dd.root+'redshift-for-sources.fits', 1, /silent)


   print, 'This should not be run again without change as it compares to OLD redshifts'
   stop
   
   openw, lun, outdir+'refined-redshifts.txt', /get_lun
   for i=0L, n_elements(t_z)-1 do begin
      if (t_z[i].confidence lt 2) then continue
      
      t_g =  jwst_1st_load_one_gaussfit(t_z[i].object)

      ;;
      ;; Calculate the redshifts for each line.
      ;;
      ;; HOWEVER we do not use doublets as they have tied wavelengths 
      ;;
      n_lines = n_elements(t_g.line)

      z_lines = []

      skip = ['OII_3726', 'OII_3729', 'SII_6717', 'SII_6731', $
              'CIII_1907', 'CIII_1909', 'NEIV_2422', 'NEIV_2424', $
              'NII_6548', 'H_ALPHA']

      z_med = []
      z_p16 = []
      z_p84 = []
      i_use = []
      for i_line=0L, n_lines-1 do begin
         this_line = t_g.line[i_line]

         if any(skip eq this_line) then $
          continue
         
         i_match = where(strupcase(lname) eq t_g.line[i_line], n_match)
         if n_match eq 0 then stop
         i_match = i_match[0]


         ;; Monte Carlo this one.
         x_c = t_g.lcenter[i_line]+randomn(sss, n_mc)*t_g.dlcenter[i_line]
         x_z = x_c/l_line[i_match]-1.0
         vals = quantile(x_z, [0.16, 0.5, 0.84])
         z_lines = [z_lines, x_z]

         z_med = [z_med, vals[1]]
         z_p16 = [z_p16, vals[0]]
         z_p84 = [z_p84, vals[2]]
         i_use = [i_use, i_line]

      endfor

      ;; Now calculate statistics of these.
      vals = quantile(z_lines, [0.16, 0.5, 0.84])

      print, format='(A12," Z(med) = ",F0.7," [",F0.7,", ",F0.7,"]")', $
             t_z[i].object, vals[1], vals[0], vals[2]

      z_best = vals[1]
      z_low = vals[0]
      z_high = vals[2]


      ;; I'll also print out a symmetrised error.
      printf, lun, format='(A12,2X,5(F0.7,2X))', $
              t_z[i].object, vals[1], 0.5*(vals[2]-vals[0]), $
              vals[0]-vals[1], vals[2]-vals[1]
      

      old = 0

      if (old) then begin
         ;;
         ;; Create a simple visualisation of this.
         ;;
         ;; I will show the redshift as a function of the S/N of the line.
         ;;
         yrange = minmax([z_lines, t_z[i].redshift])
         sn = t_g.flux[i_use]/t_g.dflux[i_use]
         l_used = t_g.lcenter[i_use]/1e4 ; micron
         
         if keyword_set(ps) then begin
            ps_on, outdir+t_z[i].object+'_z_ill.ps', /landscape, /color, $
                   /times
            !p.font = 0
            !x.thick = 2
            !y.thick = 2
         endif
         
         pos = get_position_arr(0, nx=2, ny=1, xmin=xmin, xmax=xmax, $
                                ymin=ymin, ymax=ymax, tickformat=tf, xgap=0.01)
         plot, sn, z_med, /nodata, xtitle='S/N', ytitle='Redshift', $
               charsize=1.6, yrange=yrange, xtickformat=tf[0], $
               ytickformat=tf[1], position=pos
         symbols, 2, 1.2
         myoploterr2, sn, z_med, sn, sn, z_p16, z_p84, psym=8
         
         ;; And the best fit.
         oplot, !x.crange, z_best+[0, 0], color=jb_colour('red'), $
                thick=3
         oplot, !x.crange, z_low+[0, 0], color=jb_colour('red'), $
                thick=3, linestyle=2
         oplot, !x.crange, z_high+[0, 0], color=jb_colour('red'), $
                thick=3, linestyle=2
         
         
         ;; And my original redshift overplotted
         oplot, !x.crange, t_z[i].redshift+[0, 0], linestyle=2, $
                color=jb_colour('CornFlowerBlue'), thick=3
         
         pos = get_position_arr(1, tickf=tf)
         plot, l_used, z_med, /nodata, xtitle='Wavelength [micron]', $
               charsize=1.6, yrange=yrange, xtickformat=tf[0], $
               ytickformat=tf[1], position=pos, /noerase, xrange=[1.7, 5.2]
         symbols, 2, 1.2
         myoploterr2, l_used, z_med, l_used, l_used, z_p16, z_p84, psym=8
         
         ;; And the best fit.
         oplot, !x.crange, z_best+[0, 0], color=jb_colour('red'), $
                thick=3
         oplot, !x.crange, z_low+[0, 0], color=jb_colour('red'), $
                thick=3, linestyle=2
         oplot, !x.crange, z_high+[0, 0], color=jb_colour('red'), $
                thick=3, linestyle=2
         
      
         ;; And my original redshift overplotted
         oplot, !x.crange, t_z[i].redshift+[0, 0], linestyle=2, $
                color=jb_colour('CornFlowerBlue'), thick=3
         
         xyouts, 0.5*(xmax+xmin), ymax+0.02, t_z[i].object, charsize=2, /normal, $
                 align=0.5
         
         
         if keyword_set(ps) then begin
            ps_off
            my_cleanplot, /silent
         endif else stop
      endif

   endfor
   
   free_lun, lun
   

end
