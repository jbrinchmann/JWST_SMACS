pro jwst_1st_compare_abundances, ps=ps
   ;;
   ;; Compare my O/H estimates with the literature. I keep this in a
   ;; Numbers document which I then export to CSV. 
   ;;
   
   jwst_1st_dirs, dd
   t = r_read_table(dd.root+'TeModeling/oxygen-abundances-literature.csv', sep=',', /header)

   xpos = [1, 2, 3, 4, 5]       ; Locations for each object.
   shift = [-0.2, -0.1, 0.0, 0.1, 0.2, 0.15]/2.

   ymax = 0.85
   pos = get_position_arr(0, nx=2, ny=1, xgap=0.06, ymax=0.85, $
                          xmin=xmin, xmax=xmax, ymin=ymin)


   if keyword_set(ps) then begin
      ps_on, dd.figdir+'abundance-comp-to-literature.ps', aspect=0.4, xlen=20, $
             /times, /color
      !p.font = 0
      !x.thick = 2
      !y.thick = 2
      !p.thick = 4
      cz = 1.6
      scl = 0.8
   endif else begin
      !p.background = 'ffffff'x
      !p.color = '000000'x
      cz = 1.3
      scl = 1.0
   endelse
   
   plot, [0, 1], [0, 1], xrange=[0, 6], xstyle=1+4, $
         yrange=[1, 4], ytitle=TeXtoIDL('T_e [10^4 K]'), $
         charsize=cz, /nodata, position=pos


   ;; open circle, box, diamond, triangle, filled circle, open box
   sym = [1, 30, 31, 32, 20, 2]
   symsize = [2, 0.6, 0.7, 0.8, 0.6, 2]*scl
   ;; Colour-blind-friendly
   colors = ['#F71467', '#1E88E5', '#FFC107', '#009E73', '#88CCEE', '#004D40']
   order = ['TE_TRUMP', 'TE_SCHAERER', 'TE_RHOADS', 'TE_CURTI', 'TE_AC', 'TE_MINE']
   
   for i=0L, n_elements(t.object)-1 do begin

      for j=0l, n_elements(order)-1 do begin
         symbols, sym[j], symsize[j]
         y = struct_var(t, order[j])
         y = y[i]/1e4
         if y gt 0 then $
          plots, xpos[i]+shift[j], y, psym=8, color=jb_colour(colors[j]), symsize=2

         print, 'Doing '+order[j]+' color='+colors[j]
      endfor
      if i lt 5 then $
       oplot, [xpos[i], xpos[i]]+0.5, [!y.crange[0]+0.1, !y.crange[1]], linestyle=2, color=jb_colour('gray')
      
      xyouts, xpos[i], 1., string(format='(I0)', t.object[i]), align=0.5, $
              charsize=cz
              
      
   endfor

   pos = get_position_arr(1)
   plot, [0, 1], [0, 1], xrange=[0, 6], xstyle=1+4, $
         yrange=[6.5, 9], ytitle='12 + Log O/H', charsize=cz, $
         /nodata, /noerase, position=pos


   ;; I'll handle mine a bit differently
   order = ['OH_TRUMP', 'OH_SCHAERER', 'OH_RHOADS', 'OH_CURTI', 'OH_AC']
   
   for i=0L, n_elements(t.object)-1 do begin

      for j=0l, n_elements(order)-1 do begin
         symbols, sym[j], symsize[j]
         y = struct_var(t, order[j])
         y = y[i]
         if y gt 0 then $
          plots, xpos[i]+shift[j], y, psym=8, color=jb_colour(colors[j]), symsize=2 $
         else if y gt -10 then begin
            ;; Upper limits
            plots, xpos[i]+shift[j], abs(y), psym=8, color=jb_colour(colors[j]), symsize=2
            symbols, 6, 1
            plots, xpos[i]+shift[j], abs(y), psym=8, color=jb_colour(colors[j])
         endif
      endfor
      
      ;; Deal with mine a bit differently because I want to illustrate
      ;; the effect of the T([O II]) calibration adopted.
      symbols, sym[5], symsize[5]
      yvals = [t.oh_mine_o2Izotov[i], t.oh_mine[i], t.oh_mine_Pilyugin[i]]
      oplot, xpos[i]+shift[j]+[0, 0, 0], yvals
      plots, xpos[i]+shift[j], yvals[0], psym=8, color=jb_colour(colors[j])
      plots, xpos[i]+shift[j], yvals[1], psym=8, color=jb_colour(colors[j]), symsize=2
      symbols, 1, symsize[4]
      plots, xpos[i]+shift[j], yvals[2], psym=8, color=jb_colour(colors[j])
      
      
      
      if i lt 5 then $
       oplot, [xpos[i], xpos[i]]+0.5, [!y.crange[0]+0.1, !y.crange[1]], linestyle=2, color=jb_colour('gray')
      
      xyouts, xpos[i], 6.5, string(format='(I0)', t.object[i]), align=0.5, $
              charsize=cz
              
      
   endfor

   ;; Insert a legend in a plot on the top.
;;   plot, [0, 1], [0, 1], /nodata, xstyle=4, ystyle=4, xrange=[0, 5]
   

   

   
   if keyword_set(ps) then ps_off
   my_cleanplot, /silent
   
end
   
