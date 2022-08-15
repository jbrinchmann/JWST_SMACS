pro jwst_1st_show_line_ratio_variations, lname, values, dvalues, norm_index, $
 labels=labels, _extra=_extra


   ;; Show line ratio variations

   ;; This is the flux we want to normalise against. 
   f_norm = reform(values[norm_index, *])
   df_norm = reform(dvalues[norm_index, *])

   n_mc = 5000

   
   
   plot, [-1, n_elements(lname)], [-5, 5], xstyle=4, /nodata, $
         _extra=_extra

   
   dims = size(values, /dimen)
   n_pfs = dims[1]
   c = model_grid_colors(n_grid=n_pfs+2, table=33)
   for i=0L, n_elements(lname)-1 do begin
      l_r = reform(values[i, *]/f_norm)
      l_r = alog10(l_r)

      
      

      for j=0L, n_elements(l_r)-1 do begin
         if (values[i, j] lt -100) or (values[i, j] eq 0.0) then continue


         r_x_norm = f_norm[j]+randomn(sss, N_mc)*df_norm[j]
         

         
         r_y = values[i, j]+randomn(sss, n_mc)*dvalues[i, j]
         r_log_ratio = alog10(r_y/r_x_norm)
         ok = where(finite(r_log_ratio) eq 1, n_ok)

;         stop
         if (n_ok lt 100) then continue
         vals = quantile(r_log_ratio[ok], [0.16, 0.5, 0.84])
;         stop
      
         l_low = vals[0]
         l_med = vals[1]
         l_high = vals[2]

         xpos = i+randomu(sss, 1)*0.3
         
         plots, xpos, l_med, color=c.g[j+1], psym=8, symsize=1.5
;         stop
         myoploterr2, xpos, l_med, xpos, xpos, l_low, l_high, psym=8, $
                      errcolor=c.g[j+1], color=c.g[j+1]
      endfor

;      print, minmax(l_r)
      xyouts, i, !y.crange[0], lname[i], orient=90

      
   endfor

   
   if n_elements(labels) gt 0 then begin
      ;; Label what the colours mean
      ystart = 0.6
      dy = 0.2
      ddy = -0.05
      xstart = -4
      dx = 0.2
      for i=0L, n_elements(l_r)-1 do begin
         plots, xstart, ystart-dy*i, psym=8, color=c.g[i+1], symsize=2
         xyouts, xstart+dx, ystart-dy*i+ddy, labels[i], charsize=1.6
      endfor
   endif
 

end
