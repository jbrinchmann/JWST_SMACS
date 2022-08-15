pro jwst_1st_h2_lines_vs_it2016
   ;;
   ;; Show where the H2 lines fall wrt. Izotov & Thuan
   ;;

   jwst_1st_dirs, dd

   g = jwst_1st_load_one_gaussfit('02736_09483')

   i1 = where(strtrim(g.line) eq 'H2_10S1')

   x_ref = g.flux[i1[0]]
   dx_ref = g.dflux[i1[0]]
   
   get = ['H2_10S3', 'H2_10S2', 'H2_10S1']
   for i=0L, n_elements(get)-1 do begin
      ii = where(strtrim(g.line) eq get[i])
      x = g.flux[ii[0]]
      dx = g.dflux[ii[0]]

      r = x/x_ref
      dr = r*error_on_fraction(x, dx, x_ref, dx_ref)

      print, 'The ratio of '+get[i]+' to H2_10S1 is ', r, dr


   endfor

   
   

end
