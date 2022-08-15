pro jwst_show_spectrum_with_filters, sp, filters=filters, $
                                     scale=scale, _extra=_extra
   ;;
   ;; OVerplot filters on top of a NIRSpec spectrum
   ;;
   
   if n_elements(scale) eq 0 then scale = 1

   if N_elements(filters) eq 0 then $
    filters = ['F090W', 'F150W', 'F200W', 'F277W', 'F356W', 'F444W']

   f = sp.flux
   f = f/quantile(f, 0.99)
   
   plot, sp.wavelength, f, _extra=_extra

   c = model_grid_colors(n_grid=n_elements(filters)+2, table='RdBu')
   for i=0L, n_elements(filters)-1 do begin
      f = jwst_get_filter_data(filters[i], 'NIRCAM')
      x = f.wavelength
      y = f.throughput*scale

      
      oplot, x, y, color=c.g[i], thick=4

      ;; Find a location to label
      high = where(y gt 0.7*max(y))
      i_left = min(high)
      i_right = max(high)
      i_top = ceil(0.5*(i_left+i_right))
      
      xyouts, x[i_top], y[i_top]+0.15, filters[i], $
              color=c.g[i], charsize=2
      
      
   endfor
   
   

end
