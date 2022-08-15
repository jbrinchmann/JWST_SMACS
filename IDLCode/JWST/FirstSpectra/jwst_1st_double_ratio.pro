function jwst_1st_double_ratio, t, line1, line2, line3, line4, $
                                n_mc=n_mc, err_scale=err_scale, $
                                norm_line=norm_line, min_sn=min_sn

   ;;
   ;; Calculate alog10((x1/x2)/(x3/x4)) using an MC routine 
   ;;
   if n_elements(err_scale) eq 0 then err_scale = 2.0
   if n_elements(n_mc) eq 0 then n_mc = 10001
   if n_elements(norm_line) eq 0 then norm_line = line1
   if (n_elements(min_sn) eq 0) then min_sn = 5.0

   lines = [line1, line2, line3, line4]
   
   i_norm = where(lines eq norm_line, n_norm)
   if n_norm eq 0 then stop
   
   fluxes = dblarr(n_elements(t), n_elements(lines))
   dfluxes = fluxes
   for i=0, n_elements(lines)-1 do begin
      if lines[i] eq 'OII3727' then begin
         fluxes[*, i] = struct_var(t, 'oii_3726_flux')+struct_var(t, 'oii_3729_flux')
         dfluxes[*, i] = err_scale*(sqrt(struct_var(t, 'oii_3726_flux_err')^2+ $
                                         struct_var(t, 'oii_3729_flux_err')^2))
      endif else begin
         fluxes[*, i] = struct_var(t, lines[i]+'_flux')
         dfluxes[*, i] = struct_var(t, lines[i]+'_flux_err')*err_scale
      endelse
   endfor

   ;; Go through and calculate the results

   qs = [0.025, 0.16, 0.5, 0.84, 0.84, 0.975]
   ratio = fltarr(n_elements(t), n_elements(qs))-9999
   r_full = fltarr(n_elements(t), n_mc)
   for i=0L, n_elements(t)-1 do begin
      ;; Skip low-SN cases
      if fluxes[i, i_norm]/dfluxes[i, i_norm] lt min_sn then continue

      t_x = dblarr(n_mc, n_elements(lines))
      for i_line=0L, n_elements(lines)-1 do $
       t_x[*, i_line] = fluxes[i, i_line]+randomn(sss, n_mc)*dfluxes[i, i_line]

      r = alog10((t_x[*, 0]/t_x[*, 1])/(t_x[*, 2]/t_x[*, 3]))
      ratio[i, *] = quantile(r, qs)
      r_full[i, *] = r
   endfor

   return, {summarised: ratio, full: r_full, quantiles: qs}
end
