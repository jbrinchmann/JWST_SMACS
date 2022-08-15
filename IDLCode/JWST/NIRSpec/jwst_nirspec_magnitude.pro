function simple_magnitude, l, f, filtername, filter=filter
   ;;
   ;; Given a wavelength and a flux, estimate the normalised flux sum
   ;;
   ;; I am assuming input in F_lambda

   ;; Do all calculations in AAngstrom 
   if max(l) lt 1000 then $
    l_AA = l*1e4 $
   else $
    l_AA = l


   ;;
   ;; Get the filter - we will set everything outside 1e-4 throughput
   ;; to zero. This is not always ok but I'll use it here. 
   ;;
   if N_elements(filter) eq 0 then $
    filter = jwst_get_filter_data(filtername, 'NIRCAM')
   useful_filter = where(filter.T gt 1e-4*max(filter.T), n_useful)
   limits_filter = minmax(filter.lambdaAA)

   ;;
   ;; Interpolate onto the wavelength grid of the spectrum
   ;;
   filter_on_L = interpol(filter.T, filter.wavelength, l_AA/1e4)
   outside = where(l_AA lt limits_filter[0] or l_AA gt limits_filter[1], n_outside)
   if (n_outside gt 0) then filter_on_L[outside] = 0.0

   ;;
   ;; I will use a sum approximation because it makes it easier to
   ;; deal with missing data in the middle.
   ;;
   missing = where(finite(f) eq 0 or f eq 0.0, n_missing, complement=not_missing)
   
   ;; Scale for the flux count - assuming uniform source - given the
   ;; large gaps this is probably reasonable.
   scale = total(filter_on_L)/total(filter_on_L[not_missing])

   ;; Summed normalisation integral:
   dl = l-shift(l, 1)
   dl[0] = dl[1]
   cspeed = 2.99792458d18       ; AA/s
   norm = totaL(filter_on_L*dl/l_AA^2)*cspeed
   int = total(filter_on_l[not_missing]*f[not_missing]*dl[not_missing])*scale

   return, {flux: int/norm, mag: -2.5*alog10(int/norm), norm: norm}
   

end


function jwst_nirspec_magnitude, sp, filternames
   ;;
   ;; Given a filter structure, calculate the magnitude through the
   ;; given filter. The magnitude is in AB.
   ;;

   if n_elements(filternames) eq 1 then $
    filters = [filternames] $
   else $
    filters = filternames

 
   ;; And make a spectrum object and a variance one
   if max(sp.wavelength) lt 1000 then begin
      l = sp.wavelength*1e4
   endif else begin
      l = sp.wavelength
   endelse


   ;; Set zeros to NaNs
   missing = where(sp.flux eq 0.0, n_missing)
   if (n_missing eq 0) then stop
   sp.flux[missing] = !values.d_nan
   sp.dflux[missing] = !values.d_nan

   
   ;; For this we would like to fill the NaNs or zeros.
;   spec = obj_new('spectrum', l, sp.flux)
;   var_spec = obj_new('spectrum', l, sp.dflux^2)



   
;   stop

   N_f = n_elements(filters)
   mags = fltarr(N_f)
   var_mags = mags
   flux = mags
   dflux = mags
   filternorm = mags
   for i=0L, N_f-1 do begin
      filter = jwst_get_filter_data(filters[i], 'NIRCAM')
;      mags[i] = spec->magnitude(f, 0.0, /ab, /quiet, flux=tmp_flux, dflux=sp.dflux)
;      var_mags[i] = var_spec->magnitude(f, 0.0, /ab, /quiet, flux=tmp_var_flux)
;      flux[i] = tmp_flux
;      var_flux[i] = tmp_var_flux

      tmp = simple_magnitude(sp.wavelength, sp.flux, filters[i], filter=filter)

      ;; Let's try to Monte Carlo it all..
      n_Mc = 5000
      xr = dblarr(n_mc)
      for i_mc=0L, n_MC-1 do begin
         l = sp.wavelength
         f = sp.flux+randomn(ss, n_elements(sp.flux))*sp.dflux
         tmp = simple_magnitude(l, f, filters[i], filter=filter)
         xr[i_mc] = tmp.flux
      endfor
      v= quantile(xr, [0.16, 0.5, 0.84])
      flux[i] = v[1]
      dflux[i] = 0.5*(v[2]-v[1])
      filternorm[i] = tmp.norm
      
   endfor


   return, {mag: -2.5*alog10(flux), dmag: 1.08574*dflux/flux, flux: flux, $
            dflux: dflux, filternorm: filternorm}
   
end
