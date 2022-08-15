
function jwst_1st_fit_gaussian, sp, l_line, pf=pf, adjust=adjust
   ;;
   ;; Fit a single Gaussian on a spectrum. The user must also give a
   ;; platefit-like structure with information. 
   ;;

   x = sp.wavelength
   y = sp.flux
   dy = sp.dflux
   
   x_rest = x/(1.0+pf.z)

   i_close = where(abs(x_rest-l_line) lt 150 and finite(y) eq 1, n_close)
   x = x[i_close]
   y = y[i_close]
   dy = dy[i_close]

   if keyword_set(adjust) then begin
      dl = mkarr(-10, 10, 0.01)
      fluxes = dblarr(n_elements(dl))
      for i=0L, n_elements(dl)-1 do begin
         g = fit_one_gaussian_mp(x, y, a, /return_flux, errors=dy, $
                           center = l_line+dl[i], covar = da)
         fluxes[i] = a[3]
      endfor
      stop
   endif else begin
      g = fit_one_gaussian_mp(x, y, a, /return_flux, errors=dy, $
                              center=l_line, covar=da)
;      stop
   endelse
   
end
