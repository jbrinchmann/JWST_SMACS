function empirical_sigma, l_in, f_in, dw=dw

   if (n_elements(dw) eq 0) then dw = 3

   ok = where(finite(f_in) eq 1)
   f = f_in[ok]
   l = l_in[ok]

   
   ;; Estimate a noise spectrum
   df=f-smooth(f, 10)
   sig = df*0.0
   
   dw = 3
   n_l = n_elements(l)
   for i=0L, n_elements(l)-1 do $
    sig[i]=stdev(df[((i-dw)>0):((i+dw)<(n_l-1))])   

   
   sig_bak = sig
   ;; Check very strong outliers.
   med = median(sig)
   flag = bytarr(n_elements(f))
   for i=0L, n_elements(f)-1 do begin
      if f[i] gt 3*med then begin
         i_low = (i-(dw+3))>0
         i_high = (i+(dw+3)) < (n_elements(f)-1)
        
         flag[i_low:i_high] = 1
      endif
   endfor
   high = where(flag eq 1, n_high)
   
   if (n_high gt 0) then $
    sig[high] = med

   sigout = dblarr(n_elements(f_in))
   sigout[ok] = sig

   return, sigout
end
