function jwst_1st_load, name, min_restwl=min_restwl
   ;;
   ;; Load one of A, B C or D.
   ;;

   ;; By default keep all
   if (n_elements(min_restwl) eq 0) then min_restwl = 0.0
   
   jwst_1st_dirs, dd

   redshifts = {A: 7.6601913, B: 2.7411191, C: 8.4949753, $
                D: 5.2755742}


   z = struct_var(redshifts, name)
   
   fname = dd.root+'Spectra_Target'+name+'.csv'
   r = r_read_table(fname, /header, sep=',')

   ;; Reformulate
   l = 1e4*r.wavelength_in_microns
   f = r.relative_brightness*1e-17
   l_mu = l/1e4

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
      if f[i] gt 3*med then $
       flag[i-(dw+3):i+(dw+3)] = 1
   endfor
   high = where(flag eq 1, n_high)

   
   
   if (n_high gt 0) then $
    sig[high] = med

   ;; Finally, make the wavelength array unique
   si = sort(l)
   ui = uniq(l[si])

   f = f[si[ui]]
   l = l[si[ui]]
   l_mu = l_mu[si[ui]]
   sig = sig[si[ui]]


   ;; There are also big gaps between the bands. These are best
   ;; handled with NaNs as that will play ok with platefit.
   diff=l-shift(l,1)
   big_jumps = where(diff gt 50, n_big)


   lnew = l
   fnew = f
   signew = sig
   if (n_big gt 0) then begin
      offset = 0L
      for igap=0L, n_big-1 do begin
         ;; These two bridge the gap
         i_left = big_jumps[igap]-1
         i_right = big_jumps[igap]

         ;; Create a flux array full of NaNs with wavelength steps
         ;; equal to the media of the 10 bins before. I do this to
         ;; avoid accidentally lots of lots of small steps
         dlambda = median(diff[big_jumps[igap]-11:big_jumps[igap]-1])
         infill_l = mkarr(l[i_left]+dlambda, l[i_right]-dlambda, dlambda)
         infill_f = infill_l*0.0+!values.f_nan

         lnew = [lnew[0:i_left+offset], infill_l, lnew[i_right+offset:*]]
         fnew = [fnew[0:i_left+offset], infill_f, fnew[i_right+offset:*]]
         signew = [signew[0:i_left+offset], infill_f, signew[i_right+offset:*]]
         offset = offset + n_elements(infill_l)
      endfor
   endif

   l_mu = double(l/1e4)
   f = double(fnew)
   l = double(lnew)
   sig = double(signew)
   restwl = l/(1+z)


   keep = where(restwl gt min_restwl)
   l = l[keep]
   l_mu = l_mu[keep]
   f = f[keep]
   sig = sig[keep]
   restwl = restwl[keep]

   ;; ui = uniq(l)
   ;; l = l[ui]
   ;; f = f[ui]
   ;; l_mu = l_mu[ui]
   ;; sig = sig[ui]

;   stop
   
   
   sp = {flux: f, wave: l, wave_mu: l_mu, dflux: sig, $
        restwl: restwl, z: z}

   return, sp
end
