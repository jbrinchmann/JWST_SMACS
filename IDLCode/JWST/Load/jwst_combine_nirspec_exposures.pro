function jwst_combine_nirspec_exposures, d1, d2, filter=filter, factor=factor
   ;;
   ;; Combine two NIRSpec exposures. This works on 1D data only
   ;; 

   if (n_elements(factor) eq 0) then factor = 10.0
   if (n_elements(filter) eq 0) then $
    filter = 'combined'
   
   ;;
   ;;
   sp1 = d1[filter].oned
   sp2 = d2[filter].oned

   sp1_2d = d1[filter].twod
   sp2_2d = d2[filter].twod

   
   ;;
   ;; I'll keep the d1 as reference
   ;;
   lout = sp1.wavelength
   spout = sp1
   spout_2d = sp1_2d

   f1 = sp1.flux
   df1 = sp1.flux_error
   f2 = interpol(sp2.flux, sp2.wavelength, lout)
   df2 = interpol(sp2.flux_error, sp2.wavelength, lout)

   ;; The 2D image also needs t be interpolated on.
   dims = size(sp1_2d.flux, /dimen)
   im = sp1_2d.flux*0.0
   for i=0L, dims[1]-1 do begin
      im[*, i] = interpol(sp2_2d.flux[*, i], sp2.wavelength, lout)   
   endfor
   
   
   ;; Go through and check for big outliers
   fout = fltarr(n_elements(lout))
   dfout = fout
   for i=0L, n_elements(fout)-1 do begin
      if not finite(f1[i]) and finite(f2[i]) then begin
         fout[i] = f2[i]
         dfout[i] = df2[i]
      endif else if finite(f1[i]) and not finite(f2[i]) then begin
         fout[i] = f1[i]
         dfout[i] = df1[i]
      endif else if f2[i] gt 0 and f1[i] gt 0  then begin
         ;; Both positive so check for big outliers
         if f1[i] gt f2[i]+factor*df2[i] then begin
            ;; Chuck out f1 
            fout[i] = f2[i]
            dfout[i] = df2[i]
         endif else if f2[i] gt f1[i]+factor*df1[i] then begin
            ;; Chuck out f2
            fout[i] = f1[i]
            dfout[i] = df1[i]
         endif else begin
            fout[i] = 0.5*(f1[i]+f2[i])
            dfout[i] = 0.5*sqrt((df1[i]^2+df2[i]^2))
         endelse         
      endif else begin
         ;; At least one negative flux
         fout[i] = 0.5*(f1[i]+f2[i])
         dfout[i] = 0.5*sqrt((df1[i]^2+df2[i]^2))
      endelse

      ;; Now for the images.
      
   endfor

   ;; Then mimick the usual format.
   sp1.flux = fout
   sp1.flux_error = dfout
   if tag_exist(sp1, 'DFLUX') then $
    sp1.dflux = dfout $
   else $
    sp1 = create_struct(sp1, 'dflux', dfout)

   out = {oned: sp1, twod: d1[filter].twod}
           
   
   
   return, out
   
   
end
