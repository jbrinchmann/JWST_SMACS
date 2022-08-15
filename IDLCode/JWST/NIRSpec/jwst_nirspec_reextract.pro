function jwst_nirspec_reextract, d, centre, narrow_sum=narrow_sum, $
                                 ap_background=ap_background, bad=bad
   ;;
   ;; Given a spectrum structure, re-extract a 1D spectrum and return
   ;; this updated spectrum. 
   ;;

   ;;
   ;; A couple of different ways exist:
   ;;
   ;;   NARROW_SUM: This sums simply over a rectangle centred on
   ;;               CENTRE with a width equal to NARROW_SUM
   ;; 
   ;;   Gaussian: this fits a Gaussian to each bin in x (is this feasible?) 
   ;;

   if n_elements(narrow_sum) gt 0 then begin
      im_orig = d.twod.flux
      dim_orig = d.twod.dflux
      var_im_orig = dim_orig^2
      
      width = narrow_sum
      iy_min = floor(centre-width/2)
      iy_max = ceil(centre+width/2)
      
      im = im_orig[*, iy_min:iy_max]
      dim = dim_orig[*, iy_min:iy_max]
      var_im = dim^2
      
      dims = size(im, /dimen)
      nx = dims[0]
      ny = dims[1]

      dims_orig = size(im_orig, /dimen)
      nx_orig = dims[0]
      ny_orig = dims[1]
      
      ;; If it is requested to do a adjacent background subtraction
      ;; (this increases noise of course)
      if keyword_set(ap_background) then begin
         ii_max = iy_max+width
         ii_min = iy_min-width
         if (ii_min lt 0) then ii_min = 0
         if (ii_max gt dims_orig[1]-1) then ii_max = dims_orig[1]-1
         
         side1 = im_orig[*,ii_min:iy_min-1]
         side2 = im_orig[*,iy_max+1:ii_max]
         side = 0.5*(side1+side2)
         dims = size(side, /dimen)
         sky_spec = total(side, 2)/dims[1]

         v_side1 = var_im_orig[*,ii_min:iy_min-1]
         v_side2 = var_im_orig[*,iy_max+1:ii_max]
         v_side = 0.5*(v_side1+v_side2)
         var_sky_spec = total(v_side, 2)/dims[1]

         
         for i=0L, ny-1 do begin
            im[*, i] = im[*, i]-sky_spec
            var_im[*,i] = var_im[*, i]+var_sky_spec
         endfor

      endif
      
      summed = sigma_clip_sum(im, var_im, clipsig=3.0, dim=1, varsum=var_summed)

   endif

   if (n_elements(bad) gt 0) then begin
      summed[bad.bad] = !values.f_nan
      var_summed[bad.bad] = !values.f_nan
   endif
   dout = d
   dout.oned.flux = summed
   dout.oned.dflux = sqrt(var_summed)
   dout.oned.flux_error = sqrt(var_summed)

   
   return, dout
   
end
