pro jwst_1st_convolve_stamps_3042
   ;;
   ;;  Just F090W, F150W and F200W are convolved.
   ;;

   jwst_1st_dirs, dd

   PSFDIR = dd.root+'Images/PSFs/'

   filters = ['f090w', 'f150w', 'f200w']
   images = dd.stampdir+'2736_3042_'+filters+'.fits'
   convolved_images = dd.stampdir+'convolved_2736_3042_'+filters+'.fits'

   for i=0L, n_elements(filters)-1 do begin
      kernel_file = PSFDIR+'kernel_'+strupcase(filters[i])+'_to_F200W.fits'
      kernel = mrdfits(kernel_file)
      im = mrdfits(images[i], 0, hdr)

      dims = size(im, /dimen)
      hcongrid, im, hdr, newim, newhd, outsize=dims*4
;      im = congrid(im, dims[0]*4, dims[1]*4)
      
      imc = convolve(newim, kernel)

      mwrfits, imc, convolved_images[i], newhd, /create
      
   endfor
   

end
