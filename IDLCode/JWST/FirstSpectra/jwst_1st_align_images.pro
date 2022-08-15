pro jwst_1st_align_images
   ;;
   ;; Align all images to F200W. 
   ;;
   jwst_1st_dirs, dd
   
   IMDIR = dd.root+'Images/'
   OUTDIR = IMDIR+'Aligned/'
   readcol, IMDIR+'alignment-positions.txt', filter, ra, dec, skip=1, $
            format='A,F,F'
 
   i_ref = where(filter eq 'F200W')
   dra = ra-ra[i_ref]
   ddec = dec-dec[i_ref]

   jwst_filters = 'f'+['090', '150', '200', '277', '356', '444']+'w'
   hst_filters = 'f'+['450', '606', '814']+'w'

   jwst_files = IMDIR+'jw02736-o001_t001_nircam_clear-'+jwst_filters+'_i2d.fits'
   hst_files = dd.root+'HST/hlsp_relics_hst_acs-30mas_smacs0723-73_'+ $
               hst_filters+'_v1_drc.fits'

   filters = [hst_filters, jwst_filters]
   imfiles = [hst_files, jwst_files]

   ;;
   ;; Load the reference image.
   ;;
   imref_file = imfiles[where(filters eq 'f200w')]
   imref_o = mrdfits(imref_file, 'SCI', hdr_im_ref_o)
   getrot, hdr_im_ref_o, rot, cdelt
   hrot, imref_o, hdr_im_ref_o, imref, hdr_im_ref, rot, -1, -1, 0
   outfile = OUTDIR+'aligned-to-f200W-f200w.fits'
   sxdelpar, hdr_im_ref, 'XTENSION'
   mwrfits, imref, outfile, hdr_im_ref, /create


   ;;
   ;; Go through each image, read it in and apply a shift to the
   ;; header. Then match to the reference image.
   ;;
   for i=0L, n_elements(filters)-1 do begin

      print, 'Doing filter '+filters[i]


      if filters[i] eq 'f450w' then continue

      if i le 2 then $
       im_orig_o = mrdfits(imfiles[i], 0, hdr_im_orig_o) $
      else $
       im_orig_o = mrdfits(imfiles[i], 'SCI', hdr_im_orig_o) 
      getrot, hdr_im_orig_o, rot, cdelt
      ;; Rotate to get N & E aligned properly
      hrot, im_orig_o, hdr_im_orig_o, im_orig, hdr_im_orig, rot, -1, -1, 0

      ;; Apply shift.
      ra_this = sxpar(hdr_im_orig, 'CRVAL1')
      dec_this = sxpar(hdr_im_orig, 'CRVAL2')
      print, "CRVAL1 before=", ra_this
      print, "CRVAL2 before=", dec_this
      ra_this = ra_this[0]-dra[0]
      dec_this = dec_this[0]+ddec[0]

;      print, 'Shifting by Delta Ra=', dra, format='(A,2X,F12.7)'
;      print, 'Shifting by Delta Ra=', ddec, format='(A,2X,F12.7)'

      sxaddpar, hdr_im_orig, 'CRVAL1', ra_this
      sxaddpar, hdr_im_orig, 'CRVAL2', dec_this
;      mwrfits, im_orig, OUTDIR+'check.fits', hdr_im_orig, /create

;      print, "CRVAL1 after=", sxpar(hdr_im_orig, 'CRVAL1'), format='(A,2X,F12.7)'
;      print, "CRVAL2 after=", sxpar(hdr_im_orig, 'CRVAL2'), format='(A,2X,F12.7)'
      


      

      
      hastrom, im_orig, hdr_im_orig, im, hdr_im, hdr_im_ref
      sxdelpar, hdr_im, 'XTENSION'
      
      outfile = OUTDIR+'aligned-to-f200W-'+filters[i]+'.fits'
      
      mwrfits, im, outfile, hdr_im, /create
   endfor

end
