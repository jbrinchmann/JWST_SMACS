pro jwst_1st_extract_images, width=width, no_align=no_align
   ;;
   ;; Extract images for each object here.
   ;;
   if (n_elements(width) eq 0) then width = 5.0 ; arcsec
   

   SPECDIR = '/data2/jarle/JWST/ERO/SMACS/Spectra/'
   OUTDIR = '/data2/jarle/JWST/ERO/SMACS/Stamps/'
   IMDIR = '/Volumes/AstroData/JWST/ERO/SMACS/jw02736/L3/t/'
   HSTDIR = '/data2/jarle/JWST/ERO/SMACS/HST/'

   ;; Read in the locations for NIRCam
   readcol, SPECDIR+'../matches-with-comments.txt', specid, ra_nirspec, dec_nirspec, $
            racorr, deccorr, idNIRCAM_f200w, ra_nc, dec_nc, sep, $
            format='A,F,F,F,F,L,F,F,F', /silent

   t = {specid: specid, ra_nc: ra_nc, dec_nc: dec_nc, dist_to_cat: sep}
   t = struct_of_arr_to_arr_of_struct(t)
   mwrfits, t, SPECDIR+'../matches-with-comments.fits', /create
   
   ;; And the Matched RELICS catalogue
   t_relics = mrdfits(HSTDIR+'../matched-to-relics.fits', 1)

   ;; We use the NIRCAM locations below but we also need the HST
   ;; coordinates 
   ra_hst = t_relics.ra_wf3
   dec_hst = t_relics.dec_wf3

   ;; In some cases the match is non-existent
   for i=0L, n_elements(ra_hst)-1 do begin
      if (ra_hst[i] lt 0 and t_relics[i].ra_acs gt 0) then begin
         ra_hst[i] = t_relics[i].ra_acs
         dec_hst[i] = t_relics[i].dec_acs
      endif
   endfor
   
   ;;; Add some extra coordinates. We apply the average offset between
   ;;; NIRCam and ACS of 1" in ra.

   ii = where(ra_hst lt 0, n_ii)
   if (n_ii gt 0) then begin
      ra_hst[ii] = ra_nc[ii]-1.0/3600.
      dec_hst[ii] = dec_nc[ii]
   endif
      
   ;;
   ;; I want all stamps to be aligned to a reference image for ease of
   ;; later display. I choose this to be the f150w NIRCam image for no
   ;; particular reason
   ;;
   ;; However, for photometry I also want the original so that is
   ;; stored as well.
   ;;
   imref_file = IMDIR+'jw02736-o001_t001_nircam_clear-f150w_i2d.fits'
   imref = mrdfits(imref_file, 'SCI', hdr_im_ref)
   
   

   ;; Focus on NIRCAM first and then HST
   filters = ['f090w', 'f150w', 'f200w', 'f277w', $
              'f356w', 'f444w', 'HST_f435w', 'HST_f606w', 'HST_f814w']

   for i=0L, n_elements(filters)-1 do begin

      print, 'Doing filter '+filters[i]

      if strmid(filters[i], 0, 3) eq 'HST' then begin
         tmpfilter = strmid(filters[i], 4)
         imfile = HSTDIR+'hlsp_relics_hst_acs-30mas_smacs0723-73_'+tmpfilter+'_v1_drc.fits'
         im_orig = mrdfits(imfile, 0, hdr_im_orig)

         if keyword_set(no_align) then begin
            im = im_orig
            hdr_im = hdr_im_orig
         endif else begin
            hastrom, im_orig, hdr_im_orig, im, hdr_im, hdr_im_ref
         endelse
         ra = ra_hst
         dec = dec_hst
         is_hst = 1
         
      endif else begin
         root = IMDIR+'jw02736-o001_t001_nircam_clear-'+filters[i]
         
         ;; We want image and segmentation map for the area
         segfile = root+'_segm.fits'
         imfile = root+'_i2d.fits'
         catfile = root+'_cat.ecsv.fits'

         ;; Actual loading
         im_orig = mrdfits(imfile, 'SCI', hdr_im_orig)
         seg_orig = mrdfits(segfile, 'SCI', hdr_seg_orig)
         cat = mrdfits(catfile, 1)

         ;; Alignment
         if keyword_set(no_align) then begin
            im = im_orig
            hdr_im = hdr_im_orig
         endif else begin
            hastrom, im_orig, hdr_im_orig, im, hdr_im, hdr_im_ref
         endelse
         
         ra = ra_nc
         dec = dec_nc

         is_hst = 0
      endelse

      
      ;; Calculate x & y positions
      adxy, hdr_im, ra, dec, xpos, ypos
      
      ;; Get platescale
      getrot, hdr_im, rot, cdelt

      width_pixel = width/(abs(cdelt[0])*3600)
      dx = ceil(width_pixel/2.) ; At least the width
      dy = dx

      dims = size(im, /dimen)
      nx = dims[0]
      ny = dims[1]
      
;      stop
      ;; Loop over all objects and extract stamps
      for i_obj=0L, n_elements(ra)-1 do begin

         if is_hst eq 0 then begin
            if idNIRCAM_f200w[i_obj] eq -99 then begin
               print, 'Skipping '+specid[i_obj]+' due to lack of association'
               continue
            endif
         endif else begin
            if ra[i_obj] lt 0 then begin
               print, 'Skipping '+specid[i_obj]+' due to lack of association'

               continue
            endif
         endelse
         
         x0 = xpos[i_obj]-dx
         x1 = xpos[i_obj]+dx
         y0 = ypos[i_obj]-dy
         y1 = ypos[i_obj]+dy

         adjust = 0 ; No adjustment needed by default
         if (x0 lt 0) or  (x1 gt nx-1) $
          or (y0 lt 0) or (y1 gt ny-1) then begin
            ;;
            ;; Adjust edges
            ;;
            adjust = 1
;            stop
            ;; This a ctually only happens for y1 for these images
            y1 = ny-1
         endif
         
         hextract, im, hdr_im, stamp, hdr_stamp, x0, x1, y0, y1

         if adjust then begin
            stampfull = fltarr(dx*2+1, dy*2+1)
            stampfull[*,*] = 0 ; Fill all with zeros
            ;; only happens on the top edge.
            yorig = ypos[i_obj]+dy

            ;; Insert the data
            dims = size(stamp, /dimen)
            stampfull[*, 0:dims[1]-1] = stamp
            
            sxaddpar, hdr_stamp, 'NAXIS2', dy*2
            
            ;; Do we need to adjust the header beyond nAXIS?
            stamp = stampfull

            print, 'Adjusted '+specid[i_obj]+" for filter "+filters[i]
         endif

         if keyword_set(no_align) then $
          suffix = '-orig' $
         else $
          suffix = ''
         
         outfile = OUTDIR+specid[i_obj]+'_'+filters[i]+suffix+'.fits'
         sxdelpar, hdr_stamp, 'XTENSION'
         mwrfits, stamp, outfile, hdr_stamp, /create

;         if adjust then stop
      endfor

   endfor
   

end



