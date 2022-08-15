pro jwst_1st_aper_photometry_all, debug=debug
   ;;
   ;; Go through all galaxies and put a 1" aperture on the
   ;; centre to get aperture photometry
   ;;

   ROOT = '/data2/jarle/JWST/ERO/SMACS/'
   filters = ['f090w', 'f150w', 'f200w', 'f277w', 'f356w', 'f444w']

   nc = hdf5_to_struct(ROOT+'NIRCam-photometry-NIRSpec.h5')
   s = nc.sources

   ;; Note that the error analysis is likely quite flakey so it is
   ;; recalibrated later
   fluxes = dblarr(n_elements(s), n_elements(filters))-9999.
   dfluxes = dblarr(n_elements(s), n_elements(filters))-9999.

   
   for i=0L, n_elements(s)-1 do begin
      for i_f=0L, n_elements(filters)-1 do begin
         stamp = ROOT+'Stamps/'+strtrim(s[i])+'_'+filters[i_f]+'.fits'

         if not file_test(stamp) then begin
            print, 'No stamp for '+s[i]+' in filter '+filters[i_f]
            continue
         endif
         im = mrdfits(stamp, 0, hdr, /silent)
         
         getrot, hdr, rot, cdelt
         pixscale = abs(cdelt[0])*3600
         dims = size(im, /dimen)
         xc = dims[0]/2.0
         yc = dims[1]/2.0

         ;; We also want some random positions
         x_sky = randomu(sss, 150)*dims[0]*0.8+dims[0]*0.1
         y_sky = randomu(sss, 150)*dims[1]*0.8+dims[1]*0.1

;         stop
         ap_size = 1.0/pixscale
         sky_ap_size = [1.5, 2.5]/pixscale

;         print, 'I will use an aperture scale of ', ap_size

         aper, im, xc, yc, flux, errflux, sky, skyerr, $
               1.0, [ap_size], sky_ap_size, [0, 0], /flux, $
               /silent
;         stop

         sky_value = sky[0]
         ;; And then the sky regions
         aper, im, x_sky, y_sky, flux_sky, errflux_sky, sky, skyerr, $
               1.0, [ap_size], sky_ap_size, [0, 0], /flux, $
               /silent, setskyval=sky_value

;         if i eq 2 then stop

         
         fluxes[i, i_f] = flux
         dfluxes[i, i_f] = errflux
         if keyword_set(debug) then begin
            tvim_true, im
            tvcircle, ap_size/2.0, xc, yc, /data
            tvcircle, sky_ap_size[0]/2.0, xc, yc, /data, color=jb_colour('red')
            tvcircle, sky_ap_size[1]/2.0, xc, yc, /data, color=jb_colour('red')
         endif
      endfor
      
      
   endfor

   phot = {source: '', f090w: 0.0d0, f150w: 0.d0, f200w: 0.0d0, f277w: 0.0d0, $
          f356w: 0.0d0, f444w: 0.0d0}
   phot = {source: s}
   for i_f=0L, n_elements(filters)-1 do begin
      phot = create_struct(phot, filters[i_f], reform(fluxes[*, i_f]), $
                           'd'+filters[i_f], reform(dfluxes[*, i_f]))
   endfor
   phot = struct_of_arr_to_arr_of_struct(phot)

   
   mwrfits, phot, ROOT+'Stamps/AperturePhotometry/photometry_idl.fits', /create
   
   
end
