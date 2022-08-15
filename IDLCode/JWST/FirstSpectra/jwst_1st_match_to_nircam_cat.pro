function _get_cat, filter
   ROOT = '/data2/jarle/JWST/ERO/SMACS/'
   NIRCAM_DIR = ROOT+'Images/'

   fname = 'jw02736-o001_t001_nircam_clear-'+strlowcase(filter)+'_cat.ecsv.fits'

   return, mrdfits(NIRCAM_DIR+fname, 1)
   
end

pro jwst_1st_match_to_nircam_cat, min_dist=min_dist
   ;;
   ;; Match up the NIRSpec to NIRCam catalogues.
   ;;

   if n_elements(min_dist) eq 0 then min_dist = 0.1
   
   ROOT = '/data2/jarle/JWST/ERO/SMACS/'
   NIRCAM_DIR = ROOT+'Images/'

   ;; Get the match file
   m = mrdfits(ROOT+'match-ids.fits', 1)

   ;; Then go through the filters
   filters = ['f090w', 'f150w', 'f200w', 'f277w', 'f356w', 'f444w']

   ;; First get the reference positions (from F200W)
   ra_ref = m.nircam_ra
   dec_ref = m.nircam_dec

   ;; Columns to get from the photometry file - the filter name will
   ;; be appended
   to_get = ['aper50_flux', 'aper50_flux_err', 'aper_total_flux', $
             'aper_total_flux_err', 'aper50_abmag', $
             'aper50_abmag_err', 'aper_total_abmag', $
             'aper_total_abmag_err', 'ellipticity', $
             'orientation', 'semimajor_sigma', 'semiminor_sigma']
   

   results = {filters: filters, sources: m.specid}
   for i=0L, n_elements(filters)-1 do begin
      t = _get_cat(filters[i])


      ;; Find closest match.
      i_match = lonarr(n_elements(ra_ref))
      
      for i_obj=0L, n_elements(ra_ref)-1 do begin
         if m[i_obj].NIRCAM_F200W_ID lt 0 then begin
            i_match[i_obj] = -9
            continue
         endif
         
         gcirc, 2, ra_ref[i_obj], dec_ref[i_obj], t.sky_centroid_ra, $
                t.sky_centroid_dec, dis
         close = min(dis, i_close)

         if (close gt min_dist) then begin
            print, 'No match in filter '+filters[i]+' for '+m[i_obj].specid
            i_match[i_obj] = -9
         endif else begin
            ;; Ok, a match!
            i_match[i_obj] = i_close
         endelse
      endfor

      ;; Ok, time to assemble this.
      for j=0L, n_elements(to_get)-1 do begin
         x = struct_var(t, to_get[j])
         xout = fltarr(n_elements(ra_ref))-9999.
         use = where(i_match ge 0)
         xout[use] = x[i_match[use]]

         key = to_get[j]+'_'+filters[i]
         results = create_struct(results, key, xout)
      endfor
   endfor
   
   struct_to_hdf5, results, ROOT+'NIRCam-photometry-NIRSpec.h5'

   stop

end
