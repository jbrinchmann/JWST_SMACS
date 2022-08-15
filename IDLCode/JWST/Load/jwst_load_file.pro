;;
;; The JWST files appear to have an organised format so they can
;; be loaded straightforwardly letting the code here dispatch
;; to specific loaders depending on the filename
;;

function jwst_load_x1d, fname
   ;;
   ;; Load x1d spectral files
   ;;

   ;; This is stored in a table
   t = mrdfits(fname, 'EXTRACT1D', hdr, /silent)

   dum = mrdfits(fname, 0, primary_hdr, /silent)

   tn = tag_names(t)
   tout = {header: hdr, primary_header: primary_hdr}
   for i=0L, n_elements(tn)-1 do begin
      x = struct_var(t, tn[i])
      tout = create_struct(tout, tn[i], x)
   endfor

   type = sxpar(hdr, 'SRCTYPE')
   RA = sxpar(hdr, 'SRCRA')
   DEC = sxpar(hdr, 'SRCDEC')
   id = sxpar(hdr, 'SOURCEID')
   slit_id = sxpar(hdr, 'SLITID')
   name = sxpar(hdr, 'SRCNAME')
   slit_RA = sxpar(hdr, 'SLIT_RA')
   slit_DEC = sxpar(hdr, 'SLIT_DEC')

   tout = create_struct(tout, 'type', type, 'ra', ra, 'dec', dec, 'id', id, $
                        'slit_id', slit_id, 'name', name, 'slit_ra', slit_ra, $
                        'slit_dec', slit_dec)
   
   return, tout
end

function jwst_load_s2d, fname
   ;;
   ;; Load 2D spectral files
   ;;

   ;; This is a multi-extension FITS file and we
   ;; load only flux and uncertainty array
   
   flux = mrdfits(fname, 'SCI', hdr, /silent)
   dflux = mrdfits(fname, 'ERR', /silent)

   t_exp = sxpar(hdr, 'XPOSURE')
   unit = sxpar(hdr, 'BUNIT')
   RA = sxpar(hdr, 'SRCRA')
   DEC = sxpar(hdr, 'SRCDEC')
   id = sxpar(hdr, 'SOURCEID')
   name = sxpar(hdr, 'SRCNAME')
   
   return, {flux: flux, dflux: dflux, header: hdr, exptime: t_exp, $
            unit: unit, ra: ra, dec: dec, id: id, name: name}
end



function jwst_load_file, fname, status=status
   ;;
   ;; Dispatcher
   ;;

   tail = strmid(fname, strlen(fname)-8, 3)
   status = 1                   ; Assume ok
   case tail of
      'x1d': res = jwst_load_x1d(fname)
      'cal': res = jwst_load_cal(fname)
      's2d': res = jwst_load_s2d(fname)
      'crf': res = jwst_load_crf(fname)
      else: begin
         print, 'Unknown (or unimplemented) extension!'
         res = -1
         status = -1
      end
   endcase

   return, res
end
