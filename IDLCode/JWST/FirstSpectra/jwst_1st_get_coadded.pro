function jwst_1st_get_coadded, number, filter=filter, blank=blank, $
                               keep_bad=keep_bad, d1=d1, d2=d2

   if n_elements(blank) eq 0 then blank = 1

   
   DIR = '/data2/jarle/JWST/ERO/SMACS/Spectra/'
   
   propid = string(format='(I5.5)', 2736)
   objid = string(format='(I5.5)', number)
   
   rootnames = 'jw'+propid+'-o00'+['7', '8']+'_s'+objid
   ;; Get the data

   d1 = jwst_load_nirspec_data(DIR, rootnames[0], blank=blank, $
                               keep_bad=keep_bad)

   
   d2 = jwst_load_nirspec_data(DIR, rootnames[1], blank=blank, $
                               keep_bad=keep_bad)

;   stop
   if number eq 4590 then begin
      ;;
      ;; Following Curti et al we only use observation 8 - they also
      ;; used part of 7 but as I am not re-reducing the data that is
      ;; not possible
      ;;
      d = d2['combined']

   endif else if number eq 5144 then begin
      ;;
      ;; Here the line ratios in the 008 observations are wonky.
      ;;
;      stop
      d = d1['combined']
   endif else if number eq 4580 then begin
      ;;
      ;; Visit 007 has a major problem near Ha.
      ;;
;      stop
      d = d2['combined']
   endif else begin
      d = jwst_combine_nirspec_exposures(d1, d2, factor=10., filter=filter)
   endelse
   

   return, d
end
