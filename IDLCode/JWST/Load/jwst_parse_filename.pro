function jwst_parse_filename, full_fname, level=level
   ;;
   ;; Parse a Level X filename. At the moment only level 3 is
   ;; implemented and anything else will give an error.

   if n_elements(level) eq 0 then $
    level = 3


   if level ne 3 then begin
      print, 'Only level 3 products are supported at the moment'
      return, -1
   endif


   ;; Remove the directory.
   fname = file_basename(full_fname)

   
   ;;
   ;; The naming convention is taken from
   ;;    https://jwst-pipeline.readthedocs.io/en/latest/jwst/data_products/file_naming.html
   ;;
   ;;
   ;;     jw<ppppp>-<AC_ID>_[<”t”TargID | “s”SourceID>](-<”epoch”X>)_<instr>_<optElements>(-<subarray>)_<prodType>(-<ACT_ID>).fits
   
   ;; where
   ;;         ppppp: Program ID number
   ;;         AC_ID: Association candidate ID
   ;;         TargID: 3-digit Target ID (either TargID or SourceID must be present)
   ;;         SourceID: 5-digit Source ID
   ;;         epochX: The text “epoch” followed by a single digit epoch number (optional)
   ;;         instr: Science instrument name (e.g. ‘nircam’, ‘miri’)
   ;;         optElements: A single or hyphen-separated list of optical elements (e.g. filter, grating)
   ;;         subarray: Subarray name (optional)
   ;;         prodType: Product type identifier (e.g. ‘i2d’, ‘s3d’, ‘x1d’)
   ;;         ACT_ID: 2-digit activity ID (optional)
   ;; An example Stage 3 product FITS file name is:
   ;;     jw87600-a3001_t001_niriss_f480m-nrm_amiavg.fits


   

   ;; Regexp for parsing
  re = '^jw([0-9]+)-([a-zA-Z0-9]+)_(s|t)([0-9]+)(-epoch[0-9])?_([a-zA-Z0-9]+)_([-a-zA-Z0-9]+)_([0-9A-Za-z]+)(-[0-9]+)?.fits'

   m = stregex(fname, re, /subexpr, /extract)

   if strlen(m[0]) eq 0 then begin
      print, 'The filename does not seem to conform to JWST naming conventions'
;      stop
      return, -1
   endif

   ;; We now need to process this to create a handy return information
   progID = long(m[1])
   progIDstr = m[1]
   ac_id = m[2]
   if m[3] eq 's' then begin
      sourceID = long(m[4])
;      stop
      targetID = -999
      single_target = 0
   endif else begin
      targetID = long(m[4])
      stop
      sourceID = -999
      single_target = 1
   endelse

   if strlen(m[5]) gt 0 then begin
      epoch = long(strmid(m[5], 6, 1))
      multi_epoch = 1
   endif else begin
      multi_epoch = 0
      epoch = 0
   endelse

   instrument = m[6]
   
   ;; The separation between subarray and optical elements is not made
   ;; in the regexp above so needs to be dealt with here.
   optical_elements = m[7]
   list = strsplit(optical_elements, '-', /extract)

   filter = 'None'
   grating = 'None'
   subarray = 'None'
   for i=0L, n_elements(list)-1 do begin
      first_char = strmid(strlowcase(list[i]), 0, 1)
      if first_char eq 'f' or first_char eq 'c' then begin
         filter = list[i]
      end else if first_char eq 'g' or first_char eq 'p' then begin
         grating = list[i]
      endif else $
       subarray = list[i] ;;; ?? Untested!
   endfor

   optical_elements = filter+'-'+grating
   
   product_type = m[8]

   info = {progID: progID, progIDstr: progIDstr, ac_id: ac_id, sourceID: sourceID, $
           targetID: targetID, single_target: single_target, epoch: epoch, $
           multi_epoch: multi_epoch, instrument: instrument, $
           optical_elements: optical_elements, subarray: subarray, $
           product_type: product_type, filter: filter, grating: grating}

   return, info


   
end
