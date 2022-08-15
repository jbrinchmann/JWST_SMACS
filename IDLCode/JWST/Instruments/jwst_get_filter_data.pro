function jwst_get_filter_data, filter, instrument, object=object
   ;;
   ;; Look up basic filter information

   if strlowcase(instrument) eq 'nircam' then begin
      ;;
      ;; NIRCAM transmission files
      ;;

      NIRCAM_ROOT = '/data2/jarle/JWST/InstrumentData/NIRCAM/'
      ;; We use the mean AB filters total
      FILTER_DIR = NIRCAM_ROOT+'nircam_throughputs/modAB_mean/nrc_plus_ote/'

      fname = FILTER_DIR+filter+'_NRC_and_OTE_ModAB_mean.txt'

      if not file_test(fname) then begin
         print, 'I could not find a file for filter '+filter
         stop
         return, -1
      endif

      rdfloat, fname, l_f, throughput, skip=1

      result = {wavelength: l_f, throughput: throughput, $
                lambdaAA: l_f*1e4, T: throughput}
      
   endif else stop

   
   ;;

;;    NIRSPEC = DICTIONARY(

;;    G140M/F070LP	~1,000



;; 	0.70–1.27
;; G140M/F100LP	0.97–1.84
;; G235M/F170LP	1.66–3.07
;; G395M/F290LP	2.87–5.10
;; G140H/F070LP	~2,700

;; 	0.81–1.27
;; G140H/F100LP	0.97–1.82
;; G235H/F170LP	1.66–3.05
;; G395H/F290LP	2.87–5.14
;; PRISM/CLEAR	~100	0.60-5.30

   if keyword_set(object) then begin
      ;;
      ;; In this case get a filter object out.
      ;;
      result = obj_new('filter', filter, result.lambdaAA, result.throughput, $
                       lambda_units='Angstrom')

   endif
   
   return, result

end
