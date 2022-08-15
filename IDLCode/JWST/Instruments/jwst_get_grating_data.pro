function jwst_get_grating_data, grating, instrument, filter=filter, info=info


   NIRSPEC = {PRISMCLEAR: {lmin: 0.6, lmax: 5.3, Rmin: 30, Rmax: 330}, $
              G140MF070LP: {lmin: 0.7, lmax: 1.3, Rmin: 500, Rmax: 890}, $ 
              G140MF100LP: {lmin: 1.0, lmax: 1.9, Rmin: 700, Rmax: 1340}, $ 
              G235MF170LP: {lmin: 1.7, lmax: 3.2, Rmin: 720, Rmax: 1340}, $ 
              G395MF290LP: {lmin: 2.9, lmax: 5.2, Rmin: 730, Rmax: 1315}, $ 
              G140HF070LP: {lmin: 0.7, lmax: 1.3, Rmin: 1320, Rmax: 2395}, $ 
              G140HF100LP: {lmin: 1.0, lmax: 1.9, Rmin: 1850, Rmax: 3675}, $ 
              G235HF170LP: {lmin: 1.7, lmax: 3.2, Rmin: 1910, Rmax: 3690}, $ 
              G395HF290LP: {lmin: 2.9, lmax: 5.2, Rmin: 1930, Rmax: 3615}}

   instruments = {NIRSPEC: NIRSPEC}


   if (n_elements(info) gt 0) then begin
      grating = info.grating
      instrument = info.instrument
      filter = info.filter
   endif
   

   table = struct_var(instruments, instrument)

   if strlowcase(instrument)  eq 'nirspec' then begin
      if n_elements(filter) eq 0 then begin
         print, 'For NIRSpec you need to provide both grating and filter!'
         return, -1
      endif
      key = grating+filter
   endif else key = grating

   ;; 
   ;; Get grating curve if possible
   ;;
   DATADIR_ROOT = '/data2/jarle/JWST/InstrumentData/'
   DATADIR = DATADIR_ROOT+strupcase(instrument)+'/'
   case strlowcase(instrument) of
      'nirspec': begin
         fname = 'jwst_nirspec_'+grating+'_disp.fits'
         curve = mrdfits(DATADIR+fname, 1, /silent)
      end
      else: curve = -1
   endcase

   res = struct_var(table, key)
   res = create_struct(res, 'CURVE', curve)
   return, res
end
