
pro jwst_1st_feii_to_snr

   ;; Those objects showing Fe II in emission
   todo = ['03042', '09239', '09483']
   z = [1.9934264, 2.4624202, 1.1615779]
   
   Lsun = 3.826d33              ; erg

   for i=0L, n_elements(todo)-1 do begin
      g = jwst_1st_load_one_gaussfit('02736_'+todo[i])

      ii = where(g.line eq 'FEII_1P257', n_ii)

      
      ;; Flux in erg/s/cm^2
      if (g.flux[ii] gt 1e-15) then $
       scale = 1d-17 $
      else scale = 1.0
         
      feii_1p257 = g.flux[ii]*scale
      dfeii_1p257 = g.dflux[ii]*scale

      
      ;; Convert to Fe II 1.644 - using a fixed 1.359 (Alonso-Herrero
      ;; et al 1997) ratio
      feii_1p644 = feii_1p257/1.359
      dfeii_1p644 = feii_1p257/1.359

      ;; Convert to luminosity in erg/s
      L_1p644 = L_from_flux(feii_1p644, z[i], /ergs)
      dL_1p644 = L_from_flux(dfeii_1p644, z[i], /ergs)
      
      ;; Convert to W.
      L_1p644_W = L_1p644/1e7
      dL_1p644_W = dL_1p644/1e7

      ;; The relation proposed by Alonso-Herrero et al (2003)
      SNr = 0.8*L_1p644_W*1d-34

      print, 'SNR for '+todo[i]+' = ', SNr

      ;; Alternatively - how many SNRs?
      print, "N(SNR) = ", L_1p644/1d36

;      stop
      
   endfor
end

