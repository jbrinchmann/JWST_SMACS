pro jwst_1st_estimate_te
   ;;
   ;; Estimate the electron temperature from double ratios
   ;;

   todo = ['04590', '05144', '06355', '08140', '10612']


   for i=0L, n_elements(todo)-1 do begin
      
      jwst_1st_line_luminosities, todo[i], ['OIII_4363', 'H_GAMMA', 'OIII_5007', 'H_BETA'], $
                                  lum=lum, dlum=dlum, /silent

      ;; For comparison I would also like to do this on the
      ;; non-corrected fluxes.
      jwst_1st_line_luminosities, todo[i], ['OIII_4363', 'H_GAMMA', 'OIII_5007', 'H_BETA'], $
                                  lum=lum_nc, dlum=dlum_nc, /silent, touse='full'
      

      res_double = te_pyneb_o3_mc(lum, dlum)
      res_standard = te_pyneb_o3_mc(lum[[0,2]], dlum[[0,2]])

      res_double_nc = te_pyneb_o3_mc(lum, dlum)
      res_standard_nc = te_pyneb_o3_mc(lum[[0,2]], dlum[[0,2]])
      
      if i eq 0 then begin
         results_d = replicate(res_double, n_elements(todo))
         results_s = replicate(res_standard, n_elements(todo))
         results_d_nc = replicate(res_double, n_elements(todo))
         results_s_nc = replicate(res_standard, n_elements(todo))
      endif

      results_d[i] = res_double
      results_s[i] = res_standard
      results_d_nc[i] = res_double_nc
      results_s_nc[i] = res_standard_nc
      
      ;; Print the results.
      f = '(A5,2X,4(F7.1,2X))'
      print, format=f, todo[i], res_double.p50, res_standard.p50, res_double_nc.p50, $
             res_standard_nc.p50
      
   endfor

   ;;
   ;; The effect of the renormalisation is modest
   ;;
   jwst_1st_dirs, dd

   mwrfits, results_d, dd.root+'TeModeling/Te_double_ratio.fits', /create
   mwrfits, results_s, dd.root+'TeModeling/Te_standard_ratio.fits', /create

   


   
end
  
