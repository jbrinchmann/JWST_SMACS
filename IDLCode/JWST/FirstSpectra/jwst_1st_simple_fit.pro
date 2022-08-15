pro jwst_1st_simple_fit, sp, result


   ii4363 = where(abs(sp.restwl-4363) lt 15)
   iihg = where(abs(sp.restwl-4340) lt 15)
   ii5007 = where(abs(sp.restwl-5007) lt 15)


   g4363 = fit_one_gaussian_mp(sp.restwl[ii4363], sp.flux[ii4363], a4363, $
                               err=sp.dflux[ii4363], $
                               perror=perror4363, /return_flux)

   ghg = fit_one_gaussian_mp(sp.restwl[iihg], sp.flux[iihg], ahg, $
                               err=sp.dflux[iihg], $
                               perror=perrorhg, /return_flux)

   g5007 = fit_one_gaussian_mp(sp.restwl[ii5007], sp.flux[ii5007], a5007, $
                               err=sp.dflux[ii5007], $
                               perror=perror5007, /return_flux)


   N = 5000
   r4363 = a4363[0]+randomn(SS, N, /normal)*perror4363[0]
   r5007 = a5007[0]+randomn(SS, N, /normal)*perror5007[0]
   rhg = ahg[0]+randomn(ss, N)*perrorhg[0]

   ratio = r5007/r4363
   v = quantile(ratio, [0.16, 0.5, 0.84])

;   stop
   
   ratio = r4363/rhg
   w = quantile(ratio, [0.16, 0.5, 0.84])

   result = {r_5007_4363: v[1], r_5007_4363_low: v[0], r_5007_4363_high: v[2], $
             r_4363_hg: w[1], r_4363_hg_low: w[0], r_4363_hg_high: w[2], $
             f_hg: ahg[0], df_hg: perrorhg[0], f_4363: a4363[0], $
             df_4363: perror4363[0], f5007: a5007[0], df5007: perror5007[0]}

   
end
