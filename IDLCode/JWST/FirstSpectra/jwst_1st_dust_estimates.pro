pro _print_tv, s, what
   print, format='("TauV(",A,")=",F6.2,"^{+",F6.3,"}_{-",F6.3,"}")', $
          what, s.p50, s.p50-s.p16, s.p84-s.p50
end

pro jwst_1st_dust_estimates
   ;;
   ;; Estimate dust for galaxies. This is the simple approach. 
   ;;

   jwst_1st_dirs, dd
   t_z = mrdfits(dd.root+'redshift-for-sources.fits', 1, /silent)


   ;;
   ;; 1917 - this has HI 3-6 (weakly), 3-5, 3-4
   ;; 
   flx = jwst_1st_load_flux_overview(1917)
   t = jwst_1st_load_one_gaussfit('02736_01917')
   i_hi_34 = where(t.line eq 'HI_3_4')
   i_hi_35 = where(t.line eq 'HI_3_5')
   i_hi_36 = where(t.line eq 'HI_3_6')

   i_renorm = where(strtrim(flx.labels) eq 'refull')
   
   ;; This has poor-ish platefit fits
   r_34_35_pf = [flx.hi_3_4[i_renorm], flx.hi_3_5[i_renorm]]
   r_34_35_g = [t.flux[i_hi_34], t.flux[i_hi_35]]
   r_34_36_g = [t.flux[i_hi_34], t.flux[i_hi_36]]
   
   dr_34_35_pf = [flx.dhi_3_4[i_renorm], flx.dhi_3_5[i_renorm]]
   dr_34_35_g = [t.dflux[i_hi_34], t.dflux[i_hi_35]]
   dr_34_36_g = [t.dflux[i_hi_34], t.dflux[i_hi_36]]
   
   tV_34_35_pf = tau_from_hydrogen(5500., 3, 4, 3, 5, r_34_35_pf, $
                                   dfluxes=dr_34_35_pf, /steep)
   tV_34_35_g = tau_from_hydrogen(5500., 3, 4, 3, 5, r_34_35_g, $
                                  dfluxes=dr_34_35_g, /steep)

   tV_34_36_g = tau_from_hydrogen(5500., 3, 4, 3, 6, r_34_36_g, $
                                  dfluxes=dr_34_36_g, /steep)

   print, 'SpecID: 1917'
   print, '----------------'
   _print_tv, tV_34_35_pf, '3-4,3-5,Pfit'
   _print_tv, tV_34_35_g, '3-4,3-5,Gaussian'
   _print_tv, tV_34_36_g, '3-4,3-6,Gaussian'


   ;;----------------------------------
   ;;   3042
   ;;   This has Ha and Paschen lines
   ;;----------------------------------
   flx = jwst_1st_load_flux_overview(3042)
   t = jwst_1st_load_one_gaussfit('02736_03042')
   i_hi_23 = where(t.line eq 'H_ALPHA')
   i_hi_35 = where(t.line eq 'HI_3_5')
   i_hi_37 = where(t.line eq 'HI_3_7') ; dubious
   
   
   ;; This has poor-ish platefit fits
   r_35_23_pf = [flx.hi_3_5[i_renorm], flx.h_alpha[i_renorm]]
   r_35_37_pf = [flx.hi_3_5[i_renorm], flx.hi_3_7[i_renorm]]
   r_35_23_g = [t.flux[i_hi_35], t.flux[i_hi_23]]
   r_35_37_g = [t.flux[i_hi_35], t.flux[i_hi_37]]

   dr_35_23_pf = [flx.dhi_3_5[i_renorm], flx.dh_alpha[i_renorm]]
   dr_35_37_pf = [flx.dhi_3_5[i_renorm], flx.dhi_3_7[i_renorm]]
   dr_35_23_g = [t.dflux[i_hi_35], t.dflux[i_hi_23]]
   dr_35_37_g = [t.dflux[i_hi_35], t.dflux[i_hi_37]]

   
   tV_35_23_pf = tau_from_hydrogen(5500., 3, 5, 2, 3, r_35_23_pf, $
                                   dfluxes=dr_35_23_pf, /steep)
   tV_35_23_g = tau_from_hydrogen(5500., 3, 5, 2, 3, r_35_23_g, $
                                  dfluxes=dr_35_23_g, /steep)

   tV_35_37_pf = tau_from_hydrogen(5500., 3, 5, 3, 7, r_35_37_pf, $
                                   dfluxes=dr_35_37_pf, /steep)
   tV_35_37_g = tau_from_hydrogen(5500., 3, 5, 3, 7, r_35_37_g, $
                                  dfluxes=dr_35_37_g, /steep)
   
   print, 'SpecID: 3042'
   print, '----------------'
   _print_tv, tV_35_23_pf, '3-5,2-3,Pfit'
   _print_tv, tV_35_23_g, '3-5,2-3,Gaussian'
   _print_tv, tV_35_37_pf, '3-5,3-7,Pfit'
   _print_tv, tV_35_37_g, '3-5,3-7,Gaussian'

   
   ;;----------------------------
   ;;    4590
   ;;    This has Hd, Hg and Hb
   ;;----------------------------
   flx = jwst_1st_load_flux_overview(4590)
   t = jwst_1st_load_one_gaussfit('02736_04590')
   i_hi_24 = where(t.line eq 'H_BETA')
   i_hi_25 = where(t.line eq 'H_GAMMA')
   i_hi_26 = where(t.line eq 'H_DELTA')


   ;; This has poor-ish platefit fits
   r_hb_hg_pf = [flx.h_beta[i_renorm], flx.h_gamma[i_renorm]]
   r_hb_hd_pf = [flx.h_beta[i_renorm], flx.h_delta[i_renorm]]
   r_hb_hg_g = [t.flux[i_hi_24], t.flux[i_hi_25]]
   r_hb_hd_g = [t.flux[i_hi_24], t.flux[i_hi_26]]

   dr_hb_hg_pf = [flx.dh_beta[i_renorm], flx.dh_gamma[i_renorm]]
   dr_hb_hd_pf = [flx.dh_beta[i_renorm], flx.dh_delta[i_renorm]]
   dr_hb_hg_g = [t.dflux[i_hi_24], t.dflux[i_hi_25]]
   dr_hb_hd_g = [t.dflux[i_hi_24], t.dflux[i_hi_26]]

   
   tV_hb_hg_pf = tau_from_hydrogen(5500., 2, 4, 2, 5, r_hb_hg_pf, $
                                   dfluxes=dr_hb_hg_pf, /steep)
   tV_hb_hg_g = tau_from_hydrogen(5500., 2, 4, 2, 5, r_hb_hg_g, $
                                  dfluxes=dr_hb_hg_g, /steep)

   tV_hb_hd_pf = tau_from_hydrogen(5500., 2, 4, 2, 6, r_hb_hd_pf, $
                                   dfluxes=dr_hb_hd_pf, /steep)
   tV_hb_hd_g = tau_from_hydrogen(5500., 2, 4, 2, 6, r_hb_hd_g, $
                                  dfluxes=dr_hb_hd_g, /steep)


   print, 'SpecID: 4590'
   print, '----------------'
   _print_tv, tV_hb_hg_pf, 'Hb/Hg,Pfit'
   _print_tv, tV_hb_hg_g, 'Hb/Hg,Gaussian'
   _print_tv, tV_hb_hd_pf, 'Hb/Hd,Pfit'
   _print_tv, tV_hb_hd_g, 'Hb/Hd,Gaussian'


      
   ;;----------------------------
   ;;    5144
   ;;    This has Hd, Hg,  Hb and Ha 
   ;;----------------------------
   flx = jwst_1st_load_flux_overview(5144)
   t = jwst_1st_load_one_gaussfit('02736_05144')
   i_hi_23 = where(t.line eq 'H_ALPHA')
   i_hi_24 = where(t.line eq 'H_BETA')
   i_hi_25 = where(t.line eq 'H_GAMMA')
   i_hi_26 = where(t.line eq 'H_DELTA')


   ;; This has poor-ish platefit fits
   r_hb_hg_pf = [flx.h_beta[i_renorm], flx.h_gamma[i_renorm]]
   r_hb_hd_pf = [flx.h_beta[i_renorm], flx.h_delta[i_renorm]]
   r_ha_hb_pf = [flx.h_alpha[i_renorm], flx.h_beta[i_renorm]]
   r_hb_hg_g = [t.flux[i_hi_24], t.flux[i_hi_25]]
   r_hb_hd_g = [t.flux[i_hi_24], t.flux[i_hi_26]]
   r_ha_hb_g = [t.flux[i_hi_23], t.flux[i_hi_24]]

   dr_hb_hg_pf = [flx.dh_beta[i_renorm], flx.dh_gamma[i_renorm]]
   dr_hb_hd_pf = [flx.dh_beta[i_renorm], flx.dh_delta[i_renorm]]
   dr_ha_hb_pf = [flx.dh_alpha[i_renorm], flx.dh_beta[i_renorm]]
   dr_hb_hg_g = [t.dflux[i_hi_24], t.dflux[i_hi_25]]
   dr_hb_hd_g = [t.dflux[i_hi_24], t.dflux[i_hi_26]]
   dr_ha_hb_g = [t.dflux[i_hi_23], t.dflux[i_hi_24]]

   
   tV_hb_hg_pf = tau_from_hydrogen(5500., 2, 4, 2, 5, r_hb_hg_pf, $
                                   dfluxes=dr_hb_hg_pf, /steep)
   tV_hb_hg_g = tau_from_hydrogen(5500., 2, 4, 2, 5, r_hb_hg_g, $
                                  dfluxes=dr_hb_hg_g, /steep)

   tV_hb_hd_pf = tau_from_hydrogen(5500., 2, 4, 2, 6, r_hb_hd_pf, $
                                   dfluxes=dr_hb_hd_pf, /steep)
   tV_hb_hd_g = tau_from_hydrogen(5500., 2, 4, 2, 6, r_hb_hd_g, $
                                  dfluxes=dr_hb_hd_g, /steep)

   tV_ha_hb_pf = tau_from_hydrogen(5500., 2, 3, 2, 4, r_ha_hb_pf, $
                                   dfluxes=dr_ha_hb_pf, /steep)
   tV_ha_hb_g = tau_from_hydrogen(5500., 2, 3, 2, 4, r_ha_hb_g, $
                                  dfluxes=dr_ha_hb_g, /steep)


   print, 'SpecID: 5144'
   print, '----------------'
   _print_tv, tV_hb_hg_pf, 'Hb/Hg,Pfit'
   _print_tv, tV_hb_hg_g, 'Hb/Hg,Gaussian'
   _print_tv, tV_hb_hd_pf, 'Hb/Hd,Pfit'
   _print_tv, tV_hb_hd_g, 'Hb/Hd,Gaussian'
   _print_tv, tV_ha_hb_pf, 'Ha/Hb,Pfit'
   _print_tv, tV_ha_hb_g, 'Ha/Hb,Gaussian'


   ;;---------------------
   ;;  5735
   ;;
   ;; This has Pb & Pa.
   ;;---------------------
   flx = jwst_1st_load_flux_overview(5735)
   t = jwst_1st_load_one_gaussfit('02736_05735')
   i_hi_35 = where(t.line eq 'HI_3_5')
   i_hi_34 = where(t.line eq 'HI_3_4')
   
   
   ;; I'll just use the Gaussfits for these objects
   r_34_35_g = [t.flux[i_hi_34], t.flux[i_hi_35]]
   dr_34_35_g = [t.dflux[i_hi_35], t.dflux[i_hi_37]]
   
   tV_34_35_g = tau_from_hydrogen(5500., 3, 4, 3, 5, r_34_35_g, $
                                  dfluxes=dr_34_35_g, /steep)
   
   print, 'SpecID: 5735'
   print, '----------------'
   _print_tv, tV_34_35_g, '3-4,3-5,Gaussian'

   
   ;;----------------------------
   ;;    6355
   ;;    This has Hd, Hg,  Hb
   ;;----------------------------
   flx = jwst_1st_load_flux_overview(6355)
   t = jwst_1st_load_one_gaussfit('02736_06355')
   i_hi_24 = where(t.line eq 'H_BETA')
   i_hi_25 = where(t.line eq 'H_GAMMA')
   i_hi_26 = where(t.line eq 'H_DELTA')


   ;; We use platefit as reference
   r_hb_hg_pf = [flx.h_beta[i_renorm], flx.h_gamma[i_renorm]]
   r_hb_hd_pf = [flx.h_beta[i_renorm], flx.h_delta[i_renorm]]
   r_hb_hg_g = [t.flux[i_hi_24], t.flux[i_hi_25]]
   r_hb_hd_g = [t.flux[i_hi_24], t.flux[i_hi_26]]

   dr_hb_hg_pf = [flx.dh_beta[i_renorm], flx.dh_gamma[i_renorm]]
   dr_hb_hd_pf = [flx.dh_beta[i_renorm], flx.dh_delta[i_renorm]]
   dr_hb_hg_g = [t.dflux[i_hi_24], t.dflux[i_hi_25]]
   dr_hb_hd_g = [t.dflux[i_hi_24], t.dflux[i_hi_26]]

   
   tV_hb_hg_pf = tau_from_hydrogen(5500., 2, 4, 2, 5, r_hb_hg_pf, $
                                   dfluxes=dr_hb_hg_pf, /steep)
   tV_hb_hg_g = tau_from_hydrogen(5500., 2, 4, 2, 5, r_hb_hg_g, $
                                  dfluxes=dr_hb_hg_g, /steep)

   tV_hb_hd_pf = tau_from_hydrogen(5500., 2, 4, 2, 6, r_hb_hd_pf, $
                                   dfluxes=dr_hb_hd_pf, /steep)
   tV_hb_hd_g = tau_from_hydrogen(5500., 2, 4, 2, 6, r_hb_hd_g, $
                                  dfluxes=dr_hb_hd_g, /steep)


   print, 'SpecID: 6355'
   print, '----------------'
   _print_tv, tV_hb_hg_pf, 'Hb/Hg,Pfit'
   _print_tv, tV_hb_hg_g, 'Hb/Hg,Gaussian'
   _print_tv, tV_hb_hd_pf, 'Hb/Hd,Pfit'
   _print_tv, tV_hb_hd_g, 'Hb/Hd,Gaussian'



   ;;----------------------------
   ;;    8140
   ;;    This has Hd, Hg,  Hb, Ha - indeed also H_9, H_10 and H_11 in
   ;;    the Gaussian fit
   ;;----------------------------
   flx = jwst_1st_load_flux_overview(8140)
   t = jwst_1st_load_one_gaussfit('02736_08140')
   i_hi_23 = where(t.line eq 'H_ALPHA')
   i_hi_24 = where(t.line eq 'H_BETA')
   i_hi_25 = where(t.line eq 'H_GAMMA')
   i_hi_26 = where(t.line eq 'H_DELTA')
   i_hi_29 = where(t.line eq 'H9')


   r_hb_hg_pf = [flx.h_beta[i_renorm], flx.h_gamma[i_renorm]]
   r_hb_hd_pf = [flx.h_beta[i_renorm], flx.h_delta[i_renorm]]
   r_ha_hb_pf = [flx.h_alpha[i_renorm], flx.h_beta[i_renorm]]
   r_hb_hg_g = [t.flux[i_hi_24], t.flux[i_hi_25]]
   r_hb_hd_g = [t.flux[i_hi_24], t.flux[i_hi_26]]
   r_ha_hb_g = [t.flux[i_hi_23], t.flux[i_hi_24]]
   r_hb_h9_g = [t.flux[i_hi_24], t.flux[i_hi_29]]

   dr_hb_hg_pf = [flx.dh_beta[i_renorm], flx.dh_gamma[i_renorm]]
   dr_hb_hd_pf = [flx.dh_beta[i_renorm], flx.dh_delta[i_renorm]]
   dr_ha_hb_pf = [flx.dh_alpha[i_renorm], flx.dh_beta[i_renorm]]
   dr_hb_hg_g = [t.dflux[i_hi_24], t.dflux[i_hi_25]]
   dr_hb_hd_g = [t.dflux[i_hi_24], t.dflux[i_hi_26]]
   dr_ha_hb_g = [t.dflux[i_hi_23], t.dflux[i_hi_24]]
   dr_hb_h9_g = [t.dflux[i_hi_24], t.dflux[i_hi_29]]

   
   tV_hb_hg_pf = tau_from_hydrogen(5500., 2, 4, 2, 5, r_hb_hg_pf, $
                                   dfluxes=dr_hb_hg_pf, /steep)
   tV_hb_hg_g = tau_from_hydrogen(5500., 2, 4, 2, 5, r_hb_hg_g, $
                                  dfluxes=dr_hb_hg_g, /steep)

   tV_hb_hd_pf = tau_from_hydrogen(5500., 2, 4, 2, 6, r_hb_hd_pf, $
                                   dfluxes=dr_hb_hd_pf, /steep)
   tV_hb_hd_g = tau_from_hydrogen(5500., 2, 4, 2, 6, r_hb_hd_g, $
                                  dfluxes=dr_hb_hd_g, /steep)

   tV_ha_hb_pf = tau_from_hydrogen(5500., 2, 3, 2, 4, r_ha_hb_pf, $
                                   dfluxes=dr_ha_hb_pf, /steep)
   tV_ha_hb_g = tau_from_hydrogen(5500., 2, 3, 2, 4, r_ha_hb_g, $
                                  dfluxes=dr_ha_hb_g, /steep)

   tV_hb_h9_g = tau_from_hydrogen(5500., 2, 4, 2, 9, r_hb_h9_g, $
                                  dfluxes=dr_hb_h9_g, /steep)


   print, 'SpecID: 8140'
   print, '----------------'
   _print_tv, tV_hb_hg_pf, 'Hb/Hg,Pfit'
   _print_tv, tV_hb_hg_g, 'Hb/Hg,Gaussian'
   _print_tv, tV_hb_hd_pf, 'Hb/Hd,Pfit'
   _print_tv, tV_hb_hd_g, 'Hb/Hd,Gaussian'
   _print_tv, tV_ha_hb_pf, 'Ha/Hb,Pfit'
   _print_tv, tV_ha_hb_g, 'Ha/Hb,Gaussian'
   _print_tv, tV_hb_h9_g, 'Hb/H9,Gaussian'


   ;;----------------------------------
   ;;   8506
   ;;   This has Ha and Paschen lines
   ;;----------------------------------
   flx = jwst_1st_load_flux_overview(8506)
   t = jwst_1st_load_one_gaussfit('02736_08506')
   i_hi_23 = where(t.line eq 'H_ALPHA')
   i_hi_35 = where(t.line eq 'HI_3_5')
   i_hi_36 = where(t.line eq 'HI_3_6')
   i_hi_37 = where(t.line eq 'HI_3_7')
   
   
   r_35_23_pf = [flx.hi_3_5[i_renorm], flx.h_alpha[i_renorm]]
   r_35_37_pf = [flx.hi_3_5[i_renorm], flx.hi_3_7[i_renorm]]
   r_35_23_g = [t.flux[i_hi_35], t.flux[i_hi_23]]
   r_35_36_g = [t.flux[i_hi_35], t.flux[i_hi_36]]
   r_35_37_g = [t.flux[i_hi_35], t.flux[i_hi_37]]

   dr_35_23_pf = [flx.dhi_3_5[i_renorm], flx.dh_alpha[i_renorm]]
   dr_35_37_pf = [flx.dhi_3_5[i_renorm], flx.dhi_3_7[i_renorm]]
   dr_35_23_g = [t.dflux[i_hi_35], t.dflux[i_hi_23]]
   dr_35_36_g = [t.dflux[i_hi_35], t.dflux[i_hi_36]]
   dr_35_37_g = [t.dflux[i_hi_35], t.dflux[i_hi_37]]
   
   
   tV_35_23_pf = tau_from_hydrogen(5500., 3, 5, 2, 3, r_35_23_pf, $
                                   dfluxes=dr_35_23_pf, /steep)
   tV_35_23_g = tau_from_hydrogen(5500., 3, 5, 2, 3, r_35_23_g, $
                                  dfluxes=dr_35_23_g, /steep)

   tV_35_36_g = tau_from_hydrogen(5500., 3, 5, 3, 6, r_35_36_g, $
                                  dfluxes=dr_35_36_g, /steep)
   tV_35_37_pf = tau_from_hydrogen(5500., 3, 5, 3, 7, r_35_37_pf, $
                                   dfluxes=dr_35_37_pf, /steep)
   tV_35_37_g = tau_from_hydrogen(5500., 3, 5, 3, 7, r_35_37_g, $
                                  dfluxes=dr_35_37_g, /steep)
   
   print, 'SpecID: 8506'
   print, '----------------'
   _print_tv, tV_35_23_pf, '3-5,2-3,Pfit'
   _print_tv, tV_35_23_g, '3-5,2-3,Gaussian'
   _print_tv, tV_35_36_g, '3-5,3-6,Gaussian'
   _print_tv, tV_35_37_pf, '3-5,3-7,Pfit'
   _print_tv, tV_35_37_g, '3-5,3-7,Gaussian'


   
   stop
   
end
   
   

