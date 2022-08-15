pro jwst_1st_create_g16_mugas_pdf, object
   ;;
   ;; Create a Mugas PDF for the G16 fits because at the moment
   ;; I can't calculate it directly for technical reasons. 
   ;;


   p_Z = jwst_1st_load_pifit_pdf(object, 'LOGZ')
   p_xsi = jwst_1st_load_pifit_pdf(object, 'XI')
   p_tauv = jwst_1st_load_pifit_pdf(object, 'TAUV')
