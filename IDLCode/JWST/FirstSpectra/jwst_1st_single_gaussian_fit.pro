pro jwst_1st_single_gaussian_fit, sp, z
   ;;
   ;; Given a 1D spectrum, fit single gaussians to each line, fixing
   ;; the location 
   ;;

   ;; Use all the lines in the linelist
   common linelist, lname, l_line, type

   N_lines = n_elements(lname)

   restwl = sp.wavelength/(1.0+z)

   
   for i_line=0L, N_lines-1 do begin

      use = where(sp.wavelength
