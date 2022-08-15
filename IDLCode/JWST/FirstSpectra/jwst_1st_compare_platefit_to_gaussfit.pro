pro jwst_1st_compare_platefit_to_gaussfit, pf_r, pf_o
   ;;
   ;; See how they compare
   ;;
   jwst_1st_dirs, dd

   
   ;; Get the platefit results both on the renormalised and normalised
   ;; versions.
   if n_elements(pf_r) eq 0 then $
    pf_r = mrdfits(dd.specdir+'Platefit/platefit-results-lines-fixedlinelist-renorm.fits', 1, /silent)
   if n_elements(pf_o) eq 0 then $
    pf_o = mrdfits(dd.specdir+'Platefit/platefit-results-lines-fixedlinelist.fits', 1, /silent)

   outdir = dd.specdir+"Redshifts/"
   t_z = mrdfits(dd.root+'redshift-for-sources.fits', 1, /silent)

   ;; I will ignore the provenance.
   flux_gaussian = []
   dflux_gaussian = []
   flux_pf_r = []
   dflux_pf_r = []
   flux_pf_o = []
   dflux_pf_o = []
   z = []
   origin = []
   line = []

   z_pf = pf_r.z
   
   for i=0L, n_elements(t_z)-1 do begin
      if (t_z[i].confidence lt 2) then continue

      ;; Match this object to the PF file.
      tmp = strsplit(t_z[i].object, '_', /extract)
      pf_r_key = string(format='(I5.5," (renorm)")', long(tmp[1]))
      pf_o_key = string(format='(I5.5)', long(tmp[1]))

      i_obj_r = where(strtrim(pf_r.specid) eq pf_r_key, n_match)
      if (n_match eq 0) then continue ; Skip these
      i_obj_o = where(strtrim(pf_o.specid) eq pf_o_key, n_match)
      if (n_match eq 0) then continue ; Skip these
      
      t_g =  jwst_1st_load_one_gaussfit(t_z[i].object)

      for i_line=0L, n_elements(t_g.line)-1 do begin
         key = t_g.line[i_line]+'_FLUX'
         xr = struct_var(pf_r, key, /silent, status=status)
         if status eq 0 then continue ; Skip this line
         dxr = struct_var(pf_r, key+'_ERR', /silent)
         xo = struct_var(pf_o, key, /silent, status=status)
         if status eq 0 then continue ; Skip this line
         dxo = struct_var(pf_o, key+'_ERR', /silent)

         print, 'Using '+key+' for '+t_z[i].object
         flux_gaussian = [flux_gaussian, t_g.flux[i_line]]
         dflux_gaussian = [dflux_gaussian, t_g.dflux[i_line]]

         flux_pf_r = [flux_pf_r, xr[i_obj_r]]
         dflux_pf_r = [dflux_pf_r, dxr[i_obj_r]]

         flux_pf_o = [flux_pf_o, xr[i_obj_o]]
         dflux_pf_o = [dflux_pf_o, dxr[i_obj_o]]

         z = [z, t_z[i].redshift]
         origin = [origin, t_z[i].object]
         line = [line, t_g.line[i_line]]
      endfor
      

   endfor
   
   stop
   

end
