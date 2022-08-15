pro jwst_1st_match_to_relics_cat
   ;;
   ;; Match the NIRSpec data to the RELICS catalogues.
   ;;
   
   DIR = '/data2/jarle/JWST/ERO/SMACS/'

   ;; Get the two RELICS catalogues
   t_acs = mrdfits(DIR+'hlsp_relics_hst_acs-wfc3ir_smacs0723-73_multi_v1_cat-noheader.fits', 1)
   t_wfc3 = mrdfits(DIR+'hlsp_relics_hst_wfc3ir_smacs0723-73_multi_v1_cat-noheader.fits', 1)

   ;; Then the match id catalogue
   m = mrdfits(DIR+'match-ids.fits', 1)
   N_nirspec = n_elements(m)
   
   ;; We now want to make one monster catalogue putting _ACS and _WF3
   ;; as suffices to each tag
   tags_acs = tag_names(t_acs)
   tags_wf3 = tag_names(t_wfc3)

   tmp = {specid: m.specid, nircam_ra: m.nircam_ra, $
          nircam_dec: m.nircam_dec, i_acs: m.relics_acs_id, $
          i_wf3: m.relics_wf3_id}

   ;; Lookups - I want an index array for looking into the ACS table
   is_acs_match = where(m.relics_acs_id ge 0, n_acs)
   print, format='("For ",I0," I have a match to the ACS catalogue")', n_acs
   i_acs = lonarr(n_acs)
   for i=0L, n_acs-1 do $
    i_acs[i] = where(m[is_acs_match[i]].relics_acs_id eq t_acs.id)

   is_wf3_match = where(m.relics_wf3_id ge 0, n_wf3)
   print, format='("For ",I0," I have a match to the WFC3 catalogue")', n_wf3
   i_wf3 = lonarr(n_wf3)
   for i=0L, n_wf3-1 do $
    i_wf3[i] = where(m[is_wf3_match[i]].relics_wf3_id eq t_wfc3.id)


   
   for i=0L, n_elements(tags_acs)-1 do begin
      xall = struct_var(t_acs, tags_acs[i])
      x = make_array(N_nirspec, type=size(xall, /type), value=-999.)
      x[is_acs_match] = xall[i_acs]
      tmp = create_struct(tmp, tags_acs[i]+'_ACS', x)
   endfor

   
   for i=0L, n_elements(tags_wf3)-1 do begin
      xall = struct_var(t_wfc3, tags_wf3[i])
      x = make_array(N_nirspec, type=size(xall, /type), value=-999.)
      x[is_wf3_match] = xall[i_wf3]
      tmp = create_struct(tmp, tags_wf3[i]+'_WF3', x)
   endfor


   tmp = struct_of_arr_to_arr_of_struct(tmp)
   
   mwrfits, tmp, DIR+'matched-to-relics.fits', /create

   
end
