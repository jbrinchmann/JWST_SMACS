pro jwst_1st_paschen_sfrs, pf_r
   ;; 
   ;; Calculate Paschen SFRs for the galaxies at z<3
   ;;

   jwst_1st_dirs, dd

   if n_elements(pf_r) eq 0 then $
    pf_r = mrdfits(dd.specdir+'Platefit/platefit-results-lines-fixedlinelist-renorm.fits', 1, /silent)

   
   t = mrdfits(dd.root+'redshift-for-sources.fits', 1)
   todo = ['01917', '03042', '05735', '08506', '09239', $
           '09483', '09721', '09922']

   line_to_use = ['HI_3_4', 'HI_3_6', 'HI_3_4', 'HI_3_5', $
                  'HI_3_5', 'HI_3_4', 'HI_3_5', 'HI_3_5']

;   stop
   for i=0L, n_elements(todo)-1 do begin
      iz = where(strtrim(t.object) eq '02736_'+todo[i], n)
      if n eq 0 then stop
      z = t[iz[0]].redshift

      if z lt 0 or z gt 5 then continue ; Skip these
      
      g = jwst_1st_load_one_gaussfit('02736_'+todo[i])
      
      ii = where(g.line eq line_to_use[i], n_ii)
      if (n_ii eq 0) then stop

      ;; I made a mess with units - this tries to fix this.
      med_f = median(g.flux[ii])
      if med_f lt 1d-18 and med_f gt 0 then $
       f = g.flux[ii] $
      else if med_f lt 1e-10 then $
       f = g.flux[ii]*1e-6 $
      else $
       f = g.flux[ii]*1d-17     ; Convert to ergs
      L = l_from_flux(f, z, /erg)


      upper = long(strmid(line_to_use[i], 0, /reverse))
      sfr = paschen_to_sfr(L, n_upper=upper)

      print, format='(A6,"  SFR=",F)', todo[i], sfr
;      stop
      
   endfor
   

   


end
