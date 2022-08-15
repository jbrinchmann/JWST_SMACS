pro jwst_1st_convert_match_catalogue

   ROOT = '/data2/jarle/JWST/ERO/SMACS/'
   readcol, ROOT+'match-ids.dat', specid, relicsid, relicsacs, $
            nircamID, RaNirCam, DecNIRCam, format='A,L,L,L,D,D'
   tmp = {specid: '', relics_wf3_id: 0L, relics_acs_id: 0L, $
          nircam_f200w_id: 0L, nircam_ra: 0.0d0, nircam_dec: 0.0d0}

   cat = replicate(tmp, n_elements(specid))
   cat.specid = specid
   cat.relics_wf3_id = relicsid
   cat.relics_acs_id = relicsacs
   cat.nircam_f200w_id = nircamID
   cat.nircam_ra = RaNIRCam
   cat.nircam_dec = DecNIRCam

   mwrfits, cat, ROOT+'match-ids.fits', /create
   
end
