function jwst_1st_get_each_obs, number, filter=filter, blank=blank, $
                                keep_bad=keep_bad

   if n_elements(blank) eq 0 then blank = 1

   
   DIR = '/data2/jarle/JWST/ERO/SMACS/Spectra/'
   
   propid = string(format='(I5.5)', 2736)
   objid = string(format='(I5.5)', number)
   
   rootnames = 'jw'+propid+'-o00'+['7', '8']+'_s'+objid
   ;; Get the data

   d1 = jwst_load_nirspec_data(DIR, rootnames[0], blank=blank, $
                              keep_bad=keep_bad)
   d2 = jwst_load_nirspec_data(DIR, rootnames[1], blank=blank, $
                               keep_bad=keep_bad)


   return, {obs07: d1, obs08: d2}
end
