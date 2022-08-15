pro jwst_1st_get_redshifts, i, data=data, _extra=_extra

   DIR = '/data2/jarle/JWST/ERO/SMACS/Spectra/'

   
   t_relics = mrdfits(DIR+'../matched-to-relics.fits', 1)
   
   source_list = DIR+'source_list_g235m.txt' ; Same for f395m
   readcol, source_list, fnames, srcnames, format='A,A', skip=1, /silent

   todo = srcnames[i]
   if strmid(todo, 0, 4) eq 'back' then begin
      print, 'Warning: you are trying to look at a background spectrum -continue?'
      ans = ''
      read, ans
      if strlowcase(ans) ne 'y' then return
   endif

   
   ;;
   ;; Create name so we can load both exposures.
   ;;
   tmp = strsplit(todo, '_', /extract)
   propid = string(format='(I5.5)', long(tmp[0]))
   objid = string(format='(I5.5)', long(tmp[1]))
   
   rootnames = 'jw'+propid+'-o00'+['7', '8']+'_s'+objid
;   stop
   ;; Get the data

   if (n_elements(data) eq 0) then begin
      d1 = jwst_load_nirspec_data(DIR, rootnames[0], /blank)
      d2 = jwst_load_nirspec_data(DIR, rootnames[1], /blank)
      d = jwst_combine_nirspec_exposures(d1, d2, factor=10.)
   endif else begin
      d = data
   endelse
   print, 'Photometric redshift estimates: '
   print, format='("  ZB(ACS)=",F0.2,"  ([",F0.2,"-",F0.2,"])")', $
          t_relics[i].zb_acs, t_relics[i].zbmin_acs, t_relics[i].zbmax_acs
   print, format='("  ZB(WF3)=",F0.2,"  ([",F0.2,"-",F0.2,"])")', $
          t_relics[i].zb_wf3, t_relics[i].zbmin_wf3, t_relics[i].zbmax_wf3
   
   jwst_view_spec_findz, d['combined'], result=res, _extra=_extra


   openw, lun, DIR+'../redshift-log.dat', /get_lun, /append
   f = '(A10,2X,I3,2X,F9.5)'
   printf, lun, format=f, todo, i, res.z
   free_lun, lun

end
