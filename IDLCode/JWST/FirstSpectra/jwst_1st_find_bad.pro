pro _fix_files
   ;;
   ;; First time I chose bad names This fixes it.
   ;;
   ROOT = '/data2/jarle/JWST/ERO/SMACS/'
   SPECDIR = ROOT+'Spectra/'

   badfiles = findfile(SPECDIR+'BadNames/bad-*.h5', count=n_bad)
   for i=0L, n_bad-1 do begin
      s = hdf5_to_struct(badfiles[i])

      fn = file_basename(badfiles[i])
      Len = strlen(fn)

      tmp = strmid(fn, 0, len-3) ; Without extension
      outfile1 = SPECDIR+tmp+'-f170lp-g235m.h5'
      outfile2 = SPECDIR+tmp+'-f290lp-g395m.h5'

      struct_to_hdf5, {bad: s.f170lp}, outfile1
      struct_to_hdf5, {bad: s.f290lp}, outfile2

   endfor


end

pro _find_bad_in_all

   ROOT = '/data2/jarle/JWST/ERO/SMACS/'
   SPECDIR = ROOT+'Spectra/'

   t = mrdfits(ROOT+'redshift-for-sources.fits', 1)

   for i=0L, n_elements(t)-1 do begin
      tmp = strsplit(t[i].object, '_', /extract)
      propid = string(format='(I5.5)', long(tmp[0]))
      objid = string(format='(I5.5)', long(tmp[1]))

      if objid ne 4590 then continue
      
      rootnames = 'jw'+propid+'-o00'+['7', '8']+'_s'+objid

      for j=0L, n_elements(rootnames)-1 do begin
         ;; Check whether this has been done.
         badfile = 'bad-'+rootnames[j]+'.h5'
         if file_test(badfile) then continue ; Skip
         
         d = jwst_load_nirspec_data(SPECDIR, rootnames[j], blank=0)

         ;; Do this for both filters.
         jwst_1st_find_bad, d['f170lp'].oned, bad_x_f170lp, $
                            title=rootnames[j]
         jwst_1st_find_bad, d['f290lp'].oned, bad_x_f290lp, $
                            title=rootnames[j]


         Len = strlen(badfile)
         tmp = strmid(badfile, 0, len-3) ; Without extension
         outfile1 = SPECDIR+tmp+'-f170lp-g235m.h5'
         outfile2 = SPECDIR+tmp+'-f290lp-g395m.h5'

;         stop
         struct_to_hdf5, {bad: bad_x_f170lp}, outfile1
         struct_to_hdf5, {bad: bad_x_f290lp}, outfile2

         ;; Old and bad naming 
         ;; s = {f170lp: bad_x_f170lp, f290lp: bad_x_f290lp}
         ;; struct_to_hdf5, s, badfile

      endfor
      

   endfor
   

   
end

pro jwst_1st_find_bad, sp, bad_x, title=title
   ;;
   ;; Given a spectrum, get the user to identify the bad pixels
   ;;

   ans = ''
   bad_x = []

   xrange_orig = minmax(sp.wavelength)
   plot, sp.wavelength, sp.flux, title=title

   y = sp.flux                  ; This is modified
   while (ans ne 'q') do begin
      print, 'a: add bad pixel, q: quit, h: halt, x: xrange'
      read, ans

      case ans of
         'a': begin
            print, 'Click on the line!'
            cursor, xx, yy
            dum = min(abs(sp.wavelength-xx), mi)

            bad_x = [bad_x, mi-1, mi, mi+1]
            y[mi-1:mi+1] = !values.f_nan
            plot, sp.wavelength, sp.flux, xrange=xrange, /nodata, title=title
            oplot, sp.wavelength, sp.flux,color='888888'x
            oplot, sp.wavelength, y
         end
         'h': stop
         'xr': begin
            xrange = xrange_orig
         end
         'x': begin
            print, 'Xrange to show: '
            read, x1, x2
            xrange = [x1, x2]
            plot, sp.wavelength, sp.flux, xrange=xrange, title=title
         end
         'r': begin
            plot, sp.wavelength, y, xrange=xrange, title=title
;            stop
         end
         'q': begin
            if (n_elements(bad_x) gt 0) then begin
               si = sort(bad_x)
               ui = uniq(bad_x[si])
               bad_x = bad_x[si[ui]]
               print, 'You marked ', n_elements(bad_x), ' pixels as bad.'
            endif else bad_x = -1
         end
      endcase
   endwhile

end
