pro _assemble_all_ratios_hiz
   jwst_1st_dirs, dd

   t = mrdfits(dd.ROOT+'redshift-for-sources.fits', 1, /silent)

   lines = ['OII3727', 'H_GAMMA', 'OII3727', 'NEIII_3869', 'CIII1909', $
            'H_GAMMA', 'OIII_4363', 'H_GAMMA', 'OIII_4363', 'OIII_5007', $
            'OIII_5007', 'H_BETA', 'H_ALPHA', 'H_BETA', 'H_GAMMA', 'H_BETA', $
            'NII_6584', 'H_ALPHA']

   use = where(t.confidence ge 2 and t.redshift gt 3, n_use)
   

   x = dblarr(n_use, n_elements(lines)/2)
   xlo = x
   xhi = x
   
   for itmp=0L, n_use-1 do begin
      i = use[itmp]
      tmp = strsplit(t[i].object, '_', /extract)
      srcid = tmp[1]
      print, 'Doing ', srcid
      jwst_1st_line_ratios, srcid, lines, /log, /silen, $
                            result=tmp
      x[itmp, *] = tmp.ratio
      xlo[itmp, *] = tmp.ratio_low
      xhi[itmp, *] = tmp.ratio_high
   endfor

   result = {object: t[use].object}
   for i=0L, n_elements(lines)-1, 2 do begin
      tag = lines[i]+'_'+lines[i+1]
      xtmp = dblarr(n_use)
      result = create_struct(result, tag, x[*, i/2], tag+'_LOW', xlo[*, i/2], $
                             tag+'_HIGH', xhi[*, i/2])
   endfor

   result = arr_of_struct_to_struct_of_arr(result)

   mwrfits, result, dd.specdir+'Platefit/all-ratios.fits', /create
   
;   stop

end

pro jwst_1st_line_ratios, object, lines, log=log, n_mc=n_mc, $
                          for_jupyter=for_jupyter, silent=silent, $
                          result=result

   if N_elements(n_mc) eq 0 then n_mc = 5000
   
   ;; Give ratios of adjacent lines
   if (object eq '05144') then touse = 'full' else touse = 'refull'
;   print, touse
   jwst_1st_line_luminosities, object, lines, lum=lum, dlum=dlum, /silent, touse=touse


   ;;
   ;; Use MC to get uncertainties
   ;;
   n_lines = n_elements(lines)
   n_ratios = n_lines/2
   ratios = dblarr(n_ratios)
   ratios_low = ratios
   ratios_high = ratios
   counter = 0L
   for i=0L, n_lines-1, 2 do begin
      r_1 = lum[i]+randomn(sss, n_Mc)*dlum[i]
      r_2 = lum[i+1]+randomn(sss, n_Mc)*dlum[i+1]
      if keyword_set(log) then $
       ratio = alog10(r_1/r_2) $
      else $
       ratio = r_1/r_2

      q = quantile(ratio, [0.16, 0.5, 0.84])

      ratios[counter] = q[1]
      ratios_low[counter] = q[0]
      ratios_high[counter] = q[2]

      counter = counter+1
   endfor

   if not keyword_set(silent) then begin
      
      if keyword_set(for_jupyter) then begin
         ;; Very specialised simple format
         print, format='("x=[",E12.4,"]")', ratios[0]
         print, format='("y=[",E12.4,"]")', ratios[1]
         print, format='("xlo=[",E12.4,"]")', ratios[0]-ratios_low[0]
         print, format='("ylo=[",E12.4,"]")', ratios[1]-ratios_low[1]
         print, format='("xhi=[",E12.4,"]")', ratios_high[0]-ratios[0]
         print, format='("yhi=[",E12.4,"]")', ratios_high[1]-ratios[1]
         
      endif else begin
         f = '(A,'+string(format='(I0)', n_ratios)+'(E12.4,", ")'+")"
         print, format=f, ' Ratios low=', ratios_low
         print, format=f, '    Ratios =', ratios
         print, format=f, 'Ratios high=', ratios_high
      endelse
   endif

   result = {ratio: ratios, ratio_low: ratios_low, ratio_high: ratios_high}


end

pro jwst_1st_line_luminosities, object, lines, lum=lum, dlum=dlum, silent=silent, $
                                touse=touse
   ;;
   ;; Print line luminosities in the given lines for the given object
   ;;

   if n_elements(touse) eq 0 then $
    touse = 'refull'            ; Default to renormalised spectra
   
   jwst_1st_dirs, dd
   DIR = dd.specdir+'Platefit/FluxOverviews/'

;   stop
   suffix = ''
   if size(object, /tname) eq 'STRING' then begin
      ;; Check if the last is a 'b'
      if strmid(object, 0, /reverse) eq 'b' then $
       suffix = '-second'
   endif

   bit = string(format='(I5.5)', object)
   file = DIR+'luminosities-'+bit+suffix+'.fits'

;   print, 'Reading from '+file
   all_lum = mrdfits(file, 1, /silent)

   ;; WE use renormalised fluxes
   use = where(strtrim(all_lum.labels) eq touse, n_use)
   if (n_use eq 0) then stop
   use = use[0]
   
   n_lines = n_elements(lines)
   l = dblarr(n_lines)
   dl = dblarr(n_lines)
   for i=0L, n_lines-1 do begin
      if lines[i] eq 'OII3727' then begin
         x1 = struct_var(all_lum, 'OII_3726')
         x2 = struct_var(all_lum, 'OII_3729')
         dx1 = struct_var(all_lum, 'DOII_3726')
         dx2 = struct_var(all_lum, 'DOII_3729')
         x = x1+x2
         dx = sqrt(dx1^2+dx2^2)
      endif else if lines[i] eq 'CIII1909' then begin
         x1 = struct_var(all_lum, 'CIII_1907')
         x2 = struct_var(all_lum, 'CIII_1909')
         dx1 = struct_var(all_lum, 'DCIII_1907')
         dx2 = struct_var(all_lum, 'DCIII_1909')
         x = x1+x2
         dx = sqrt(dx1^2+dx2^2)
      endif else begin
         x = struct_var(all_lum, lines[i])
         dx = struct_var(all_lum, 'D'+lines[i])
      endelse
      l[i] = x[use]
      dl[i] = dx[use]
;      stop
   endfor

   if not keyword_set(silent) then begin
      print, "Luminosities: ", l
      print, "Uncertainties: ", dl ; Scaled
   endif

   lum = l
   dlum = dl
end

   
