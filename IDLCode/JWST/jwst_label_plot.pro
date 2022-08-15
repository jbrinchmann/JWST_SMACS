pro jwst_label_plot, sp, log=log, asinh=asinh, $
                     z=z, observed=observed, rest=rest, $
                     reread=reread, lcolor=lcolor, _extra=_extra
   ;;
   ;; Simple labelled plot for a spectrum

   common linelist, lname, l_line, type

   if n_elements(lcolor) eq 0 then lcolor = 'yellow'
   
   if (n_elements(l_line) eq 0 or n_elements(type) eq 0) or keyword_set(reread) then begin
      jwst_initialize_labels
   endif

;   stop
   if max(sp.wavelength) lt 1000 then $
    micron = 1  $
   else micron = 0
   
   if keyword_set(observed) then begin
      x = sp.wavelength
   endif else begin
      if tag_exist(sp, 'restwl') then $
       x = sp.restwl $
      else if tag_exist(sp, 'z') then $
       x = sp.wavelength/(1.0+sp.z) $
      else if n_elements(z) gt 0 then $
       x = sp.wavelength/(1.0+z) $
      else begin
         print, 'Rest wavelength not found!'
         return
      endelse
   endelse

   if micron then $
    x = x*1e4

   y = sp.flux

   if keyword_set(log) then y = alog10(y) $
   else if keyword_set(asinh) then y = asinh(y)
   my_lineid_plot, x, y, l_line, lname, /extend, elinestyle=1, $
                   redshift=z, $
                   lcolor=jb_colour(lcolor), _extra=_extra
;   stop   

end

