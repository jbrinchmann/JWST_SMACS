pro jwst_vis_spec_1d2d, d, filter=filter, xrange=xrange, yrange=yrange, $
                        smooth=smooth, reextract=reextract, $
                        redshift=redshift, label=label, snimage=snimage, $
                        _extra=_extra
   ;;
   ;; Show both the 1d and 2d spectra for inspection.
   ;;

   
   
   
   ;; Get the relevant data
   if size(d, /tname) eq 'STRUCT' then $
    sp = d $
   else if n_elements(filter) eq 0 then $
    sp = d['combined'] $
   else $
    sp = d[filter]

   ;;
   ;; These should now have the same format. 
   ;;
   x = sp.oned.wavelength
   y = sp.oned.flux

   im2d = sp.twod.flux
   dim2d = sp.twod.dflux

   ;;
   ;; Layout
   ;;
   xmin = 0.1
   xmax = 0.95
   ymin = 0.1
   if keyword_set(label) then ymax = 0.85 else $
    ymax = 0.95
   ymid = 0.3
   
   pos_im = [xmin, ymin, xmax, ymid]
   pos_spec = [xmin, ymid, xmax, ymax]

   ;;
   ;; Subset
   ;;
   if n_elements(xrange) eq 0 then $
    xrange = minmax(x)


   ;; Limit to the data range
   if xrange[0] lt min(x) then xrange[0] = min(x)
   if xrange[1] gt max(x) then xrange[1] = max(x)
   
   keep = where(x ge xrange[0] and x le xrange[1], n_keep)
   if (n_keep eq 0) then stop

   x = x[keep]
   y = y[keep]
   im2d = im2d[keep, *]
   dim2d = dim2d[keep, *]
   
   ;; We now adjust the xranges to be exact as otherwise the image and
   ;; spectrum can be slightly misaligned
   xrange = minmax(x)
   
   ;; A tricky step now is that the sampling might be different. This
   ;; is fine for the line plot but the image will get distorted so we
   ;; sample onto a finer wavelength grid for the display.

   ;; This sampling can lead to slight offsets in image and spectrum
   ;; display 
;   dl = quantile(abs(x-shift(x, 1)), 0.01)
;  xnew = mkarr(min(x), max(x), dl)
   ;; So I use this instead now.
   xnew = min(x)+(max(x)-min(x))*findgen(n_elements(x))/(n_elements(x)-1)
   
   ynew = interpol(y, x, xnew)

   ;; Do the same for the image.
   dims = size(im2d, /dimen)
   imnew = fltarr(n_elements(xnew), dims[1])
   dimnew = imnew
   for i=0L, dims[1]-1 do begin
      imnew[*, i] = interpol(im2d[*, i], x,xnew)
      dimnew[*, i] = interpol(dim2d[*, i], x,xnew)
   endfor
   im2d = imnew
   dim2d = dimnew
   x = xnew
   y = ynew
   

   if n_elements(smooth) gt 0 then begin
      y = wr_smooth(x, y, smooth=smooth)
      for i=0L, dims[1]-1 do begin
         im2d[*, i] = wr_smooth(x, im2d[*, i], smooth=smooth)
         dim2d[*, i] = wr_smooth(x, dim2d[*, i], smooth=smooth)
      endfor
   endif
   

   ;; Show the image first so that subsequent plotting is on the
   ;; spectrum

   if keyword_set(snimage) then $
    im_to_display = im2d/dim2d $
   else im_to_display = im2d
   
   tvim_true, im_to_display, position=pos_im, xrange=xrange, $
              /xs, _extra=_extra
;   stop

   if keyword_set(label) then begin
      jwst_label_plot, {wavelength: x, flux: y}, $
                       observed=1, position=pos_spec, $
                       z=redshift, $
                       yrange=yrange, /xs, xtickformat='noticks', $
                       /noerase, _extra=_extra
   endif else begin
      plot, x, y, /xs, xtickformat='noticks', $
            yrange=yrange, position=pos_spec, /noerase, $
            _extra=_extra
   endelse
   

;   stop

end
