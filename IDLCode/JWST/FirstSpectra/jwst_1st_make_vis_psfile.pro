pro jwst_1st_make_vis_psfile, ps=ps
   ROOT = '/data2/jarle/JWST/ERO/SMACS/'
   t = mrdfits(ROOT+'redshift-for-sources.fits', 1)

   if keyword_set(ps) then begin
      ps_on, ROOT+'Figures/vis-all.ps', /times, /color, aspect=0.4, xlen=12
      !p.font = 0
      !x.thick = 2
      !y.thick = 2
   endif
   for i=0L, n_elements(t)-2 do begin ; The last is duplicated
      tmp = strsplit(t[i].object, '_', /extract)
      propid = long(tmp[0])
      objid = long(tmp[1])

      d = jwst_1st_get_coadded(objid)

      if t[i].confidence ge 2 then begin
         jwst_vis_spec_1d2d, d,  charsize=1.1, pcharsize=1.1, /reread, $
                             lcharsize=1.1, /zscale, redshift=t[i].redshift, $
                             /label, range=range_this, lcolor='orange'
      endif else begin
         jwst_vis_spec_1d2d, d,  charsize=1.1, pcharsize=1.1, /reread, $
                             lcharsize=1.1, /zscale, range=range_this, $
                             lcolor='orange'
      endelse

      jb_text, -0.05, -0.08, t[i].object, charsize=1.5, /frac, $
               /relative, color=jb_colour("CornFlowerBlue"), $
               align=1
      

      if not keyword_set(ps) then stop

   endfor

   if keyword_set(ps) then begin
      ps_off
      my_cleanplot, /silent
   endif
end
