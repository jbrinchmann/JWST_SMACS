pro jwst_1st_compare_hdhg_to_curti, label=label, $
                                    ps=ps
   ;;
   ;; Compare Hd/Hb and Hg/Hb to the Curti et al results. 
   ;; and contrast this to the non-renormalised spectrum

   if (n_elements(label) eq 0) then $
    label = 'refull'
   
   objects = [4590, 6355, 10612]
   hd_curti = [0.16, 0.26, 0.17]
   dhd_curti = [0.02, 0.02, 0.03]
   hg_curti = [0.41, 0.46, 0.44]
   dhg_curti=[0.03, 0.02, 0.03]

   ;; Renorm results

   jwst_1st_dirs, dd
   DIR = dd.specdir+'Platefit/FluxOverviews/'
   
   
   hdhb_my = fltarr(n_elements(objects))
   hghb_my = hdhb_my
   dhdhb_my = hdhb_my
   dhghb_my = hdhb_my

   hdhb_my_ref = hdhb_my
   hghb_my_ref = hdhb_my
   dhdhb_my_ref = hdhb_my
   dhghb_my_ref = hdhb_my

   for i=0L, n_elements(objects)-1 do begin
      fname = string(format='("fluxes-",I5.5,"-norm.fits")', objects[i])
      t = mrdfits(DIR+fname, 1)
      use = where(strtrim(t.labels) eq label)
      use_ref = where(strtrim(t.labels) eq 'full')
   
      x = struct_var(t, 'H_DELTA')
      hdhb_my[i] = x[use]
      x = struct_var(t, 'DH_DELTA')
      dhdhb_my[i] = x[use]
      x = struct_var(t, 'H_GAMMA')
      hghb_my[i] = x[use]
      x = struct_var(t, 'DH_GAMMA')
      dhghb_my[i] = x[use]

      x = struct_var(t, 'H_DELTA')
      hdhb_my_ref[i] = x[use_ref]
      x = struct_var(t, 'DH_DELTA')
      dhdhb_my_ref[i] = x[use_ref]
      x = struct_var(t, 'H_GAMMA')
      hghb_my_ref[i] = x[use_ref]
      x = struct_var(t, 'DH_GAMMA')
      dhghb_my_ref[i] = x[use_ref]
      
   endfor
   
   print, "I found Hd/Hb=", hdhb_my
   print, "I found Hg/Hb=", hghb_my


   if keyword_set(ps) then begin
      ps_on, dd.figdir+'balmer_ratio_versus_curti.ps', aspect=0.9, xlen=12, $
             /times, /color
      !p.font = 0
      !x.thick = 2
      !y.thick = 2
      cz = 1.6
   endif else begin
      cz = 1.5
   endelse
   
   plot, [0, 1], [0, 1], /nodata, xrange=[0, 0.7], $
         yrange=[0, 0.9], /xs, /ys, $
         xtitle='Balmer ratio (Curti et al 2022)', $
         ytitle='Balmer ratio (renorm L3)', $
         charsize=cz
   abline, 0, 1, linestyle=2, color=jb_colour('gray')


   hd_colour = 'CornFlowerBlue'
   hg_colour = 'orange'
   
   symbols, 2, 2
   x = hd_curti & dx = dhd_curti
   y = hdhb_my & dy = dhdhb_my
   myoploterr2, x, y, x-dx, x+dx, y-dy, y+dy, psym=8, $
                color=jb_colour(hd_colour), errcolor=jb_colour(hd_colour)

   symbols, 1, 2
   x = hd_curti & dx = dhd_curti
   y = hdhb_my_ref & dy = dhdhb_my_ref
   myoploterr2, x, y, x-dx, x+dx, y-dy, y+dy, psym=8, $
                color=jb_colour(hd_colour), errcolor=jb_colour(hd_colour)

   
   symbols, 30, 0.7
   x = hg_curti & dx = dhg_curti
   y = hghb_my & dy = dhghb_my
   myoploterr2, x, y, x-dx, x+dx, y-dy, y+dy, psym=8, $
                color=jb_colour(hg_colour), errcolor=jb_colour(hg_colour)


   symbols, 20, 0.7
   x = hg_curti & dx = dhg_curti
   y = hghb_my_ref & dy = dhghb_my_ref
   myoploterr2, x, y, x-dx, x+dx, y-dy, y+dy, psym=8, $
                color=jb_colour(hg_colour), errcolor=jb_colour(hg_colour)


   xyouts, 0.2, 0.05, TeXtoIDL('H\delta/H\beta'), color=jb_colour('blue'), $
           align=0.5, charsize=cz
   
   xyouts, 0.44, 0.05, TeXtoIDL('H\gamma/H\beta'), color=jb_colour('DarkOrange'), $
           align=0.5, charsize=cz


   ytop = 0.85
   xstart = 0.03
   dy = 0.05
   ddy = -0.01
   dx = 0.03
   colors = [hd_colour, hd_colour, hg_colour, hg_colour]
   symbols = [1, 2, 20, 30]
   symsize = [2, 2, 0.7, 0.7]
   text = ['L3 unmodified', 'L3 renormalised', 'L3 unmodified', $
           'L3 renormalised']
   for i=0L, n_elements(symbols)-1 do begin
      symbols, symbols[i], symsize[i]
      plots, xstart, ytop-i*dy, psym=8, symsize=1.5, $
             color=jb_colour(colors[i])
      xyouts, xstart+dx, ytop-i*dy+ddy, text[i], charsize=cz
   endfor

   
   if keyword_set(ps) then begin
      ps_off
      my_cleanplot, /silent
   endif
end
