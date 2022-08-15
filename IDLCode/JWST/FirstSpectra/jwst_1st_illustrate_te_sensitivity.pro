pro jwst_1st_illustrate_te_sensitivity, ps=ps
   ;;
   ;;
   ;; Create a simple illustration of T_e sensitivity
   ;;
   ;; We do this by changing TauV using a \lambda^(-1.3) dust-law but
   ;; the same conclusions hold for a scaling change etc. 
   ;;

   common pyneb_te, pn

   ;; First get a few model predictions.

   use = lindgen(5)*(500L/5)        ; from 0 to 500 - the pyneb grid has 500 elements

   fluxes = fltarr(4, n_elements(use))
   fluxes[0, *] = pn[use].o3_4363
   fluxes[1, *] = pn[use].o3_5007
   fluxes[2, *] = pn[use].h1_4340
   fluxes[3, *] = pn[use].h1_4861
   

   l = [4363, 5007, 4340, 4861]
   l_norm = 4700.0


   ;;
   ;; Dust attenuation grid
   ;;
   tV = mkarr(0.0, 4.0, 0.01)
   n_tV = n_elements(tV)
   n_models = n_elements(use)
   T_e_standard = dblarr(n_tV, n_models)
   T_e_double= dblarr(n_tV, n_models)


   for i_mod=0L, n_models-1 do begin
      for i=0L, n_tV-1 do begin
         ;;
         ;; Attenuate the fluxes
         ;; 
         tau = tV[i]*(l/l_norm)^(-1.3)
         fobs = fluxes[*, i_mod]*exp(-tau)
         
         ratio_standard = alog10(fobs[0]/fobs[1])
         ratio_double = alog10((fobs[0]/fobs[2])/(fobs[1]/fobs[3]))
         
         T_e_standard[i, i_mod] = te_pyneb_o3(ratio_standard, what='standard')
         T_e_double[i, i_mod] = te_pyneb_o3(ratio_double, what='doubleratio')
         
      endfor
   endfor


   ;;
   ;; Show trends. We use percentage change in temperature
   ;;

   jwst_1st_dirs, dd
   if keyword_set(ps) then begin
      ps_on, dd.figdir+'Te_sensitivity_plot.ps', /color, aspect=0.6, xlen=12, /times
      !p.font = 0
      !x.thick = 2
      !y.thick = 2
      tk = 5
      cz = 1.5
   endif

   loadct, 0
   plot, [0, 5], [-40, 0], /nodata, xtitle=TeXtoIDL('\tau_V'), $
         ytitle=TeXtoIDL('100 x (T_e-T_{e,true})/T_{e,true}'), charsize=cz
   
   c = model_grid_colors(n_grid=n_models+2, table=33, /reverse)
   for i_mod=0L, n_models-1 do begin

      y1 = 100*(T_e_standard[*, i_mod]-T_e_standard[0, i_mod])/T_e_standard[0, i_mod]
      y2 = 100*(T_e_double[*, i_mod]-T_e_double[0, i_mod])/T_e_double[0, i_mod]

      oplot, tV, y1, color=c.g[i_mod+1], thick=tk
      oplot, tV, y2, color=c.g[i_mod+1], linestyle=2, thick=tk

      Te_true = T_e_double[0, i_mod]
      if (Te_true lt 1e4) then $
       label = string(format='("T_e=",F3.1,"x10^3 K")', Te_true/1e3) $
      else $
       label = string(format='("T_e=",F3.1,"x10^4 K")', Te_true/1e4)
      
      xyouts, max(tV)+0.05, y1[n_tV-1]-0.3, TeXtoIDL(label), charsize=cz
      
   endfor

   dy = -3
   ddy = -0.4
   ypos = -35
   loadct, 0
   oplot, [0.2, 0.5], [ypos, ypos], thick=tk
   xyouts, 0.55, ypos+ddy, charsize=cz, $
           TeXtoIDL('[O III]4363/[O III]5007')
   oplot, [0.2, 0.5], [ypos+dy, ypos+dy], thick=tk, linestyle=2
   xyouts, 0.55, ypos+dy+ddy, charsize=cz, $
           TeXtoIDL('[O III]4363/H\gamma/[O III]5007/H\beta')
   
   
   
   if keyword_set(ps) then begin
      ps_off
      my_cleanplot, /silent
   endif
   
end
