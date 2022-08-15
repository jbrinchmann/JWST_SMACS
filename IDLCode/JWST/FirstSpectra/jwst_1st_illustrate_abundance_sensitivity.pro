pro jwst_1st_illustrate_abundance_sensitivity, grid, ps=ps, $
 powerlaw=powerlaw
   ;;
   ;;
   ;; Create a simple illustration of the sensitivity of the abundnace
   ;; estimator to slope changes/ 
   ;;
   ;; We do this by changing TauV using a \lambda^(-1.3) dust-law but
   ;; the same conclusions hold for a scaling change etc. 
   ;;

   common pyneb_te, pn
   
   ;; First get a few model predictions. The way I do this is to focus
   ;; on temperatures above 10^4 K. Then for each of these I create a
   ;; few abundances and predict [O III]5007 and [O II]3727 fluxes
   ;; from this. 


   ;; Find closest models to a grid
   T_e_wanted = mkarr(1.0, 2.5, 0.1)*1e4
   use = lonarr(n_elements(T_e_wanted))
   T_e_used = T_e_wanted
   for i=0L, n_elements(T_e_wanted)-1 do begin
      dum = min(abs(T_e_wanted[i]-pn.Te), i_close)
      T_e_used[i] = pn[i_close].Te
      use[i] = i_close
   endfor
   n_Te = n_elements(use)

;   stop
   ;; Calculate emissivities for crucial lines
   ;;
   ;;  [O II]3727, [O III]4363, [O III]5007
   ;;  Hd, Hg, Hb
   emiss = dblarr(6, n_Te)
   emiss[0, *] = pn[use].o2_3726+pn[use].o2_3729
   emiss[1, *] = pn[use].o3_4363
   emiss[2, *] = pn[use].o3_5007
   emiss[3, *] = pn[use].h1_4101
   emiss[4, *] = pn[use].h1_4340
   emiss[5, *] = pn[use].h1_4861
   

   ;; The wavelengths of the lines are needed for dust attenuation 
   l = [3727.0, 4363, 5007, 4101, 4340, 4861]
   l_norm = 4700.0
   i_4363 = where(l eq 4363)
   i_5007 = where(l eq 5007)
   i_3727 = where(l eq 3727)
   i_hd = where(l eq 4101)
   i_hg = where(l eq 4340)
   i_hb = where(l eq 4861)
   
   ;;
   ;; Line fluxes - scaled by abundances
   ;;
   ;; First the 12+Log O/H values we want to consider
   OH12 = [7.0, 7.25, 7.5, 7.75, 8.0, 8.25, 8.5]
   ;; These are the actual total abundances
   OH = 10.0^(OH12-12.0)

   ;; Next, we need to divide this into ionic abundances. There is no
   ;; way to do this very simply. I will use the results from
   ;; Brinchmann et al (2008) and just take the ratio of n(O+)/n(O++)
   ;; vs T(O++) up to a temperature of 17kK and then constant after
   ;; that.
   ;;
   ;;   t = mrdfits('/Users/jarle/Work/Wolf-Rayet/DR7/Abundances/wr_abundances_dr7_table.fits', 1)
   ;;   ii = where(sc.oiii4363/sc.doiii4363 gt 7 and class.i_class eq
   ;;   1)
   ;;   sixlin, t[ii].te_o3, t[ii].o2h/t[ii].o3h, a, siga, b, sigb
   ;;   IDL> print, a[0], b[0]
   ;;    0.85234140     -0.39206624
   a = 0.85234140 & b = -0.39206624
   r_no2_no3 = a + b*T_e_used/1e4
   above = where(T_e_used gt 17e3, n_above)
   if n_above gt 0 then $
    r_no2_no3[above] = a+b*1.7

   n_OH = n_elements(OH)

   OH_OII = dblarr(n_Te, n_OH)
   OH_OIII = OH_OII
   
   for i=0L, n_OH-1 do begin
      OH_OII[*, i] = (r_no2_no3/(1+r_no2_no3))*OH[i]
      OH_OIII[*, i] = (1.0/(1+r_no2_no3))*OH[i]
   endfor

   
   
   ;; Create a 2D grid

   fluxes = dblarr(6, n_Te, n_OH)
   for i=0L, n_OH-1 do begin
      ;; We need to scale the O fluxes
      fluxes[0, *, i] = emiss[0, *]*OH_OII[*, i]/emiss[5, *]
      fluxes[1, *, i] = emiss[1, *]*OH_OIII[*, i]/emiss[5, *]
      fluxes[2, *, i] = emiss[2, *]*OH_OIII[*, i]/emiss[5, *]
      ;; But we leave H lines untouched
      fluxes[3, *, i] = emiss[3, *]/emiss[5, *]
      fluxes[4, *, i] = emiss[4, *]/emiss[5, *]
      fluxes[5, *, i] = emiss[5, *]/emiss[5, *]
   endfor

   
   ;;
   ;; Loop over dust attenuation and calculate abundances.
   ;;
   if keyword_set(powerlaw) then begin
      ;; Adjust with a power-law slope
      alpha = mkarr(-4, 4, 0.1)
      n_slope = n_elements(alpha)
   endif else begin
      tV = mkarr(0, 4.0, 0.1)
      n_slope = n_elements(tV)
   endelse

   
   n_models = n_Te
   T_e_standard = dblarr(n_slope, n_models, n_OH)
   T_e_double= dblarr(n_slope, n_models, n_OH)


   if n_elements(grid) eq 0 then begin
   
      OH2_s = dblarr(n_slope, n_models, n_OH)
      OH2_d = dblarr(n_slope, n_models, n_OH)
      OH3_s = dblarr(n_slope, n_models, n_OH)
      OH3_d = dblarr(n_slope, n_models, n_OH)
      
      
      for i_mod=0L, n_models-1 do begin
         for i=0L, n_slope-1 do begin
            ;;
            ;; Attenuate the emiss
            ;;
            if keyword_set(powerlaw) then begin
               scale = (l/l_norm)^alpha[i]
            endif else begin
               tau = tV[i]*(l/l_norm)^(-1.3)
            endelse
            
            for iOH=0L, n_OH-1 do begin

               if keyword_set(powerlaw) then begin
                  fobs = fluxes[*, i_mod, iOH]*scale
               endif else begin
                  fobs = fluxes[*, i_mod, iOH]*exp(-tau)
               endelse
               
               ;; First, estimate temperature
               
               ;; 4363/5007
               Te_ratio_standard = alog10(fobs[i_4363]/fobs[i_5007])
               ;; 4363/Hg/5007/Hb
               Te_ratio_double = alog10((fobs[i_4363]/fobs[i_hg])/(fobs[i_5007]/fobs[i_hb]))
               
               T_e_standard[i, i_mod, iOH] = te_pyneb_o3(Te_ratio_standard, what='standard')
               T_e_double[i, i_mod, iOH] = te_pyneb_o3(Te_ratio_double, what='doubleratio')


               ;; Then, estimate the ionic abundances
               ;;
               ;; In the standard approach we normalise all to Hb.
               r_o2hb = fobs[i_3727]/fobs[i_hb]
               r_o3hb = fobs[i_5007]/fobs[i_hb]
               OH2_s[i, i_mod, iOH] = o_ionic_abundance_pyneb(r_o2hb, ['3727', 'Hb'], $
                                                              T_e_standard[i, i_mod, iOH])
               OH3_s[i, i_mod, iOH] = o_ionic_abundance_pyneb(r_o3hb, ['5007', 'Hb'], $
                                                              T_e_standard[i, i_mod, iOH])


               ;; In the non-standard approach we use the double-ratio
               ;; temperature and normalise O2 to Hd
               
               r_o2hd = fobs[i_3727]/fobs[i_hd]
               OH2_d[i, i_mod, iOH] = o_ionic_abundance_pyneb(r_o2hd, ['3727', 'Hd'], $
                                                              T_e_double[i, i_mod, iOH])
               OH3_d[i, i_mod, iOH] = o_ionic_abundance_pyneb(r_o3hb, ['5007', 'Hb'], $
                                                              T_e_double[i, i_mod, iOH])

            endfor
         endfor
      endfor
      
      OH_d = OH2_d+OH3_d
      OH_s = OH2_s+OH3_s

      OH12_d = 12+alog10(OH_d)
      OH12_s = 12+alog10(OH_s)
      
      grid = {fluxes: fluxes, OH2_s: OH2_s, OH2_d: OH2_d, OH3_s: OH3_s, $
              OH3_d: OH3_d, T_e_s: T_e_standard, T_e_d: T_e_double, $
              OH_d: OH_d, OH_s: OH_s, OH12_d: OH12_d, OH12_s: OH12_s}
      
      if keyword_set(powerlaw) then begin
         grid = create_struct(grid, 'type', 'powerlaw', 'alpha', alpha)
      endif else begin
         grid = create_struct(grid, 'type', 'tau', 'tauV', tV)
      endelse
      
      
   endif
   
      
      

   ;;
   ;; Show trends. We use percentage change in temperature
   ;;

   jwst_1st_dirs, dd
   if keyword_set(ps) then begin
      ps_on, dd.figdir+'abundance_sensitivity_plot.ps', /color, aspect=0.6, xlen=12, /times
      !p.font = 0
      !x.thick = 2
      !y.thick = 2
      tk = 3
      cz = 1.5
   endif else begin
      tk = 1
      cz = 1.2
   endelse

   if grid.type eq 'powerlaw' then begin
      x = grid.alpha
      xrange = [-4, 4]
      xtitle = 'Power-law index'
      yrange = [-0.3, 0.3]
   endif else begin
      x = grid.tauV
      xrange = [0, 4.5]
      xtitle = '\tau_V'
      yrange = [0.0, 0.45]
   endelse

   
   loadct, 0
   plot, [-4, 4], [-0.0, 0.5], /nodata, xtitle=TeXtoIDL(xtitle), $
         ytitle=TeXtoIDL('OH-OH_{true} [dex]'), charsize=cz, $
         xrange=xrange, yrange=yrange, /ys

      
   c = model_grid_colors(n_grid=n_models+2, table=3, /reverse)
   for i_mod=0L, n_models-1 do begin
      for i_OH=0L, n_OH-1, 3 do begin
         y1 = (grid.OH12_d[*, i_mod, i_OH]-OH12[i_OH])
         y2 = (grid.OH12_s[*, i_mod, i_OH]-OH12[i_OH])

         oplot, x, y1, color=c.g[i_mod+1], thick=tk, linestyle=2
         oplot, x, y2, color=c.g[i_mod+1], thick=tk


;         stop
         
         Te_true = grid.T_e_d[0, i_mod, i_OH]
         if (Te_true lt 1e4) then $
          label = string(format='("T_e=",F3.1,"x10^3 K")', Te_true/1e3) $
         else $
          label = string(format='("T_e=",F3.1,"x10^4 K")', Te_true/1e4)
      
         xyouts, max(x)+0.05, y2[n_slope-1]-0.01, TeXtoIDL(label), charsize=cz, $
                 color=jb_colour('gray')
;         stop
      endfor
      
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
