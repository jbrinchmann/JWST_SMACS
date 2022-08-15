pro jwst_1st_ciii_in_4590

   flx = jwst_1st_load_flux_overview(4590)
   t = jwst_1st_load_one_gaussfit('02736_04590')

   ;;
   ;; In this case I'll use the non-modified spectrum
   ;;
   ;; I did:
   ;;
   ;; d = jwst_1st_get_coadded(4590)
   ;; jwst_1st_fit_gaussian_interactive, 1e4*d.oned.wavelength,
   ;; d.oned.flux, d.oned.dflux, 8.4951057, res=res, LSF_R=d.oned.R
   ;; save, res,
   ;; file='/data2/jarle/JWST/ERO/SMACS/Spectra/GaussianFits/ciii-4590.sav'

   restore, file='/data2/jarle/JWST/ERO/SMACS/Spectra/GaussianFits/ciii-4590.sav'

   x1 = res.ciii_1907.flux+res.ciii_1909.flux
   dx1 = sqrt(res.ciii_1907.dflux^2+res.ciii_1909.dflux^2)
   x2 = res.h_beta.flux
   x3 = res.oiii_5007.flux
   dx2 = res.h_beta.dflux
   dx3 = res.oiii_5007.dflux

   xr = (x1+randomn(sss, 5000)*dx1)/(x2 + randomn(sss, 5000)*dx2)
   yr = (x3+randomn(sss, 5000)*dx3)/(x2 + randomn(sss, 5000)*dx2)

   v_lin= quantile(xr, [0.16, 0.5, 0.84])
   print, 'X=', v_lin[1], v_lin[1]-v_lin[0], v_lin[2]-v_lin[0]

   v_lin= quantile(yr, [0.16, 0.5, 0.84])
   print, 'Y=', v_lin[1], v_lin[1]-v_lin[0], v_lin[2]-v_lin[0]

   ;; Need to get this for the other measurements too.

   ;; This is from the Gaussian fits on the renormalised spectrum
   ratio_x_rn = 10.0^(-0.294353)
   ratio_y_rn = 10.0^0.520498
   ratio_x_rn_low = 10.0^(-0.294353-0.131854)
   ratio_x_rn_high = 10.0^(-0.294353+0.116223)
   ratio_y_rn_low = 10.0^(0.520498-0.059615)
   ratio_y_rn_high = 10.0^(0.520498+0.0647938)

   print, 'X=', ratio_x_rn, ratio_x_rn-ratio_x_rn_low, ratio_x_rn_high-ratio_x_rn
   print, 'Y=', ratio_y_rn, ratio_y_rn-ratio_y_rn_low, ratio_y_rn_high-ratio_y_rn
   
   
end



pro compare_to_gutkin_feltre, ps=ps

   jwst_1st_dirs, dd

   f16 = mrdfits('/data2/jarle/PIModels/data/NEOGAL/AGN_NLR_nebular_feltre16/Feltre16-models.fits', 1)
   g16 = neogal_g16_original(/extended)

   neogal_g16_dim_to_array, g16.dim, adim

   x_g16 = alog10(g16.ciii1908/g16.hb)
   y_g16 = alog10(g16.oiii5007/g16.hb)
   
   f16_ciii1908 = f16.ciii1907+f16.ciii1910
   x_f16 = alog10(f16_ciii1908/f16.hb)
   y_f16 = alog10(f16.oiii5007/f16.hb)

   whist2d, x_g16, y_g16, h_g16, [-2, 1.3], [0, 1.5], 20, 20, weight=adim.logu, /average
   whist2d, x_f16, y_f16, h_f16, [-2, 1.3], [0, 1.5],20, 20, weight=f16.logu, /average

   ;; Values from the code above
   
   ;; These are unnorm ,norm
   xpos = [1.65987, 0.507747]
   xlow = xpos-[0.454874,  0.132952]
   xhi =  xpos+[0.948323, 0.155798]

   ypos = [ 3.21955,  3.31511]
   ylow = ypos-[0.163425, 0.425209]
   yhi = ypos+[0.341555, 0.533393]

   ;;
   ;; I also want a reddening vector
   ;;
   slope = -1.3
   tau_V = 1.0
   l_V = 5500.0

   dtau_V_x = tau_V*((1909.0/l_V)^slope-(4861.0/l_V)^slope)
   d_x = alog10(exp(1.0))*dtau_V_x

   dtau_V_y = tau_V*((5007.0/l_V)^slope-(4861.0/l_V)^slope)
   d_y = alog10(exp(1.0))*dtau_V_y

   
;   stop
   
   if keyword_set(ps) then begin
      ps_on, dd.figdir+'CIII_vs_G16_F16.ps', aspect=0.5, xlen=12, /times, $
             /color
      !p.font = 0
      !x.thick = 2
      !y.thick = 2
      cz = 1.5
   endif
   get_ct, 'Spectral'

   pos = get_position_arr(0, nx=2, ny=1, xgap=0.01, xmax=0.9)
   tvim_struct, h_g16, position=pos, range=[-4.5, -1], $
                xtitle=TeXtoIDL('Log C III]1907,1909/H\beta'), $
                ytitle=TeXtoIDL('Log [O III]5007/H\beta'), charsize=cz, nbot=3
   symbols, 32, 0.8
   myoploterr2, alog10(xpos[0]), alog10(ypos[0]), alog10(xlow[0]), $
                alog10(xhi[0]), alog10(ylow[0]), alog10(yhi[0]), psym=8
   
   symbols, 30, 0.8
   myoploterr2, alog10(xpos[1]), alog10(ypos[1]), alog10(xlow[1]), $
                alog10(xhi[1]), alog10(ylow[1]), alog10(yhi[1]), psym=8

   jb_text, -0.02, 0.05, 'Gutkin et al (2016)', align=1, /frac, /relative, $
            charsize=cz

   arr_col = 'gray50'
   cgArrow, -1, 1.35,-1+d_x, 1.35+d_y, /data, thick=2, /solid, $
            color=jb_colour(arr_col)
   xyouts, -1+d_x/2., 1.35+0.02, TeXtoIDL('\tau_V=1'), color=jb_colour(arr_col), $
           charsize=cz, align=0.3
   

   get_ct, 'Spectral'
   pos = get_position_arr(1, tickf=tf)
   tvim_struct, h_f16, position=pos, range=[-4.5, -1], /noerase, /scale, $
                stitle='Log U', ytickformat=tf[1], $
                xtitle=TeXtoIDL('Log C III]1907,1909/H\beta'), $
                charsize=cz, nbot=3
   symbols, 32, 0.8
   myoploterr2, alog10(xpos[0]), alog10(ypos[0]), alog10(xlow[0]), $
                alog10(xhi[0]), alog10(ylow[0]), alog10(yhi[0]), psym=8
   
   symbols, 30, 0.8
   myoploterr2, alog10(xpos[1]), alog10(ypos[1]), alog10(xlow[1]), $
                alog10(xhi[1]), alog10(ylow[1]), alog10(yhi[1]), psym=8
   

   jb_text, -0.02, 0.05, 'Feltre et al (2016)', align=1, /frac, /relative, $
            charsize=cz
   
   symbols, 32, 0.6
   xx = -1.8
   dx = 0.07
   yy = 1.43
   ddy = -0.02
   dy = 0.09
   plots, xx, yy, psym=8
   xyouts, xx+dx, yy+ddy, 'Pipeline'

   symbols, 30, 0.6
   plots, xx, yy-dy, psym=8
   xyouts, xx+dx, yy+ddy-dy, 'Re-renormalised'

   
;   stop
   
   if keyword_set(ps) then begin
      ps_off
      my_cleanplot, /silent
   endif
   
end
