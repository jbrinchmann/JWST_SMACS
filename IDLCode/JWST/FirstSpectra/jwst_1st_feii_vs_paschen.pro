pro _load_models, mshock, magn, msb
   ;;
   ;; Load three illustrative models
   ;;
   jwst_1st_dirs, dd
   DIR = dd.root+'IteraModels/'

   
   magn = load_exported_itera_model(DIR+'feii-itera_grid_Groves2004_dusty_Z1_n1e3.txt')
   mshock = load_exported_itera_model(DIR+'feii-itera_grid_Allen2008_Solar_n0p01.txt')
   msb = load_exported_itera_model(DIR+'feii-itera_grid_Levesque2010High_Z0p2_n100.txt')

   
end
   

   

pro jwst_1st_feii_vs_paschen, ps=ps
   ;;
   ;; Show the FeII lines versus Paschen lines for the sample,
   ;; compared to the Izotov & Thuan (2016) results from BCDs.
   ;;
   ;; I will also compare to 
   
   jwst_1st_dirs, dd
   
   ;;
   ;; Get the low-z comparison sample.
   ;;
   cat_IT2016_file = dd.root+'Lowz-Catalogues/IzotovThuan2016/IT2016_FeII.fits'
   IT2016 = mrdfits(cat_IT2016_file, 1)


   fe_pab_IT2016 = IT2016.feii1p257_hb_p*(IT2016.hb/100)/IT2016.pb

   ;; We want to normalise Fe II 1.257 to Pa-Beta as they are close
   ;; And FeII 1.644 to ?? We'll keep Pa-beta for this too..
   

   ;; Those objects showing Fe II in emission
   todo = ['03042', '09239', '09483']
   upper_limit = [0, 0, 0]

   ;; 
   ;; These show upper limits. I'll take the Pb error estimate
   ;; to estimate upper limits
   todo = [todo, '01917', '05735', '08506', '09721', '09922']
   upper_limit = [upper_limit, 1,1,1,1,1]

   ;; 
   feii_1p257_pab_SMACS = []
   feii_1p644_pab_SMACS = []
   hei_pag = []
   dhei_pag = []
   dfeii_1p257_pab_SMACS = []
   dfeii_1p644_pab_SMACS = []
   is_upper = []
   is_left = []
   for i=0L, n_elements(todo)-1 do begin
      g = jwst_1st_load_one_gaussfit('02736_'+todo[i])

      ii = where(g.line eq 'FEII_1P257', n_ii)
      ii_pb = where(g.line eq 'HI_3_5', n_ii_pb)

      if upper_limit[i] then begin
         ;; Then we do not care about FeII_1P257 so this is ignored
         if n_ii_pb eq 0 then continue ; No Paschen Beta so ignored.

         f_upper = 1.0*g.dflux[ii_pb]
         r = f_upper/g.flux[ii_pb]
         feii_1p257_pab_SMACS = [feii_1p257_pab_SMACS, r]
         dfeii_1p257_pab_SMACS = [dfeii_1p257_pab_SMACS, r*0]
         is_upper = [is_upper, 1]
      endif else begin
         if (n_ii+n_ii_pb gt 1) then begin
            ;; We have both lines we want
            r = g.flux[ii]/g.flux[ii_pb]
            feii_1p257_pab_SMACS = [feii_1p257_pab_SMACS,r ]
            
            dr = r*error_on_fraction(g.flux[ii], g.dflux[ii], g.flux[ii_pb], g.dflux[ii_pb])
            dfeii_1p257_pab_SMACS = [dfeii_1p257_pab_SMACS, dr]
         endif
         is_upper = [is_upper, 0]
      endelse
         
      ;; All of these have He I 1.083 as well so we can add this safely.
      ii_hei = where(g.line eq 'HEI_1P083', n_hei)
      ii_pag = where(g.line eq 'HI_3_6', n_pag)

      this_is_left = 0
      if (n_hei+n_pag ne 2) then begin
         ;; Missing?
         if (n_hei+n_pag eq 0) then begin
            this_r = !values.f_nan 
         endif else begin
            if n_pag eq 0 then begin
               f_hei = g.flux[ii_hei]
               if n_ii_pb gt 0 then begin
                  ;; Use Case B to estimate Pa-g
                  f_pg = g.flux[ii_pb]*0.554
               endif else  stop
            endif else begin
               f_pg = g.flux[ii_pag]
               this_is_left = 1
               f_hei = g.dflux[ii_pag]
            endelse

            this_r = f_hei/f_pg
            
         endelse

         is_left = [is_left, this_is_left]
         hei_pag = [hei_pag, this_r]
         dhei_pag = [dhei_pag, this_r*0]
         
      endif else begin
         r = g.flux[ii_hei]/g.flux[ii_pag]
         dr = r*error_on_fraction(g.flux[ii_hei], g.dflux[ii_hei], g.flux[ii_pag], $
                                  g.dflux[ii_pag])
         
         hei_pag = [hei_pag, r]
         dhei_pag = [dhei_pag, dr]
         is_left = [is_left, 0]
      endelse
      
      ii = where(g.line eq 'FEII_1P644', n_ii)
      ii_pb = where(g.line eq 'HI_3_5', n_ii_pb)
      if (n_ii+n_ii_pb gt 1) then begin
         r = g.flux[ii]/g.flux[ii_pb]
         feii_1p644_pab_SMACS = [feii_1p644_pab_SMACS,r ]

         dr = r*error_on_fraction(g.flux[ii], g.dflux[ii], g.flux[ii_pb], g.dflux[ii_pb])
         dfeii_1p644_pab_SMACS = [feii_1p644_pab_SMACS, dr]
      endif
   endfor

   ;;
   ;; I want some easier variables for the below.
   ;;
   rx = feii_1p257_pab_SMACS
   ry = hei_pag
   x = alog10(rx)
   y = alog10(ry)
   drx = dfeii_1p257_pab_SMACS
   dry = dhei_pag
   xlow=alog10(rx-drx)
   xhi = alog10(rx+drx)
   ylow=alog10(ry-dry)
   yhi=alog10(ry+dry)


   i_det = where(is_upper eq 0, compl=i_nondet)
   
   ;;
   ;; Load relevant models.
   ;;
   _load_models, mshock, magn, msb

   
   

   if keyword_set(ps) then begin
      ;; I will edit this in Illustrator
      ps_on, dd.figdir+'feii_hei_dd.ps', aspect=0.8, xlen=12, /color, /times
      !p.font = 0
      !x.thick = 2
      !y.thick = 2
      cz = 1.6
   endif else cz = 1.2
   
   ;; Show [Fe II]1.257/Pa-b vs He I 1.083/Pa-gamma diagram


   plot, [0, 1], [0, 1], /nodata, xrange=[-3, 1], yrange=[-1, 1], charsize=cz, $
         xtitle=TeXtoIDL('Log [Fe II] 1.257/Pa-\beta'), $
         ytitle=TeXtoIDL('Log He I 1.083/Pa-\gamma'), xstyle=8
   
   jwst_1st_feii_hei_model_grid, mshock, /overplot
   jwst_1st_feii_hei_model_grid, magn, /overplot, color1='orange', color2='Aquamarine'
   jwst_1st_feii_hei_model_grid, msb, /overplot, color1='firebrick', color2='ForestGreen'

   symbols, 30, 1
   myoploterr2, x[i_det], y[i_det], xlow[i_det], xhi[i_det], ylow[i_det], yhi[i_det], psym=8


;   symbols, 32, 1
;   myoploterr2, x[i_nondet], y[i_nondet], xlow[i_nondet], $
;                xhi[i_nondet], ylow[i_nondet], yhi[i_nondet], psym=8
   symbols, 12, 1
   myoploterr2, x[i_nondet], y[i_nondet], xlow[i_nondet], $
                xhi[i_nondet], ylow[i_nondet], yhi[i_nondet], psym=8
   
   ;;
   ;; Make a second x-axis with L(Fe II)/SFR. To do this, we need a
   ;; shifted x-axis. 
   ;;
   ;; To do this I will use that with a Kroupa IMF and a Kennicutt
   ;; (1998) calibration from Ha to SFR and assuming no dust, we have
   ;;
   ;;    SFR = 2.085*(L_Pab/1d40) -> L_Pab = (1e40/2.085)SFR
   ;;                                      = 4.7943e39 SFR
   ;;
   ;; So we have 
   ;;
   ;;    L_Fe/L_Pb = L_Fe/(4.7943d39*SFR) => L_Fe/SFR = 4.7943d39
   ;;                                                  L_Fe/L_Pb
   AXIS, xAXIS=1, xRANGE = !x.CRANGE+alog10(4.7943d39), xSTYLE = 1, $
         xTITLE = TeXtoIDL('L([Fe II]1.257)/SFR  [erg/s/M_{sun}/yr]'), charsize=cz
   
   if keyword_set(ps) then begin
      ps_off
      ps_on, dd.figdir+'feii_histograms.ps', aspect=0.6, xlen=12, /color, /times
   endif  else stop
   

   
   ;; Next, show histograms of the [Fe II]/Pa-b ratios
   range = [-3, 1]
   nbins = 25
   whist1d_alt, alog10(fe_pab_IT2016), h_IT, range, nbins
   whist1d_alt, alog10(mshock.FEII_1_2570UM/mshock.HI_1_2818UM), h_shock, range, nbins
   whist1d_alt, alog10(magn.FEII_1_2570UM/magn.HI_1_2818UM), h_agn, range, nbins
   whist1d_alt, alog10(msb.FEII_1_2570UM/msb.HI_1_2818UM), h_sb, range, nbins

   ;; Alternatives: densities
   d_shock = r_density(alog10(mshock.FEII_1_2570UM/mshock.HI_1_2818UM))
   d_agn = r_density(alog10(magn.FEII_1_2570UM/magn.HI_1_2818UM))

   vals = alog10(msb.FEII_1_2570UM/msb.HI_1_2818UM)
   ok = where(finite(vals) eq 1)
   vals = vals[ok]
   d_sb = r_density(vals)
;   stop
   
   pos = [0.1, 0.1, 0.95, 0.95]
   plot, [-3, 1], [0, 1.2], /nodata,  yrange=[0, 1.1], xrange=[-3, 1], $
         xtitle=TeXtoIDL('Log [Fe II]1.257/Pa-\beta'), $
         charsize=cz, position=pos
   plot_hstruct, h_IT, peak=1, /overplot, /dofill, $
                 fillcolor=jb_colour('gray50'), $
                 position=pos, xstyle=8, ystyle=4
   oplot, d_shock.x, d_shock.y/max(d_shock.y), color=jb_colour('CornFlowerBlue')
   oplot, d_agn.x, d_agn.y/max(d_agn.y), color=jb_colour('orange')
   oplot, d_sb.x, d_sb.y/max(d_sb.y), color=jb_colour('ForestGreen')

;   stop
   
   ;; plot_hstruct, h_shock, peak=0.3, /overplot, color=jb_colour('CornFlowerBlue'), $
   ;;               /dofill, fillcolor=jb_colour('CornFlowerBlue'), $
   ;;               position=pos
   ;; plot_hstruct, h_agn, peak=0.3, /overplot, color=jb_colour('aquamarine'), /dofill, $
   ;;               fillcolor=jb_colour('aquamarine'), $
   ;;               position=pos
   ;; plot_hstruct, h_sb, peak=0.3, /overplot, color=jb_colour('ForestGreen'), /dofill, $
   ;;               fillcolor=jb_colour('ForestGreen'), $
   ;;               position=pos

   for i=0L, n_elements(x)-1 do begin
      oplot, x[i]+[0, 0], [0, 1], linestyle=2
      oplot, xlow[i]+[0, 0], [0, 1], linestyle=1
      oplot, xhi[i]+[0, 0], [0, 1], linestyle=1
      xyouts, x[i], 1.05, todo[i], align=0.5, charsize=cz
   endfor
   
   if keyword_set(ps) then begin
      ps_off
      my_cleanplot, /silent
   endif
   
   
end
