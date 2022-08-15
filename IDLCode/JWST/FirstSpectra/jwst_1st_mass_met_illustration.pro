function _make_one_rm, x, dx, y, dy, subset, xout, boot=boot, n_boot=n_boot, $
 outfile=outfile, force=force

   if n_elements(n_boot) eq 0 then n_boot = 1001

   
   if (n_elements(subset) gt 10000L) then $
    n_per_bin = 1501 $
   else begin
      ;; We want at least 4 bins
      n_per_bin = 101
      if n_elements(subset)/float(n_per_bin) lt 4 then $
       n_per_bin = floor(n_elements(subset)/4.0+0.5)
   endelse
      
   
   if file_test(outfile) and not keyword_set(force) then $
    rr = hdf5_to_struct(outfile) $
   else begin
      if keyword_set(boot) then begin
         rr = boot_running_median(x[subset], y[subset], xout=xout, $
                                  n_per_bin=n_per_bin, n_boot=n_boot, dx=dx[subset], $
                                  dy=dy[subset], min=6, max=12, /new)
      endif else begin
         rr = running_median(x[subset], y[subset], min=6, $
                             max=12, n_per=n_per_bin)
      endelse

      if size(rr, /tname) ne 'STRUCT' then stop
      
      struct_to_hdf5, rr, outfile
   endelse
   
   return, rr

end


function _txt, x

   if x lt 0 then $
    txt = string(format='("m",F3.1)', abs(x)) $
   else $
    txt = string(format='("p",F3.1)', x) 
 
   return, txt
end

function _txt_limit, x1, x2
   return, _txt(x1)+'_'+_txt(x2)
end

pro _calculate_boot_running_medians, s, lgm, class, oh, ssfr, sfr
   ;;
   ;; This routine does the heavy lifting for the plotting routine
   ;; below. 
   ;;


   zcut = 0.01
   
   if (n_elements(oh) eq 0) then $
    oh = dr_load_x('fib-oh', /dr7, /steep)
   if (n_elements(ssfr) eq 0) then $
    ssfr = hdf5_to_struct('/data2/jarle/SDSS/Datafiles/DR7/gal_fibspecsfr_B04_dust_dr7_v5_2.h5')

   if (n_elements(sfr) eq 0) then $
    sfr = hdf5_to_struct('/data2/jarle/SDSS/Datafiles/DR7/gal_fibsfr_B04_dust_dr7_v5_2.h5')
   

   r50 = sdss_quantity('PETROR50', /dr7, /photo)
   diam_r = 2*reform(r50[2, *])
   
   jwst_1st_dirs, dd

   OUTDIR = dd.root+'MassMet/'
   if not file_test(OUTDIR) then $
    mkdir, OUTDIR

   
   
   ;; First without matching in SFR/M* - but requiring the photometric
   ;; diameter is not too larger relative to the size of the galaxies.
   use = where(class.i_class eq 1 and s.zxcor gt zcut and diam_r/3.0 lt 3)

   flag = bytarr(n_elements(s.zxcor))
   flag[use] = 1
   
   use_sn4363 = where(flag eq 1 and s.oiii4363/s.doiii4363 gt 7)

   ;; Now match in sSFR. First the simple way - the 4363 detections
   ;; all have log SFR/M* > -9
   use_ssfr = where(flag eq 1 and ssfr.median ge -9)
   use_sn4363_ssfr = where(flag eq 1 and $
                           s.oiii4363/s.doiii4363 gt 7 and  $
                           ssfr.median ge -9)
   
;   stop
                           
                           
   xout = mkarr(8, 11.4, 0.1)
   x = lgm.lgm
   dx = 0.5*(lgm.lgm_p84-lgm.lgm_p16)
   y = oh.median
   dy = 0.5*(oh.p84-oh.p16)


   outfile_root = ['all-sf', 'all-sn4363sf', 'all-sf-ssfr_gt_m9', $
                   'all-sn4363sf-ssfr_gt_m9']
   for i=0L, n_elements(outfile_root)-1  do begin
      case i of
         0: ii = use
         1: ii = use_sn4363
         2: ii = use_ssfr
         3: ii = use_sn4363_ssfr
      endcase

      outfile = OUTDIR+outfile_root[i]+'.h5'
      print, 'Creating '+outfile+' if needed.'
      rr = _make_one_rm(x, dx, y, dy, ii, xout, outfile = outfile)
      outfile = OUTDIR+outfile_root[i]+'-boot.h5'
      print, 'Creating '+outfile+' if needed.'
      rrb = _make_one_rm(x, dx, y, dy, ii, xout, /boot, outfile = outfile)
   endfor


   ;; Finally, we need to do this in SFR bins
   SFR_bins = [[-2, -1.0], [-1.0, 0], [0, 1], [1, 2]]
   dims = size(SFR_bins, /dimen)
   n_bins = dims[1]
   outfile_root = ['all-sf-sfrbin', 'all-sn4363sf-sfrbin']
   for i=0L, n_bins-1 do begin
      for j=0L, n_elements(outfile_root)-1 do begin
         lo = SFR_bins[0, i]
         hi = SFR_bins[1, i]

         if j eq 0 then begin
            ii = where(flag eq 1 and sfr.median gt lo and $
                       sfr.median le hi, n_ii)
            if (n_ii eq 0) then stop
         endif else begin
            ii = where(flag eq 1 and sfr.median gt lo and $
                       sfr.median le hi and s.oiii4363/s.doiii4363 gt 7, n_ii)
            if (n_ii eq 0) then stop
         endelse

         tail = _txt_limit(lo, hi)
         outfile = OUTDIR+outfile_root[j]+tail+'-boot.h5'
         
         print, 'Creating '+outfile+' if needed'
         rrb = _make_one_rm(x, dx, y, dy, ii, xout, /boot, outfile = outfile)
      endfor
   endfor
end

function _get_one, all=all, ssfr=ssfr, sfr=sfr, sn4363=sn4363, $
                       noboot=noboot


   jwst_1st_dirs, dd

   OUTDIR = dd.root+'MassMet/'

   if keyword_set(all) then $
    root = 'all-sf' $
   else if keyword_set(sn4363) then $
    root = 'all-sn4363sf' $
   else begin
      print, 'You need to set either ALL or SN4363'
      return, -1
   endelse

   if keyword_set(ssfr) then $
    root = root+'-ssfr_gt_m9' $
   else if n_elements(sfr) gt 0 then begin
      tail = _txt_limit(sfr[0], sfr[1])
      root = root+'-sfrbin'+tail
   endif

   outfile = OUTDIR+root

   if keyword_set(noboot) then $
    outfile = outfile+'.h5' $
   else $
    outfile = outfile+'-boot.h5'

   if not file_test(outfile) then begin
      print, outfile+' does not exist!'
      stop
   endif

   
   return, hdf5_to_struct(outfile)
   
   
end

pro jwst_1st_mass_met_illustration, ps=ps
   ;;
   ;; Illustrate 12 + Log O/H vs mass bias
   ;;


   rr = _get_one(/all)
   rr2 = _get_one(/sn4363)
   rrssfr = _get_one(/all, /ssfr)
   rr2ssfr = _get_one(/sn4363, /ssfr)





   jwst_1st_dirs, dd

   
   if keyword_set(ps) then begin
      ps_on, dd.figdir+'mass-met-bias.ps', /color, /times, aspect=0.6, xlen=12
      !p.font = 0
      !x.thick = 2
      !y.thick = 2
      cz = 1.6
      tk = 4
   endif else begin
      cz = 1.4

   endelse
   
   fill_colors = ['gray', 'orange', 'cyan']
   colors = ['black', 'FireBrick', 'blue']
   

   plot_summary_shaded, rr, /p68,  xrange=[7.8, 11], ytitle='12 + Log O/H', $
                        xtitle=TeXtoIDL('Log M_*'), yrange=[7.5, 9.3], /ys, $
                        charsize=cz, /xs
   plot_summary_shaded, rr, /p68, /overplot, shade_colors=[fill_colors[0]], /dofill, $
                        /conf_med
   
   ok = where(rr2.median gt 7.5)
   plot_summary_shaded, rr2ssfr, /p68, /overplot, color=jb_colour(colors[1]), $
                        /dofill, shade_colors=[fill_colors[1]], ok=ok, /conf_med
   
   plot_summary_shaded, rrssfr, /p68, /overplot, color=jb_colour(colors[2]), $
                        /dofill, shade_colors=[fill_colors[2]], /conf_med


   xpos = 10.1
   ypos = 7.9
   dy = 0.1
   dx = 0.2
   h = 0.04
   ddx = 0.03
   ddy = -0.01
   
   labels = ['All SF', 'S/N [O III]4363>7', 'All SF SFR/M_* > 10^{-9}']
   

   for i=0L, n_elements(labels)-1 do begin

      ythis = ypos-i*dy
      ybottom = ythis-0.5*h+[0,0]
      ytop = ythis+0.5*h+[0,0]
      xx = xpos+[0, dx]
      get_polyfill_arr, xx, ybottom, ytop, zx, zy, xrange=!x.crange, yrange=!y.crange

      polyfill, zx, zy, color=jb_colour(fill_colors[i])
      oplot, xx, ythis+[0, 0], color=jb_colour(colors[i]), thick=tk

      xyouts, xpos+dx+ddx, ythis+ddy, TeXtoIDL(labels[i]), charsize=cz
      
   endfor
   

   jwst_1st_analogues_massmet, /overplot
   
   if keyword_set(ps) then begin
      ps_off
      my_cleanplot, /silent
   endif
   

   ;;
   ;; SFR plots.
   ;;
   if keyword_set(ps) then begin
      ps_on, dd.figdir+'mass-met-bias-vs-SFR.ps', /color, /times, aspect=0.8, xlen=10
      !p.font = 0
      !x.thick = 2
      !y.thick = 2
      cz = 1.6
      tk = 4

   endif
   SFR_bins = [[-2, -1.0], [-1.0, 0], [0, 1]]
   dims = size(SFR_bins, /dimen)
   n_bins = dims[1]
   plot, [7.8, 10.3], [-0.5, 1], /nodata, xtitle=TeXtoIDL('Log M_*'), $
         /ys, ytitle=TeXtoIDL('\Delta Log O/H [all - with 4363]'), $
         charsize=cz, /xs
   colors = ['#F71467', '#1E88E5', '#FFC107', '#009E73', '#88CCEE', '#004D40']
   symbols, 2, 1

   ypos = 0.0
   xpos = 9.6
   dy = -0.05
   ddy = -0.01
   dx = 0.06
   sym = [30, 31, 32, 20, 2]
   symsize = [0.6, 0.7, 0.8, 0.6, 2]

   for i=0L, n_bins-1 do begin
      rr_sfr = _get_one(/all, sfr = SFR_bins[*, i])
      rr2_sfr = _get_one(/sn4363, sfr = SFR_bins[*, i])
      diff = rr_sfr.median-rr2_sfr.median
      ddiff = sqrt((0.5*(rr_sfr.med_p84-rr_sfr.med_p16))^2+ $
                   (0.5*(rr2_sfr.med_p84-rr2_sfr.med_p16))^2)

      ok = where(rr2_sfr.median gt 7.5)
      symbols, sym[i], symsize[i]
      myoploterr2, rr_sfr.x[ok], diff[ok], rr_sfr.x[ok], rr_sfr.x[ok], $
                   diff[ok]-ddiff[ok], diff[ok]+ddiff[ok], $
                   psym=8, color=jb_colour(colors[i]), errcol=jb_colour(colors[i])


      label = string(format='("[",I2,", ",I2,"]")', $
                     SFR_bins[*, i])
      plots, xpos, ypos+i*dy, psym=8, color=jb_colour(colors[i])
      xyouts, xpos+dx, ypos+i*dy+ddy, label, charsize=cz
      
   endfor
   xyouts, xpos+dx*0.9, ypos-dy, TeXtoIDL("Log SFR [M_{sun}/yr]"), charsize=cz
      


   if keyword_set(ps) then begin
      ps_off
      my_cleanplot, /silent
   endif
   
   
   
end
