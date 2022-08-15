pro _check_agn, inds, class
   ;;
   ;; Summarise the ionization conditions.
   ;;
   
   i_class = [1, 3, 4]
   labels = ['Star-forming', 'Composite', 'AGN']

   for i=0L, n_elements(i_class)-1 do begin
      ii = where(class.i_class[inds] eq i_class[i], n_this)

      print, format='("   N(",A,") = ",I4,"   (",F5.2,"%)")', labels[i], n_this, $
             100.0*n_this/float(n_elements(inds))

   endfor
 

end
   

pro _find_all, s, class, lgm,  mugas, oh, dist=dist, ps=ps

   todo = ['04590', '05144', '06355', '08140', '10612']
   colors = ['#F71467', '#1E88E5', '#FFC107', '#009E73', '#004D40']


   if n_elements(mugas) eq 0 then $
    mugas = dr_load_x('fib-mugas', release='dr7')
   if n_elements(oh) eq 0 then $
    oh = dr_load_x('fib-oh', release='dr7', /steep)
   
   jwst_1st_dirs, dd
   if keyword_set(ps) then begin
      ps_on, dd.figdir+'bpt-with-analogues.ps', aspect=0.6, xlen=12, /times, /color
      !p.font = 0
      !x.thick = 2
      !y.thick = 2
      cz = 1.8
   endif else cz = 1.8
   
;   colors = ['red', 'blue', 'yellow', 'ForestGreen', 'pink']
   show_bpt_contours, yrange=[0, 1.5], $
                      charsize=cz, xrange=[-2.5, 0.5]
   ypos = 1.42
   dy = 0.05
   ddy = -0.005
   xpos = -1.8
   dx = 0.05
   
   for i=0L, n_elements(todo)-1 do begin
      jwst_1st_find_analogues, todo[i], s, class, lgm, mugas, oh, dist=dist, $
                               color=colors[i], /plot
      
      symbols, 2, 1
      plots, xpos, ypos-dy*i, psym=8, color=jb_colour(colors[i]), symsize=2
      symbols, 1, 1
      plots, xpos-0.03, ypos-dy*i, psym=8, color=jb_colour(colors[i]), symsize=2
      xyouts, xpos+dx, ypos-dy*i+ddy, string(format='(I5)', long(todo[i])), $
              charsize=cz
      
   endfor

   if keyword_set(ps) then begin
      ps_off
      my_cleanplot, /silent
   endif

   
end

pro jwst_1st_find_analogues, object, s, class, lgm, mugas, oh, dist=dist, $
                             color=color, plot=plot
   ;;
   ;; Find analogues to the JWST objects. I will define them in two
   ;; ways: 4363/Hg analogues and 3869/[O II] analogues. 
   ;;
   ;; The distance is in dex...
   if (n_elements(dist) eq 0) then dist = 0.1
   
   ;;           0           1            2          3            4         5
   lines = ['OII3727', 'NEIII_3869', 'H_GAMMA', 'OIII_4363', 'H_BETA', 'OIII_5007']
   
   if object eq '05144' then touse = 'full'
   jwst_1st_line_luminosities, object, lines, $
                               lum=lum, dlum=dlum, /silent, $
                               touse=touse

   dlum = 2.0*dlum

   r_o3hb = alog10(s.oiii5007/s.hb)
   r_o3hb_jwst = alog10(lum[5]/lum[4])
   ;;; I ignore uncertainties on the [O III]/Hb ratio 
   
   flag = intarr(n_elements(r_o3hb))
   ii = where(abs(r_o3hb-r_o3hb_jwst) lt dist)
   flag[ii] = 1
;   stop

   if object eq '08140' then begin
      n_ok_4363hg = 0
      i_4363hg = []
   endif else begin
      ;; 4363/Hg friends
      r_sdss = alog10(s.oiii4363/s.hg)
      r_jwst = lum[3]/lum[2]
      dr_jwst = r_jwst*error_on_fraction(lum[3], dlum[3], lum[2], dlum[2])
      r_jwst_low = alog10(r_jwst-dr_jwst)
      r_jwst_hi = alog10(r_jwst+dr_jwst)
      r_jwst = alog10(r_jwst)
      
      i_4363hg = where(s.oiii4363/s.doiii4363 gt 7 and s.hg/s.dhg gt 7 and $
                       r_sdss gt r_jwst_low $
                       and r_sdss lt r_jwst_hi and s.zxcor gt 0.01 and flag eq 1, n_ok_4363hg)
      
      if (n_ok_4363hg eq 0) then stop
   endelse 

   ;; [Ne III]3869/[O II]3727 friends
   r_sdss = alog10(s.neiii3869/s.oii3727)
   
   r_jwst = lum[1]/lum[0]
   dr_jwst = r_jwst*error_on_fraction(lum[1], dlum[1], lum[0], dlum[0])
   r_jwst_low = alog10(r_jwst-dr_jwst)
   r_jwst_hi = alog10(r_jwst+dr_jwst)
   r_jwst = alog10(r_jwst)

   i_3869o2 = where(s.oii3727/s.doii3727 gt 7 and s.neiii3869/s.dneiii3869 gt 7 and $
                    r_sdss gt r_jwst_low and r_sdss lt r_jwst_hi and s.zxcor gt 0.01 $
                    and flag eq 1, n_ok_3869o2)

   ;; [Ne III]3869/Hg friends
   r_sdss = alog10(s.neiii3869/s.hg)
   r_jwst = alog10(lum[1]/lum[2])
   i_3869hg = where(s.hg/s.dhg gt 7 and s.neiii3869/s.dneiii3869 gt 7 and $
                    abs(r_sdss-r_jwst) lt dist and s.zxcor gt 0.01 and flag eq 1, n_ok_3869hg)


   
   print, format='("For ",A," I found ",I0," 4363/Hg neighbours,  ",I0," for 3869/3727, and ",I0," for 3869/Hg")', $
          object, n_ok_4363hg, n_ok_3869o2, n_ok_3869hg

   if n_ok_4363hg gt 0 then begin
      print, '---------------- 4363/Hg -------------'
      _check_agn, i_4363hg, class
   endif 
   if n_ok_3869o2 gt 0 then begin
      print, '---------------- 3869/3727 -------------'
      _check_agn, i_3869o2, class
   endif 
   if n_ok_3869hg gt 0 then begin
      print, '---------------- 3869/Hg -------------'
      _check_agn, i_3869hg, class
   endif 

   if keyword_set(plot) then begin
;      allinds = [i_4363hg, i_3869o2, i_3869hg]
;      allinds = [i_4363hg, i_3869o2, i_3869hg]
      x_bpt = alog10(s.nii6584/s.ha)
      y_bpt = alog10(s.oiii5007/s.hb)
      symbols, 2, 1
      if (n_ok_4363hg gt 0) then $
       oplot, x_bpt[i_4363hg], y_bpt[i_4363hg], psym=8, color=jb_colour(color)
      symbols, 1, 1
      plots, x_bpt[i_3869o2], y_bpt[i_3869o2], psym=8, color=jb_colour(color)
      
   endif

   ;;
   ;; Next, we want the M*-O/H plane
   ;;
   jwst_1st_dirs, dd
   outfile = dd.root+'CL01Fit/MOH-analogues-SF-'+string(format='(I5.5)', object)+'h5'

   if (n_ok_4363hg gt 0) then begin
      tmp = bytarr(n_elements(lgm.lgm))
      tmp[i_4363hg] = 1
      use = where(tmp eq 1 and class.i_class eq 1, n_use)
      if (n_use eq 0) then stop
      
      struct_to_hdf5, {lgm: lgm.lgm[use], oh: oh.median[use]}, outfile
   endif   
   

   ;; FInally, assemble the Mugas distributions.
   CL01FITS = '/data2/jarle/SDSS/CL01FitsPerObj/'
   dim = hdf5_to_struct(CL01FITS+'dim.h5')
   first = 1
   for i=0L, n_ok_4363hg-1 do begin
      idr7 = i_4363hg[i]
      if class.i_class[idr7] ne 1 then continue ; Skip all non-SF
      pid = s.plateid[idr7]
      mjd = s.mjd[idr7]
      fid = s.fibre[idr7]

      subdir = CL01FITS+string(format='(I4.4)', pid)+'/'
      fname = subdir+string(format='("cl01fit-",I4.4,"-",I5,"-",I3.3,".h5")', $
                            pid, mjd, fid)
      h = hdf5_to_struct(fname)
      
      if first then begin
         xMUGAS = h.xmugas
         yMUGAS = h.bin_mugas*0.0d0

         xSFR = h.xlogsfr
         ySFR = h.bin_logsfr

         xoh = h.xoh
         yOH = h.bin_oh

         xTE = h.xte
         yTE = h.bin_te

         ylogU = h.bin_U
         first = 0
      endif

      yMUGAS = yMUGAS + h.bin_mugas/total(h.bin_mugas)
      ySFR = ySFR + h.bin_logsfr/total(h.bin_logsfr)
      yOH = yOH + h.bin_oh/total(h.bin_oh)
      yTE = yTE + h.bin_te/total(h.bin_te)
      ylogU = ylogU + h.bin_U/total(h.bin_U)
      
   endfor
   

   jwst_1st_dirs, dd
   if (n_elements(xmugas) gt 0 ) then begin
      outfile = dd.root+'CL01Fit/cl01-analogues-'+string(format='(I5.5)', object)+'.h5'

      ss = {x_mugas: xMUGAS, pdf_mugas: yMUGAS, x_logsfr: xSFR, pdf_logsfr: ySFR, $
            x_oh: xOH, pdf_oh: yOH, x_TE: xTE, pdf_TE: yTE, pdf_logU: ylogU, $
            x_logU: dim.logu}
      struct_to_hdf5,ss , outfile
   endif
   
end



pro jwst_1st_analogues_massmet, overplot=overplot
   ;;
   ;; Show the mass metallicity diagram for the analogues.
   ;;

   jwst_1st_dirs, dd

   ;; We can only compare to 4590, 6355 and 10612 as those have mass
   ;; estimates. 
   todo = ['04590', '06355', '10612']
   colors = ['#F71467', '#FFC107', '#004D40']
   oh_cl01 = [7.619, 7.946, 7.789]
   mstar=[7.9, 8.7, 8.1]
   symbols = [30, 31, 32]
   symsize = [1.0, 1.0, 1.3]
   
   ;; This function is in jwst_1st_mass_met_illustration.pro 
   rr = _get_one(/all)


   if not keyword_set(overplot) then begin
      plot_summary_shaded, rr, /p68, yrange=[7.5, 9.2], /ys, xrange=[7.8, 10.5], /xs
   endif
   
   for i=0L, n_elements(mstar)-1 do begin
      outfile = dd.root+'CL01Fit/MOH-analogues-SF-'+todo[i]+'h5'
      h = hdf5_to_struct(outfile)

      symbols, symbols[i], symsize[i]
      if i eq 1 then scl = 0.3 else scl = 0.7
      oplot, h.lgm, h.oh, psym=8, color=jb_colour(colors[i]), symsize=scl
      
;      symbols, 30, 1
      plots, mstar[i], oh_cl01[i], psym=8, color=jb_colour(colors[i]), symsize=1.5

;      mveq, h.oh
      
   endfor
   


end
