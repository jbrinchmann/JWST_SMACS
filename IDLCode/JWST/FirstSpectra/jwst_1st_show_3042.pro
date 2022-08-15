pro jwst_1st_show_3042, ps=ps
   ;;
   ;; Show the spectrum, with line identifications and the stamps. 
   ;;

   jwst_1st_dirs, dd

   common linelist, lname, l_line, type

   
   infix = '-idl-'
   outfile = dd.ROOT+'Spectra/renorm-'+infix+string(format='(I5.5,"_",I5.5)', 2736, 3042)+'.fits'

   ;; Load spectrum and define plottable quantities
   sp = mrdfits(outfile, 1)

   
   z = 1.99381
   x = sp.wavelength
   y = alog10(sp.flux)
   

;   stop
   lname = TeXtoIDL(['H\alpha', '[N II]6584', '[S II]6717', '[S II]6731', $
                      '[S III]9068', '[S III]9533', 'Pa\delta', $
                     'Pa\gamma', 'Pa\beta', 'He I 1.083', '[Fe II]1.257'])


   l_line = [6564.61, 6585.28,  6718.29,  6732.67, 9068.60, $
             9533.2, 10052.13, 10941.09, 12821.59, 10833.22, 12570.24]
   type[*] = 'em'

;   stop
   
   xmin = 0.1
   xmax = 0.9

   ymin = 0.1
   dy_stamp = 0.2
   ygap = 0.1
   ymin_spec = ymin+dy_stamp+ygap
   ymax = 0.85

   ;; The spectrum position - the stamp positions are calculated next
   sp_pos = [xmin, ymin_spec, xmax, ymax]

   if keyword_set(ps) then begin
      ps_on, dd.figdir+'plot-3042a.ps', aspect=0.4, xlen=12, /times, $
             /color
      !p.font = 0
      !x.thick = 2
      !y.thick = 2
   endif
   my_lineid_plot, x, y, l_line/1e4, lname, /extend, elinestyle=1, $
                   lcolor=jb_colour('orange'), $
                   xtitle=TeXtoIDL('Observed wavelength [\mu')+'m]', $
                   yrange=[-19.8, -18.5], position=sp_pos, /xs, $
                   ytitle='Log F', /ys
   
   
   filters = ['f435w', 'f606w', 'f814w', 'f090w', 'f150w', $
              'f200w', 'f277w', 'f356w', 'f444w']

   pos = get_position_arr(0, nx=n_elements(filters), ny=1, $
                          xmin=xmin, xmax=xmax, ymin=ymin, $
                          ymax=ymin+dy_stamp, xgap=0.02)

   for i=0L, n_elements(filters)-1 do begin
      s = jwst_1st_load_stamp(3042, filters[i])
      pos = get_position_arr(i)
      
      get_ct, 0, /reverse
      tvim_true, s.im, /noframe, position=pos, /noerase
      jb_text, 0.5, -0.1, filters[i], /fraction, /relative, $
               color=jb_colour('red'), align=0.5
   endfor


   if keyword_set(ps) then begin
      ps_off
   endif else stop
   
   d = jwst_1st_get_coadded(3042)
   fsp = mrdfits(dd.specdir+'Platefit/renorm-fit-spec-02736_03042-fixedlinelist.fits', 1)

   d.oned.flux = fsp.flux-(fsp.continuum+fsp.resid_cont)

   reverse = 1
   ctab = 3
   lcolor = 'orange'
   ;; I will glue a bunch of things together in Illustrator
   if keyword_set(ps) then begin
      ps_on, dd.figdir+'plot-3042b.ps', aspect=0.9, xlen=12, /times, $
             /color
   endif
   get_ct, ctab, reverse=reverse
   jwst_vis_spec_1d2d, d,  charsize=1.6, pcharsize=1.6, lcharsize=1.4, $
                       /zscale, redshift=1.99381, /label, xrange=[1.9, 2.05], $
                       range=[0, 8e-6], ytickformat='noticks', lcolor=lcolor
   if keyword_set(ps) then begin
      ps_off
   endif else stop


   ;; I will glue a bunch of things together in Illustrator
   if keyword_set(ps) then begin
      ps_on, dd.figdir+'plot-3042b-alt.ps', aspect=0.6, xlen=12, /times, $
             /color
   endif
   get_ct, ctab, reverse=reverse
   jwst_vis_spec_1d2d, d,  charsize=1.6, pcharsize=1.6, lcharsize=1.4, $
                       /zscale, redshift=1.99381, /label, xrange=[1.7, 3.9], $
                       range=[0, 6e-6], ytickformat='noticks', lcolor=lcolor
   if keyword_set(ps) then begin
      ps_off
   endif else stop

   
   if keyword_set(ps) then begin
      ps_on, dd.figdir+'plot-3042c.ps', aspect=0.6, xlen=12, /times, $
             /color
   endif
   get_ct, ctab, reverse=reverse
   jwst_vis_spec_1d2d, d,  charsize=1.6, pcharsize=1.6, lcharsize=1.4, $
                       /zscale, redshift=1.99381, /label, xrange=[2.6, 3.1], $
                       range=[0, 3.5e-6], ytickformat='noticks', lcolor=lcolor
   if keyword_set(ps) then begin
      ps_off
   endif else stop

;   stop
   if keyword_set(ps) then begin
      ps_on, dd.figdir+'plot-3042d.ps', aspect=0.6, xlen=12, /times, $
             /color
   endif
   get_ct, ctab, reverse=reverse
   jwst_vis_spec_1d2d, d,  charsize=1.6, pcharsize=1.6, lcharsize=1.4, /zscale, $
                       redshift=1.99381, /label, xrange=[3.2, 4], range=[0, 3e-6], $
                       ytickformat='noticks', lcolor=lcolor
   if keyword_set(ps) then begin
      ps_off
   endif else stop
   
   if keyword_set(ps) then begin
      ps_on, dd.figdir+'plot-3042e.ps', aspect=0.9, xlen=12, /times, $
             /color
   endif
   get_ct, ctab, reverse=reverse
   jwst_vis_spec_1d2d, d,  charsize=1.6, pcharsize=1.6, lcharsize=1.4, /zscale, $
                       redshift=6.253, /label, xrange=[4.4, 5.1], range=[0, 2e-6], $
                       ytickformat='noticks', lcolor=lcolor
   if keyword_set(ps) then begin
      ps_off
   endif else stop
   

end
