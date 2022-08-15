pro jwst_1st_show_neiv_6355, ps=ps
   ;;
   ;; Show this line.
   ;;

   d = jwst_1st_get_coadded(6355, d1=d1, d2=d2)
   sp = jwst_nirspec_reextract(d, 19, narrow=3, /ap)
   sp1 = jwst_nirspec_reextract(d1['combined'], 19, narrow=3, /ap)
   sp2 = jwst_nirspec_reextract(d2['combined'], 19, narrow=3, /ap)

   if keyword_set(ps) then begin
      jwst_1st_dirs, dd
      ps_on, dd.figdir+'neiv_in_6355.ps', aspect=0.6, xlen=12, $
             /color, /times
      !p.font = 0
      !x.thick = 2
      !y.thick = 2
   endif

   redshift = 7.664
;   get_ct, 1, /reverse;, white_at=0, black_at=1
   loadct, 1
   jwst_vis_spec_1d2d, sp, xrange=[2, 2.2], /zscale,   $
                       range=[3e-7, 3e-6], smooth=1.5, $
                        xtitle='Wavelength [micron]', yrange=[-5e-7, 2e-6]
   ;;                     
   ;; jwst_vis_spec_1d2d, sp, xrange=[2, 2.2], /zscale,   $
   ;;                     range=[0, 4], smooth=1.5, $
   ;;                     xtitle='Wavelength [micron]', yrange=[-5e-7, 2e-6], $
   ;;                     snimage=1

;   stop
   oplot, sp1.oned.wavelength, sp1.oned.flux, color=jb_colour('#F71467')
   oplot, sp2.oned.wavelength, sp2.oned.flux, color=jb_colour('#1E88E5')

   l_neiv = (1+redshift)*[2422.55, 2425.16]/1e4 ; micron
   for i=0L, 1 do $
    oplot, l_neiv[i]+[0, 0], [1.5, 2]*1e-6, linestyle=2
   

   
   if keyword_set(ps) then begin
      ps_off
      my_cleanplot, /silent
   endif
   
end
