pro jwst_1st_show_te_values, s, class, ab=ab, ps=ps
   ;;
   ;; Show that the Te values for the JWST galaxies relative to CL01
   ;; temperatures are reasonable.
   ;;
   if (n_elements(ab) eq 0) then $
    ab = mrdfits('/Users/jarle/Work/Wolf-Rayet/DR7/Abundances/wr_abundances_dr7_table.fits', 1)

   ss = hdf5_to_struct('/data2/jarle/SDSS/Datafiles/DR7/fib-te_dr7_v5_2.h5')

   x_jwst = [2.52093, 1.4455, 1.39225, 1.53705]
   y_jwst = [1.6693357, 1.1893408, 1.4390195, 1.5593959]

   ;; The subsample 'ss' was calculated for.
   use_ss = where(class.i_class eq 1)
   sn_4363 = where(s.oiii4363[use_ss]/s.doiii4363[use_ss] gt 7)

   te_o3 = ab[use_ss].te_o3


   if keyword_set(ps) then begin
      jwst_1st_dirs, dd
      ps_on, dd.figdir+'te_comparison_sdss_jwst.ps', aspect=0.6, $
             /color, /times, xlen=12
      !p.font = 0
      !x.thick = 2
      !y.thick = 2
      cz = 1.6
   endif else cz = 1.8
   plot, [0, 1], [0, 1], /nodata, psym=8, $
         xtitle=TeXtoIDL('T_e [10^4 K, direct method]'), $
         ytitle=TeXtoIDL('T_e [10^4 K, CL01fit]'), charsize=cz, $
         xrange=[0.8, 2.6], /xs, yrange=[0.5, 2.0]
   
   oplot, te_o3[sn_4363], ss.median[sn_4363]/1e4, psym=8, $
          color=jb_colour('gray50')

   plots, x_jwst, y_jwst, psym=8, symsize=3, color=jb_colour('red')

   abline, 0, 1, color=jb_colour('orange'), linestyle=2, thick=tk

   if keyword_set(ps) then begin
      ps_off
      my_cleanplot, /silent
   endif

   
end
