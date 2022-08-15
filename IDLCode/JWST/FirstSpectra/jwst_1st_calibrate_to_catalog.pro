pro jwst_1st_calibrate_to_catalog
   ;;
   ;; I was unable to exactly reproduce the photometry in the
   ;; catalogue (although I get close) so I am doing a double-check
   ;; here. 
   ;;


   filters = ['f090w', 'f150w', 'f200w', 'f277w', 'f356w', 'f444w']

   ROOT = '/data2/jarle/JWST/ERO/SMACS/'
   nc = hdf5_to_struct(ROOT+'NIRCam-photometry-NIRSpec.h5')
   PHOTDIR = ROOT+"Stamps/SExtractor/"
   idlphot = mrdfits(ROOT+"Stamps/AperturePhotometry/photometry_idl.fits", 1)
   
   
   scale = dblarr(n_elements(filters))
   scale_idl = dblarr(n_elements(filters))
   pos = get_position_arr(0, nx=3, ny=2, xgap=0.05, ygap=0.05)
   for i_f=0L, n_elements(filters)-1 do begin
      
      phot_file = PHOTDIR+'phot-'+filters[i_f]+'-aligned.fits'
      t_my = mrdfits(phot_file, 1)

      ;; We compare half-light fluxes from the catalogue to the BEST
      ;; fluxes from SExtractor - this because the total fluxes seem
      ;; very poor in the catalogue, and then I multiply this by 2 to
      ;; get total light
      f_nircam = struct_var(nc, 'aper50_flux_'+filters[i_f])*2
      f_my = t_my.flux_best
      f_idl = struct_var(idlphot, filters[i_f])

      df_nircam = struct_var(nc, 'aper50_flux_err_'+filters[i_f])*2
      df_my = t_my.fluxerr_best
      df_idl = struct_var(idlphot, 'd'+filters[i_f])

      ratio = f_nircam/f_my
      ratio_idl = f_nircam/f_idl
      ok = where(f_my gt 0 and f_nircam gt 0 and finite(ratio) eq 1 and $
                 ratio gt 0 and ratio lt 1e-7, n_ok)
      ok_idl = where(f_idl gt 0 and f_nircam gt 0 and finite(ratio_idl) eq 1 and $
                 ratio_idl gt 0 and ratio_idl lt 1e-7, n_ok)

      scale[i_f] = mean(ratio[ok])
      scale_idl[i_f] = mean(ratio_idl[ok_idl])
      pos = get_position_arr(i_f)
      plot, ok, ratio[ok], psym=8, title=filters[i_f], position=pos, $
            noerase=(i_f ne 0), /ylog
      oplot, ok_idl, ratio_idl[ok_idl], psym=8, color=jb_colour('CornFlowerBlue')
      oplot, !x.crange, scale[i_f]+[0,0], color=jb_colour('red'), linestyle=2
      oplot, !x.crange, scale_idl[i_f]+[0,0], color=jb_colour('orange'), linestyle=2
      

      if i_f eq 2 then stop
      
   endfor

   ;; I can almost but not quite reproduce this scaling by assuming
   ;; that I should multiply by 10^6 [from MJy], pixelscale^2 and
   ;; divide by (180/pi)^2 to get pixel scale in arcsec^2 which I
   ;; think is what SEXtractor assumes. This seems right because when
   ;; I use the -aligned versus -orig I see a difference for the
   ;; reddest filters which have different pixelscales.

   ;; Oh well, I'll stick with this.
   
   print, 'For the aligned stamps the scale is ', mean(scale)
   print, 'and versus IDL photometry ', mean(scale_idl)
   stop   
end
