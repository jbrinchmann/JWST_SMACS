pro jwst_1st_illustrate_flux_normalisation, ps=ps

   to_show = 3042
   propid = 2736

   
   ;; The best results are from the default NIRCam photometry 
                                ;infix = '-idl-'
   infix = ''

   ROOT = '/data2/jarle/JWST/ERO/SMACS/'
   aper = mrdfits(ROOT+'Stamps/AperturePhotometry/photometry_idl.fits', 1)
   nc_mag = hdf5_to_struct(ROOT+'NIRCam-photometry-NIRSpec.h5')
   
   tmpid = string(format='(I0,"_",I0)', propid, to_show)
;   i_obj = where(strtrim(aper.source) eq tmpid, N)
   i_obj = where(strtrim(nc_mag.sources) eq  tmpid, N)
   if (N eq 0) then stop
   
;   stop
;   f = [aper[i_obj].f200w, aper[i_obj].f277w, aper[i_obj].f356w, aper[i_obj].f444w]*8.3359649e-9
   f = [nc_mag.aper50_flux_f200W[i_obj], nc_mag.aper50_flux_f277w[i_obj], $
        nc_mag.aper50_flux_f356w[i_obj], nc_mag.aper50_flux_f444w[i_obj]]
      
   outfile = ROOT+'Spectra/renorm-'+infix+string(format='(I5.5,"_",I5.5)', propid, to_show)+'.fits'
   sp_renorm = mrdfits(outfile, 1)

   ;; We also want the default spectrum
   d = jwst_1st_get_coadded(to_show)
   ;; Unnormalised conversion to Fnu
   fnu_orig = d.oned.flux*d.oned.wavelength^2
   
   lpivot = [1.989, 2.762, 3.568, 4.408]
   i_pivot = 1
   dum = min(abs(d.oned.wavelength-lpivot[i_pivot]), i_norm)
   l_norm = d.oned.wavelength[i_norm]
   
   ii_norm = where(abs(d.oned.wavelength-l_norm) lt 150/1e4 and $
                   finite(d.oned.flux) eq 1 and d.oned.flux ne 0.0, n_norm)
   fnu_orig_norm = median(fnu_orig[ii_norm])
   orig_scale = f[i_pivot]/fnu_orig_norm
   fnu_orig_scaled = fnu_orig*orig_scale


   if keyword_set(ps) then begin
      ps_on, ROOT+'Figures/illustration_of_flux_calibration.ps', $
             /color, /times, aspect=0.4, xlen=12
      !p.font = 0
      !x.thick = 2
      !y.thick = 2
      !p.thick = 1.8
      cz = 1.6
   endif else cz = 1.5

   old_col = 'gray'
   plot, d.oned.wavelength, 1e6*fnu_orig_scaled, ytitle=TeXtoIDL('F_{\nu} [\mu')+'Jy]', /nodata, $
         xtitle=TeXtoIDL('Wavelength [\mu'+'m]'), /xs, yrange=[0, 2], $
         charsize=cz
   oplot, d.oned.wavelength, 1e6*fnu_orig_scaled, color=jb_colour(old_col), $
          linestyle=1
   oplot, d.oned.wavelength, 1e6*sp_renorm.flux_nu


   
   symbols, 30, 1
   oplot, lpivot, 1e6*f, psym=8, symsize=1, color=jb_colour('red')
   jb_text, 0.5, -0.1,  string(to_show), charsize=1.8, /frac, /relative, align=0.5


   dy = -0.1
   oplot, [1.9, 2.1], 1.85+[0, 0]+dy, color=jb_colour(old_col), linestyle=1
   xyouts, 2.13, 1.84+dy, 'Original', charsize=cz

   oplot, [1.9, 2.1], 1.70+[0, 0]+dy
   xyouts, 2.13, 1.69+dy, 'Renormalised', charsize=cz
   
   
   if keyword_set(ps) then begin
      ps_off
      my_cleanplot, /silent
   endif

   
   
end
