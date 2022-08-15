function _get_flux, l, f, lpivot, width=width

   width = 150.0/1e4

   dum = min(abs(l-lpivot), i_norm)

   l_norm = l[i_norm]
   
   ii_norm = where(abs(l-l_norm) lt width and $
                   finite(f) eq 1 and f ne 0.0, n_norm)

   return, median(f[ii_norm])

   
end


pro jwst_1st_simple_flux_scaling
   ;;
   ;; Try a very simple flux scaling:
   ;;
   ;; Look at the spectrum near the pivot wavelength and since there
   ;; are two filters per grism, try find a solution for
   ;;
   ;;  f_obs*a + b = ftrue at each wavelength
   ;;
   ;;


      
   ROOT = '/data2/jarle/JWST/ERO/SMACS/'
   SPECDIR = ROOT+'Spectra/'

   t = mrdfits(ROOT+'match-ids.fits', 1)
   zcat = mrdfits(ROOT+'redshift-for-sources.fits', 1)

   ;; The original catalogues
   nc_mag = hdf5_to_struct(ROOT+'NIRCam-photometry-NIRSpec.h5')

   ;; My IDL aper photometry
   aper = mrdfits(ROOT+'Stamps/AperturePhotometry/photometry_idl.fits', 1)
   aper = arr_of_struct_to_struct_of_arr(aper)

   lpivot = [1.989, 2.762, 3.568, 4.408]

   nobj = n_elements(t)
   
   for i=0L, nobj-1 do begin

      dum = strsplit(t[i].specid, '_', /extract)
      propid = long(dum[0])
      objid = long(dum[1])

      if objid ne 5144 then continue

      ;; Get fluxes from NIRCAM

      if keyword_set(use_idlphot) then begin
         ;; Get the fluxes and scale to the provided catalog fluxes.
         f = [aper.f200w[i], aper.f277w[i], aper.f356w[i], aper.f444w[i]]*8.3359649e-9 
      endif else begin
         f = [nc_mag.aper50_flux_f200W[i], nc_mag.aper50_flux_f277w[i], $
              nc_mag.aper50_flux_f356w[i], nc_mag.aper50_flux_f444w[i]]
         
      endelse

      
      ;; Load the spectrum and figure out the blue and red grisms
      ;;
      d = jwst_1st_get_coadded(objid, blank=0)
      diffR=d.oned.R-shift(d.oned.R, 1)
      diffR[0] = 0
      big_jump = where(abs(diffR) gt 100, n_jump)
      if (n_jump eq 0) or (n_jump gt 1) then stop
      l_jump = d.oned.wavelength[big_jump[0]]
      blue = where(d.oned.wavelength lt l_jump, complement=red)
      
      ;; We need fNu for this normalisation - since I will renormalise
      ;; below I am leaving out the normalisation factor (units &
      ;; speed of light)
      fnu_orig = d.oned.flux*d.oned.wavelength^2*1e-8


      ;;----------------------------
      ;; Find scale for blue
      ;;----------------------------
      width = 150
      spec_f0 = _get_flux(d.oned.wavelength, fnu_orig, lpivot[0], width = width/1e4)
      spec_f1 = _get_flux(d.oned.wavelength, fnu_orig, lpivot[1], width = width/1e4)

      aB = (f[0]-f[1])/(spec_f0-spec_f1)
      bB = f[1]-aB*spec_f1


      ;;----------------------------
      ;; Find scale for red
      ;;----------------------------
      spec_f2 = _get_flux(d.oned.wavelength, fnu_orig, lpivot[2], width = width/1e4)
      spec_f3 = _get_flux(d.oned.wavelength, fnu_orig, lpivot[3], width = width/1e4)

      aR = (f[2]-f[3])/(spec_f2-spec_f3)
      bR = f[3]-aR*spec_f3

      ;; Caclulate scaled spectrum
      fnu_new = fnu_orig*0
      fnu_new[blue] = fnu_orig[blue]*aB+bB
      fnu_new[red] = fnu_orig[red]*aR+bR

;      stop
      

   endfor
end
