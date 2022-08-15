pro jwst_1st_get_magnitudes
   ROOT = '/data2/jarle/JWST/ERO/SMACS/'
   SPECDIR = ROOT+'Spectra/'

   t = mrdfits(ROOT+'redshift-for-sources.fits', 1)
   filters = ['F200W', 'F277W', 'F356W', 'F444W']

   mag = dblarr(n_elements(t), n_elements(filters))
   dmag = mag
   flux = mag
   dflux = mag
   
   for i=0L, n_elements(t)-1 do begin
      print, "Doing "+t[i].object
      tmp = strsplit(t[i].object, '_', /extract)
      propid = string(format='(I5.5)', long(tmp[0]))
      objid = string(format='(I5.5)', long(tmp[1]))

      d = jwst_1st_get_coadded(long(objid), blank=0)

      res = jwst_nirspec_magnitude(d.oned, filters)

      mag[i, *] = res.mag
      dmag[i, *] = res.dmag
      flux[i, *] = res.flux
      dflux[i, *] = res.dflux

   endfor

   struct_to_hdf5, {object: t.object, mag: mag, dmag: dmag, $
                    flux: flux, dflux: dflux, filters: filters}, $
                   ROOT+'raw-specmags.h5'
      
end
