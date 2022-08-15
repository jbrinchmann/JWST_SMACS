function jwst_1st_load_one_gaussfit, object, pipeline_error=pipeline_error

   jwst_1st_dirs, dd
   GDIR = dd.specdir+'GaussianFits/'

   fname = GDIR+'gaussfit-'+object
   if keyword_set(pipeline_error) then $
    fname = fname+'-pipe_err'
   fname = fname+'.txt'

   readcol, fname, skip=1, format='A,F,F,F,F,F,F,F,F', $
            line, flux, dflux, lcenter, dlcenter, sigma, dsigma, level ,dlevel, $
            /silent

   t = {line: line, flux: flux, dflux: dflux, lcenter: lcenter, $
        dlcenter: dlcenter, sigma: sigma, dsigma: dsigma, level: level, dlevel: dlevel}
   return,t 


end
