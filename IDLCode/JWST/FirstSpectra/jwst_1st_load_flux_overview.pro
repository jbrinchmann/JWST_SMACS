function jwst_1st_load_flux_overview, specid
   
   jwst_1st_dirs, dd
   DIR = dd.specdir+'Platefit/FluxOverviews/'

   
   if size(specid, /tname) ne 'STRING' then $
    id = string(format='(I5.5)', specid) $
   else id = specid
   
   if strmid(id, 0, /reverse) eq 'b' then begin
      fname = string(format='("fluxes-",I5.5,"-second.fits")', long(id))
   endif else begin
      fname = string(format='("fluxes-",A0,".fits")', id)
   endelse
   
   return, mrdfits(DIR+fname, 1)
end
