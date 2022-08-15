;;
;; Compare the platefits of the individual spectra. 
;;
;; In this case I ran platefit on each observation and here I want to
;; compare the re-measurements. I will assume that any discrepancy is
;; due to noise and that this is independent of the line used but
;; dependent on the noise estimate adopted. 
;;

pro jwst_1st_analyse_individual_pfits, sncut=sncut, pipeline_error=pipeline_error

   if n_elements(sncut) eq 0 then sncut = 10.0

   if keyword_set(pipeline_error) then suffix = '-pipeline-err' else suffix = ''

   jwst_1st_dirs, dd
   pf1 = mrdfits(dd.SPECDIR+'Platefit/individual-obs07'+suffix+'.fit', 1, /silent)
   pf2 = mrdfits(dd.SPECDIR+'Platefit/individual-obs08'+suffix+'.fit', 1, /silent)

   ;; The lines used
   readcol, dd.SPECDIR+'Platefit/linelist-fixed-strong.txt', $
            lname, lcentr, llow, lhigh, type, format='A,F,F,F,A', /silent

   keep = where(type eq 'em')
   lname = lname[keep]
   lcentr = lcentr[keep]
;   stop
   
   n_lines = n_elements(lname)
   
   ;; Go through and assemble lines where both are high enough S/N.
   flux1 = []
   flux2 = []
   dflux1 = []
   dflux2 = []
   lobs = []
   lrest = []
   lused = []
   for i=0L, n_lines-1 do begin
      ln = lname[i]+'_FLUX'
      x1 = struct_var(pf1, ln)
      x2 = struct_var(pf2, ln)
      dx1 = struct_var(pf1, ln+'_ERR')
      dx2 = struct_var(pf2, ln+'_ERR')
      
      good = where(x1/dx1 gt sncut and x2/dx2 gt sncut, n_good)
;      stop

      if n_good gt 0 then begin
         print, format='(A12," gives ",I2," high S/N lines")', $
                lname[i], n_good
         flux1 = [flux1, x1[good]]
         flux2 = [flux2, x2[good]]
         dflux1 = [dflux1, dx1[good]]
         dflux2 = [dflux2, dx2[good]]

         lused = [lused, replicate(lname[i], n_good)]
         this_lobs = lcentr[i]*(1+pf1[good].z)
         lobs = [lobs, this_lobs]
         lrest = [lrest, replicate(lcentr[i], n_good)]
      endif

   endfor
   mean_sn = 0.5*(flux1/dflux1+flux2/dflux2)
   ndiff = (flux1-flux2)/sqrt(dflux1^2+dflux2^2)
   ok =where(mean_sn lt 100, n_ok)

   histogauss, ndiff[ok], a;, /noplot
   print, ''
   print, 'Mean of Gaussian  (007-008)=', a[1]
   print, 'Sigma of Gaussian=', a[2]
   print, format='("I used ",I0," galaxies.")',  n_ok
   
;   stop
end
