pro _fit_all_highz
   
   
   todo = ['04590', '05144', '06355', '08140', '10612']

   d = Dictionary()
   s = Dictionary()
   
   pos = get_position_arr(5, nx=5, ny=1)
   for i=0L, n_elements(todo)-1 do begin
      jwst_1st_fit_cl01, todo[i], result, ss

      
      d['id'+todo[i]] = result
      s['id'+todo[i]] = ss 
      pos = get_position_arr(i)
      plot, result.xoh, result.bin_oh, psym=10, $
            position=pos, noerase=(i ne 0), ytickformat='noticks', $
            title=todo[i]

      print, '12 + Log O/H', ss.oh.median
      
      
   endfor

   jwst_1st_dirs, dd
   outdir = dd.root+'CL01Fit/'

   save, d, s, file=outdir+'combined-results.sav'
   
end

pro jwst_1st_fit_cl01, object, results, ss


   
   if size(object, /tname) ne 'STRING' then $
    object = string(format='(I0)', object)
   jwst_1st_dirs, dd
   t = mrdfits(dd.ROOT+'redshift-for-sources.fits', 1)

   i_obj = where(strtrim(t.object) eq '02736_'+object, n_match)
   if (n_match eq 0) then stop

   z = t[i_obj].redshift


   lines_to_use = ['LOII_3727', 'LHg_4340', 'LHb_4861', 'LOIII_5007', $
                   'LHa_6563', 'LNII_6584']
   lines = ['OII3727', 'H_GAMMA', 'H_BETA', $
                                        'OIII_5007', 'H_ALPHA', 'NII_6584']
   jwst_1st_line_luminosities, object, lines, $
                               lum=l, dlum=dl, /silent

   use = where(l gt -1e4)
   l = l[use]
   dl = dl[use]
   lines = lines[use]
   lines_to_use = lines_to_use[use]

   print, 'I will fit using these lines: ', lines
;   stop
   
   jwst_fit_one_cl01, z, l, dl, models=models, dim=dim, oh=oh, te=te, $
                      etaha=etaha, etao2=etao2, mugas=mugas, dgr=dgr, $
                      results=results, lines=lines_to_use, apply_tauv_prior=1

   ss = jwst_summarise_one_cl01(results, dim)

   
   outdir = dd.root+'CL01Fit/'
   struct_to_hdf5, dim, outdir+'dim.h5'
   
   outfile = outdir+'cl01fit-'+object+'.h5'
   struct_to_hdf5, results, outfile

   save, ss, file=outdir+'summary-cl01fit-'+object+'.h5'
   
end
