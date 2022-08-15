pro jwst_1st_redshift8717

   d = jwst_1st_get_coadded(8717)
   d2 = d

   cont = median(d.oned.flux, 151)
   
   d2.oned.flux=d2.oned.flux-cont
   for i=0L, 27-1 do d2.twod.flux[*, i]=d2.twod.flux[*, i]-median(d2.twod.flux[*, i], 151)


   zs = mkarr(0.4, 0.7, 0.01)
   for i=0L, n_elements(zs)-1 do begin
      jwst_vis_spec_1d2d, d2,  charsize=1.6, pcharsize=1.6, $
                          lcharsize=1.4, /zscale, range=[-6e-6, 1e-5], $
                          xrange=[3.2, 3.8], redshift=zs[i], /label, /reread
      ans = ''
      print, format='($,"z=",F5.2)', zs[i]
      read, ans

      if ans eq 'q' then break
   endfor

end
