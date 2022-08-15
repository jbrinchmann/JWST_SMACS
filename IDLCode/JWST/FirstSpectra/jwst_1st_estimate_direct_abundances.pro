pro _create_all
   o2temp = ['equal', 'Stasinska1990', 'Izotov2006', 'Pilyugin2009']

   for i=0l, n_elements(o2temp)-1 do begin
      jwst_1st_estimate_direct_abundances, /standard, o2temp=o2temp[i]
      jwst_1st_estimate_direct_abundances, o2temp=o2temp[i]
      jwst_1st_estimate_direct_abundances, /standard, o2temp=o2temp[i], /norm_to_gamma
      jwst_1st_estimate_direct_abundances, o2temp=o2temp[i], /norm_to_gamma

   endfor
   
end

pro jwst_1st_estimate_direct_abundances, standard=standard, o2temp=o2temp, $
 n_mc=n_mc, norm_to_gamma=norm_to_gamma
   ;;
   ;; Calculate direct abundances for the galaxies. 
   ;;

   if n_elements(o2temp) eq 0 then $
    o2temp = 'equal'            ; T(O2)=T(O3)

   if n_elements(n_mc) eq 0 then n_mc = 5000
   
   jwst_1st_dirs, dd
   

   ;; Calculate abundances - for this we need the fluxes.
   todo = ['04590', '05144', '06355', '08140', '10612']
   OHObj = dblarr(n_elements(todo))
   
   lines = ['OII3727', 'H_DELTA', 'OIII_4363', 'H_GAMMA', 'OIII_5007', 'H_BETA']
   i_o2 = where(lines eq 'OII_3727')
   i_hd = where(lines eq 'H_DELTA')
   i_4363 = where(lines eq 'OIII_4363')
   i_hg = where(lines eq 'H_GAMMA')
   i_5007 = where(lines eq 'OIII_5007')
   i_hb = where(lines eq 'H_BETA')
   
   for i=0L, n_elements(todo)-1 do begin
      if todo[i] eq '05144' then $
       touse = 'full' $
      else $
       touse = 'refull'

;      if i ne 2 then continue
      
      jwst_1st_line_luminosities, todo[i], lines, $
                                  lum=lum, dlum=dlum, /silent, touse=touse

      ;; Scale uncertainties
      dlum = dlum*2.0
      
      
      N_L = n_elements(lum)
      ;;
      ;; Create random values for the luminosities.
      ;;
      r_lum = dblarr(N_L, n_mc)
      for i_l=0L, N_L-1 do $
       r_lum[i_l, *] = lum[i_l] + randomn(sss, n_mc)*dlum[i]


      ;; Go through and calculate temperatures and abundances

      Te_O3 = dblarr(n_mc)
      Te_O2 = dblarr(n_mc)
      OH2 = dblarr(n_mc)
      OH3 = dblarr(n_mc)
      if keyword_set(standard) then begin
         ratio_o3 = alog10(r_lum[i_4363, *]/r_lum[i_5007, *])
         ratio_o2hb = reform(r_lum[i_o2, *]/r_lum[i_hb, *])
         ratio_o3hb = reform(r_lum[i_5007, *]/r_lum[i_hb, *])
         what = 'standard'
         o2lines = ['3727', 'Hb']
         ratio_o2 = ratio_o2hb
      end else begin
         ratio_o3 = alog10((r_lum[i_4363, *]/r_lum[i_hg, *])/(r_lum[i_5007, *]/r_lum[i_hb, *]))
         ratio_o2hd = reform(r_lum[i_o2, *]/r_lum[i_hd, *])
         ratio_o2hg = reform(r_lum[i_o2, *]/r_lum[i_hg, *])
         ratio_o3hb = reform(r_lum[i_5007, *]/r_lum[i_hb, *])
         o2lines = ['3727', 'Hd']
         o2lines_hg = ['3727', 'Hg']

         if keyword_set(norm_to_gamma) then begin
            o2lines = o2lines_hg
            ratio_o2 = ratio_o2hg
         endif else begin
            ratio_o2 = ratio_o2hd
         endelse
         
         what = 'doubleratio'
      endelse
      ratio_o3 = reform(ratio_o3)     ; Remove leading dimension 1

      
      for i_mc=0L, n_mc-1 do begin
         Te_O3[i_mc] = te_pyneb_o3(ratio_o3[i_mc], what=what)

         ;;
         ;; What do we do with the [O II] temperature? 
         ;; There are multiple options here.
         ;; Need to do .r temperature_diagnostics before running this..
         
         case o2temp of
            'equal': Te_O2[i_mc] = Te_O3[i_mc]
            'Stasinska1990': Te_O2[i_mc] = 1e4*te_oii(Te_O3[i_mc]/1e4)
            'Izotov2006': Te_O2[i_mc] = 1e4*izotov_etal_2006_t_oii(Te_O3[i_mc]/1e4, /low)
            'Pilyugin2009': Te_O2[i_mc] = 1e4*(0.264+0.835*Te_O3[i_mc]/1e4)
            else: stop
         endcase

         ;; Calculate ionic abundances for this source.
         OH2[i_mc] = o_ionic_abundance_pyneb(ratio_o2[i_mc], o2lines, Te_O2[i_mc])
         OH3[i_mc] = o_ionic_abundance_pyneb(ratio_o3hb[i_mc], ['5007', 'Hb'], Te_O2[i_mc])
      endfor

      ;; This is then not corrected for O+++! Needs to be done. 
      OH = OH2+OH3

;      stop
      
      ;;
      ;; Calculate basic stats
      ;;
      q = [0.16, 0.5, 0.84]
      use = where(OH gt 0)      ; Not all iterations worked.
      ss_OH = quantile(OH[use], q, /struct)
      ss_OH2 = quantile(OH2[use], q, /struct)
      ss_OH3 = quantile(OH3[use], q, /struct)
      ss_Te_O3 = quantile(Te_O3[use], q, /struct)
      ss_Te_O2 = quantile(Te_O2[use], q, /struct)
      
      
      ;;
      ;; Create a text report for convenience.
      ;;
      outdir = dd.root+'TeModeling/'
      outfile = outdir+'OH-'+todo[i]+'-'
      if keyword_set(standard) then $
       outfile = outfile+'classic-' $
      else $
       outfile = outfile+'double-'

      outfile = outfile+'O2temp='+o2temp
      if keyword_set(norm_to_gamma) then $
       outfile = outfile+'-Hgnorm' $
      else $
       outfile = outfile+'-Hdnorm'
      outfile = outfile+'.txt'


      f = '(A14,2X,3(E11.4,2X))'
      openw, lun, outfile, /get_lun
      printf, lun, format=f, '(O/H)+', ss_OH2.p16, ss_OH2.p50, ss_OH2.p84
      printf, lun, format=f, '(O/H)++', ss_OH3.p16, ss_OH3.p50, ss_OH3.p84
      printf, lun, format=f, '(O/H)_tot', ss_OH.p16, ss_OH.p50, ss_OH.p84
      printf, lun, format=f, 'T_e([OIII])', ss_Te_O3.p16, ss_Te_O3.p50, ss_Te_O3.p84
      printf, lun, format=f, 'T_e([OII])', ss_Te_O2.p16, ss_Te_O2.p50, ss_Te_O2.p84
      free_lun, lun

      OHobj[i] = 12+alog10(ss_OH.p50)

      print, todo[i], OHobj[i]
      
   endfor

      
;   stop
   
   

end
