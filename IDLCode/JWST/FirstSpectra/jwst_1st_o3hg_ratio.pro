function jwst_1st_o3hg_ratio, t, n_mc=n_mc, err_scale=err_scale
   ;;
   ;; Calculate the ([O III]4363/Hg)/([O III]5007/Hb) ratio for all
   ;; sources with Hg detected at S/N > 5
   ;;
   if n_elements(err_scale) eq 0 then err_scale = 2.0
   if n_elements(n_mc) eq 0 then n_mc = 10001



   return, jwst_1st_double_ratio(t, 'OIII_4363', 'H_GAMMA', 'OIII_5007', 'H_BETA', $
                                 n_mc=n_mc, err_scale=err_scale, norm_line='H_GAMMA', $
                                 min_sn=5.0)
end
   
;;    x1 = t.oiii_4363_flux
;;    dx1 = t.oiii_4363_flux_err*err_scale

;;    x2 = t.h_gamma_flux
;;    dx2 = t.h_gamma_flux_err*err_scale

;;    x3 = t.oiii_5007_flux
;;    dx3 = t.oiii_5007_flux_err*err_scale

;;    x4 = t.h_beta_flux
;;    dx4 = t.h_beta_flux_err*err_scale


;;    useful = where(x2/dx2 gt 5, n_useful)

;;    qs = [0.025, 0.16, 0.5, 0.84, 0.84, 0.975]
;;    ratio = fltarr(n_elements(t), n_elements(qs))
   
;;    for ii=0L, n_useful-1 do begin
;;       i = useful[ii]
 
;;       t_x1 = x1[i]+randomn(sss, n_mc)*dx1[i]
;;       t_x2 = x2[i]+randomn(sss, n_mc)*dx2[i]
;;       t_x3 = x3[i]+randomn(sss, n_mc)*dx3[i]
;;       t_x4 = x4[i]+randomn(sss, n_mc)*dx4[i]

;;       val = quantile(alog10((t_x1/t_x2)/(t_x3/t_x4)), qs)

;;       ratio[i, *] = val
      
;;    endfor

;;    return, ratio
;; end
