pro jwst_1st_feii_hei_model_grid, m, model_par1=model_par1, model_par2=model_par2, $
                                  color1=color1, color2=color2, overplot=overplot, $
                                  _extra=_extra
   ;;
   ;; Show a simple model grid given an input.
   ;;

   if (n_elements(color1) eq 0) then $
    color1 = 'red'
   if (n_elements(color2) eq 0) then $
    color2 = 'CornFlowerBlue'

   m_pars = m.model_pars
   if (n_elements(m_pars) gt 2) then begin
      stop                      ; Needs handling
   endif else begin
      par1 = m_pars[0]
      par2 = m_pars[1]
   endelse

   p1 = struct_var(m, par1)
   p2 = struct_var(m, par2)

   ;; Find the unique parameter values
   si1 = sort(p1)
   ui1 = uniq(p1[si1])
   unique_p1 = p1[si1[ui1]]
   si2 = sort(p2)
   ui2 = uniq(p2[si2])
   unique_p2 = p2[si2[ui2]]

   x = alog10(m.FEII_1_2570UM/m.HI_1_2818UM)
   y = alog10(m.HEI_1_0833UM/m.HI_1_0938UM)

   if not keyword_set(overplot) then begin
      plot, x, y, /nodata, xtitle=TeXtoIDL('Log [Fe II] 1.257/Pa-\beta'), $
            ytitle=TeXtoIDL('Log He I 1.083/Pa-\gamma'), _extra=_extra
   endif
   
   for i=0L, n_elements(unique_p1)-1 do begin
      ii = where(p1 eq unique_p1[i], n)
      oplot, x[ii], y[ii], color=jb_colour(color1)

      
      xyouts, x[ii[n-1]]+0.01, y[ii[n-1]], string(format='(A0,"=",F0.1)', m_pars[0], unique_p1[i]), $
              color=jb_colour(color1)
   endfor

   for i=0L, n_elements(unique_p2)-1 do begin
      ii = where(p2 eq unique_p2[i], n)
      oplot, x[ii], y[ii], color=jb_colour(color2)
      xyouts, x[ii[0]]+0.01, y[ii[0]], string(format='(A0,"=",F0.1)', m_pars[1], unique_p2[i]), $
              align=0.5, color=jb_colour(color2)
   endfor


end
