function jwst_summarise_one_cl01, r_in, dim
   ;;
   ;; Given the results of a CL01 fit, summarise it. 
   ;;

   todo = ['MUGAS', 'OH', 'DUSTTOGAS', 'TAUV', 'U', 'XSI', 'Z', 'LOGSFR', 'LOG_TR', 'TE']
   dimname = ['XMUGAS', 'XOH', 'XDUSTTOGAS', 'TAUV', 'LOGU', 'XSI', 'LOGZ', $
              'XLOGSFR', 'XLOG_TR', 'BIN_TE']


   r = r_in
   r = create_struct(r, dim)

   res = {todo: todo}
   for i_quantity=0L, n_elements(todo)-1 do begin
      ;;
      ;;---------------------
      ;; Get the fit results
      ;;---------------------
      print, 'Doing '+dimname[i_quantity]

      x = struct_var(r, dimname[i_quantity])
      y = struct_var(r, 'BIN_'+todo[i_quantity])

;      stop
      ss = summarise_distribution(x, y, /silent, $
                                  missing=-1.0D0)
      
      res = create_struct(res, todo[i_quantity], ss)
      
   endfor

   return, res
   
   

end
