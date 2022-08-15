function jwst_1st_load_pifit_pdf, object, what, include_ne3=include_ne3, gridded=gridded
   ;;
   ;; Load an output PDF from the PIFit runs.
   ;;

   if keyword_set(include_ne3) then $
    infix = '-w-ne3' $
   else $
    infix = ''
   


   jwst_1st_dirs, dd

   if keyword_set(gridded) then begin
      DIR = dd.root+'PIFit/G16-grid/PDFs/'
   fname = string(format='("pdf-",I5.5,"-",A,".fits")', $
                  object, what)
   endif else begin
      DIR = dd.root+'PIFit/G16/PDFs/'
      fname = string(format='("PDF-G16-",I5.5,"-",A,A,".fits")', $
                     object, what, infix)
   endelse
   

   fname = DIR+fname

   return, mrdfits(fname, 1, /silent)


end
