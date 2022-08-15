function jwst_1st_load_stamp, object, filter, orig=orig
   ;;
   ;; Get a stamp
   ;;
   jwst_1st_dirs, dd

   if keyword_set(orig) then suffix = '-orig' $
   else suffix = ''

   
   filter = strlowcase(filter)
   if filter eq 'f435w' or filter eq 'f606w' or filter eq 'f814w'  then $
    filt_name = 'HST_'+filter $
   else filt_name = filter

   fname = dd.stampdir+'2736_'+string(format='(I0)', object)+'_'+ $
           filt_name+suffix+'.fits'

   im = mrdfits(fname, 0, hdr)

   return, {im: im, header: hdr}
end
   
