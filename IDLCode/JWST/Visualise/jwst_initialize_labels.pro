pro jwst_initialize_labels

   common linelist, lname, l_line, type

   readcol, getenv('HOME')+'/IDL/platefit/etc/linelist-jwst-full.txt', lname1, $
            l_line1, llimit, ulimit, type1, format='(A,F,F,F,A)', $
               /silent
   
   readcol, getenv('HOME')+'/IDL/platefit/etc/uv-linelist-strong.txt', lname2, $
            l_line2, lower, upper, type2, ref, format='(A,F,F,F,A,I)', $
            /silent
   type2[*] = 'is'              ; NOT CORRECT!!!
   type = [type1, type2]
   l_line = [l_line1, l_line2]
   lname = [lname1, lname2]

end
