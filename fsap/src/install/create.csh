# For any file X.rp, add modules and write the complete, compileable 
# source code to X.p. X.lst contains a listing of all modules added.
# The files modcat and modlib must exist in the current 
# directory. Also, all files to be created must be in
# current directory and have the file suffix .rp
foreach i (*.rp)
   cat $i > raw
   module
   set name = $i:r
   mv sout $name.p
   mv list $name.lst
   /usr/bin/rm raw
end
