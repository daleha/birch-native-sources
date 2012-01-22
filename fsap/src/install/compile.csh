# Compile all pascal programs with names in the form X.p
# Write the final executable code to X.
# Write all compiler messages to X.msg
# By default, X is written to the directory ../../bin 
# To write executable files to another directory,
# change bindir to the name of that directory.
set bindir = ../../bin 
 
foreach i (*.p)
   set name = $i:r

# use this form if programs will only be run on a machine on which
# Pascal libraries reside (library routines automatically linked at
# runtime, produces smaller code)

#   pc -L -O4 -o $bindir/$name $i > $i.msg

#  use this form if programs will be run on workstations that do not
#  have Pascal libraries (produces larger executable codefiles that
#  have all library routines pre-linked)
    pc -L -Bstatic -O4 -o $bindir/$name $i > $i.msg

   chmod a+rx $bindir/$name
end
