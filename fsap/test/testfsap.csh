#!/bin/csh
# TESTFSAP - a c-shell script (version 29 Mar 2006)
# Run FSAP programs on test datasets.
# 
# 'testfsap.csh' with no arguments tests all FSAP programs.
# 'testfsap.csh program_name' tests the one program
# eg. 'testfsap.csh numseq' tests numseq

#####################################################################
# set default parameters
#####################################################################

# RM_CMD - command to be used for removing files and directories
if (-e /usr/bin/rm) then
   set RM_CMD = /usr/bin/rm
else
   if (-e /bin/rm) then
      set RM_CMD = /bin/rm
   else
      set RM_CMD = rm
   endif
endif
#set BINDIR = "../bin-solaris-sparc"
#set BINDIR = "../bin-solaris-amd64"
#set BINDIR = "../bin-linux-intel"


# set $testset to programs to be tested.
set NUMARGS = $#argv

if ($NUMARGS == 0) then 
   set TESTSET = (bachrest d3hom d4hom multidigest funnel gel intrest numseq\
       p1hom p2hom prostat testcode)
else
   set TESTSET = $1
endif

#####################################################################
# test programs one at a time
#####################################################################
foreach PROGRAM ($TESTSET)
  echo "Testing $PROGRAM ..."
  # run the program
  if (-e temp.$PROGRAM.out) $RM_CMD -f temp.$PROGRAM.out
  $BINDIR/$PROGRAM <$PROGRAM.inp

  # if successful, compare sample files with new output
  set SUCCESS = $status
  if ($SUCCESS == 0) then
     echo "     $PROGRAM successful"
     echo "     Comparing sample files with new output (temp.*)..."
     diff $PROGRAM.out temp.$PROGRAM.out >result
     if (-z result) then
        echo "     output files compare okay"
        $RM_CMD -f temp.$PROGRAM.out result
     else
        echo ">>>> $PROGRAM.out is different from temp.$PROGRAM.out"
        echo ">>>> see difference.$PROGRAM"
        mv result difference.$PROGRAM
     endif 
   else
     echo ">>>> $PROGRAM failed"
  endif
end #foreach

echo TEST OF FSAP COMPLETE
