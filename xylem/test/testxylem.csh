#!/bin/csh
# TESTXYLEM - a c-shell script (version 18 Aug 01)
# Run XYLEM programs on test datasets.
# 
# 'textxylem' with no arguments tests all xylem programs.
# 'testxylem program_name' tests the one program
# eg. 'testxylem splitdb' tests splitdb

#####################################################################
# set default parameters
#####################################################################
set RM = rm
set BINDIR = ../bin

# set $testset to programs to be tested.
set NUMARGS = $#argv

if ($NUMARGS == 0) then 
   set TESTSET = (splitdb getloc identify getob reform ribosome shuffle\
                  dbstat flat2phyl findkey fetch features)
else
   set TESTSET = $1
endif

#####################################################################
# test programs one at a time
#####################################################################
foreach PROGRAM ($TESTSET)
  echo "Testing $PROGRAM ..."
  switch ($PROGRAM)
    ######### splitdb #######
    case splitdb:
       # run the program
       $BINDIR/splitdb sample.gen temp.ano temp.wrp temp.ind

       # if successful, compare sample files with new output
       set SUCCESS = $status
       if ($SUCCESS == 0) then
          echo "     splitdb successful"
          echo "     Comparing sample files with new output (temp.*)..."
          diff sample.ano temp.ano >result
          if (-z result) then
             echo "     .ano files compare okay"
             $RM temp.ano result
          else
             echo ">>>> sample.ano is different from temp.ano"
             echo ">>>> see difference.ano"
             mv result difference.ano
          endif 
          diff sample.wrp temp.wrp >result
          if (-z result) then
             echo "     .wrp files compare okay"
             $RM temp.wrp result
          else
             echo ">>>> sample.wrp is different from temp.ano"
             echo ">>>> see difference.wrp"
             mv result difference.wrp
          endif
          diff sample.ind temp.ind >result
          if (-z result) then
             echo "     .ind files compare okay"
             $RM temp.ind result
          else
             echo ">>>> sample.ind is different from temp.ind"
             echo ">>>> see difference.ind"
             mv result difference.ind
          endif 
        else
          echo ">>>> splitdb failed"
       endif
      breaksw

    ######### getloc #######
    case getloc:
       # run the program
       $BINDIR/getloc sample.nam sample.ano sample.wrp sample.ind temp.gen

       # if successful, compare sample files with new output
       set SUCCESS = $status
       if ($SUCCESS == 0) then
          echo "     getloc successful"
          echo "     Comparing sample files with new output (temp.*)..."
          diff sample.gen temp.gen >result
          if (-z result) then
             echo "     .gen files compare okay"
             $RM temp.gen result
          else
             echo ">>>> sample.gen is different from temp.gen"
             echo ">>>> see difference.gen"
             mv result difference.gen
          endif 
        else
          echo ">>>> getloc failed"
       endif
      breaksw

    ######### identify #######
    case identify:
       # run the program
       egrep -n -i LOCUS sample.ano >temp.egrep
       $BINDIR/identify temp.egrep sample.ind temp.nam temp.fnd

       # if successful, compare sample files with new output
       set SUCCESS = $status
       if ($SUCCESS == 0) then
          echo "     identify successful"
          echo "     Comparing sample files with new output (temp.*)..."
          diff sample.nam temp.nam >result
          if (-z result) then
             echo "     .nam files compare okay"
             $RM temp.nam result
          else
             echo ">>>> sample.nam is different from temp.nam"
             echo ">>>> see difference.nam"
             mv result difference.nam
          endif 
          diff sample.fnd temp.fnd >result
          if (-z result) then
             echo "     .fnd files compare okay"
             $RM temp.fnd result
          else
             echo ">>>> sample.nam is different from temp.nam"
             echo ">>>> see difference.fnd"
             mv result difference.fnd
          endif
        else
          echo ">>>> identify failed"
       endif
       $RM temp.egrep
      breaksw

    ######### getob #######
    case getob:
       # run the program
       $BINDIR/getob sample.inp sample.nam sample.ano sample.wrp sample.ind\
                     temp.getob.msg temp.getob.out temp.getob.exp

       # if successful, compare sample files with new output
       set SUCCESS = $status
       if ($SUCCESS == 0) then
          echo "     getob successful"
          echo "     Comparing sample files with new output (temp.*)..."
          diff sample.msg temp.getob.msg >result
          if (-z result) then
             echo "     .msg files compare okay"
             $RM temp.getob.msg result
          else
             echo ">>>> sample.msg is different from temp.getob.msg"
             echo ">>>> see difference.getob.msg"
             mv result difference.getob.msg
          endif 
          diff sample.out temp.getob.out >result
          if (-z result) then
             echo "     .out files compare okay"
             $RM temp.getob.out result
          else
             echo ">>>> sample.out is different from temp.getob.out"
             echo ">>>> see difference.getob.out"
             mv result difference.getob.out
          endif 
          diff sample.exp temp.getob.exp >result
          if (-z result) then
             echo "     .exp files compare okay"
             $RM temp.getob.exp result
          else
             echo ">>>> sample.exp is different from temp.getob.exp"
             echo ">>>> see difference.getob.exp"
             mv result difference.getob.exp
          endif 
        else
          echo ">>>> getob failed"
       endif
      breaksw

    ######### reform #######
    case reform:
       # run the program
       $BINDIR/reform -g -p -c -fi <sample.aln > temp.ref

       # if successful, compare sample files with new output
       set SUCCESS = $status
       if ($SUCCESS == 0) then
          echo "     reform successful"
          echo "     Comparing sample files with new output (temp.*)..."
          diff sample.ref temp.ref >result
          if (-z result) then
             echo "     .ref files compare okay"
             $RM temp.ref result
          else
             echo ">>>> sample.ref is different from temp.ref"
             echo ">>>> see difference.ref"
             mv result difference.ref
          endif 
        else
          echo ">>>> reform failed"
       endif
      breaksw

    ######### ribosome #######
    case ribosome:
       # run the program
       $BINDIR/ribosome <sample.out > temp.pep

       # if successful, compare sample files with new output
       set SUCCESS = $status
       if ($SUCCESS == 0) then
          echo "     ribosome successful"
          echo "     Comparing sample files with new output (temp.*)..."
          diff sample.pep temp.pep >result
          if (-z result) then
             echo "     .pep files compare okay"
             $RM temp.pep result
          else
             echo ">>>> sample.pep is different from temp.pep"
             echo ">>>> see difference.pep"
             mv result difference.pep
          endif 
        else
          echo ">>>> ribosome failed"
       endif
      breaksw

    ######### prot2nuc #######
    case prot2nuc:
       # run the program
       $BINDIR/prot2nuc <sample.pro > temp.prot2nuc

       # if successful, compare sample files with new output
       set SUCCESS = $status
       if ($SUCCESS == 0) then
          echo "     prot2nuc successful"
          echo "     Comparing sample files with new output (temp.*)..."
          diff sample.prot2nuc temp.prot2nuc >result
          if (-z result) then
             echo "     .prot2nuc files compare okay"
             $RM temp.prot2nuc result
          else
             echo ">>>> sample.prot2nuc is different from temp.prot2nuc"
             echo ">>>> see difference.prot2nuc"
             mv result difference.prot2nuc
          endif 
        else
          echo ">>>> prot2nuc failed"
       endif
      breaksw

     ######### dbstat #######
    case dbstat:
       # run the program
       $BINDIR/dbstat <sample.pep > temp.dbstat

       # if successful, compare sample files with new output
       set SUCCESS = $status
       if ($SUCCESS == 0) then
          echo "     dbstat successful"
          echo "     Comparing sample files with new output (temp.*)..."
          diff sample.dbstat temp.dbstat >result
          if (-z result) then
             echo "     .dbstat files compare okay"
             $RM temp.dbstat result
          else
             echo ">>>> sample.dbstat is different from temp.dbstat"
             echo ">>>> see difference.dbstat"
             mv result difference.dbstat
          endif 
        else
          echo ">>>> dbstat failed"
       endif
      breaksw

     ######### flat2phyl #######
    case flat2phyl:
      
       # check results for interleaved format
       $BINDIR/flat2phyl <sample.flat > temp.interleaved.phylip

       # if successful, compare sample files with new output
       set SUCCESS = $status
       if ($SUCCESS == 0) then
          echo "     flat2phyl successful for interleaved format"
          echo "     Comparing sample files with new output (temp.*)..."

          diff sample.interleaved.phylip temp.interleaved.phylip >result
          if (-z result) then
             echo "     flat2phyl interleaved files compare okay"
             $RM temp.interleaved.phylip result
          else
             echo ">>>> sample.interleaved.phyip is different from\
                 temp.interleaved.phylip"
             echo ">>>> see difference.flat2phyl.interleaved"
             mv result difference.flat2phyl.interleaved
          endif 
        else
          echo ">>>> flat2phyl failed for interleaved format"
       endif

       # check results for sequential format
       $BINDIR/flat2phyl -s <sample.flat > temp.sequential.phylip

       # if successful, compare sample files with new output
       set SUCCESS = $status
       if ($SUCCESS == 0) then
          echo "     flat2phyl successful for sequential format"
          echo "     Comparing sample files with new output (temp.*)..."

          diff sample.sequential.phylip temp.sequential.phylip >result
          if (-z result) then
             echo "     flat2phyl sequential files compare okay"
             $RM temp.sequential.phylip result
          else
             echo ">>>> sample.sequential.phyip is different from\
                 temp.sequential.phylip"
             echo ">>>> see difference.flat2phyl.sequential"
             mv result difference.flat2phyl.sequential
          endif 
        else
          echo ">>>> flat2phyl failed for sequential format"
       endif
      breaksw

   ######### shuffle #######
    case shuffle:
       # run the program
       $BINDIR/shuffle -s3456 <sample.out > temp.shuf

       # if successful, compare sample files with new output
       set SUCCESS = $status
       if ($SUCCESS == 0) then
          echo "     shuffle successful"
          $RM temp.shuf
        else
          echo ">>>> shuffle failed"
       endif
      breaksw

    ######### findkey #######
    case findkey:
       # run the program
       echo "LOCUS" > temp.kw
       $BINDIR/findkey -G sample.ano temp.kw temp.nam temp.fnd

       # if successful, compare sample files with new output
       set SUCCESS = $status
       if ($SUCCESS == 0) then
          echo "     findkey successful"
          echo "     Comparing sample files with new output (temp.*)..."
          diff sample.nam temp.nam >result
          if (-z result) then
             echo "     .nam files compare okay"
             $RM temp.nam result
          else
             echo ">>>> sample.nam is different from temp.nam"
             echo ">>>> see difference.nam"
             mv result difference.nam
          endif 
          diff sample.fnd temp.fnd >result
          if (-z result) then
             echo "     .fnd files compare okay"
             $RM temp.fnd result
          else
             echo ">>>> sample.nam is different from temp.nam"
             echo ">>>> see difference.fnd"
             mv result difference.fnd
          endif
        else
          echo ">>>> findkey failed"
       endif
       $RM temp.kw
      breaksw

    ######### fetch #######
    case fetch:
       # run the program
       $BINDIR/fetch -G sample sample.nam temp.gen

       # if successful, compare sample files with new output
       set SUCCESS = $status
       if ($SUCCESS == 0) then
          echo "     fetch successful"
          echo "     Comparing sample files with new output (temp.*)..."
          diff sample.gen temp.gen >result
          if (-z result) then
             echo "     .gen files compare okay"
             $RM temp.gen result
          else
             echo ">>>> sample.gen is different from temp.gen"
             echo ">>>> see difference.gen"
             mv result difference.gen
          endif 
        else
          echo ">>>> fetch failed"
       endif
      breaksw

    ######### features #######
    case features:
       # run the program
       cp sample.exp temp.exp
       $BINDIR/features -u sample -E temp.exp 

       # if successful, compare sample files with new output
       set SUCCESS = $status
       if ($SUCCESS == 0) then
          echo "     features successful"
          echo "     Comparing sample files with new output (temp.*)..."
          diff sample.out temp.out >result
          if (-z result) then
             echo "     .out files compare okay"
             $RM temp.out result
          else
             echo ">>>> sample.out is different from temp.out"
             echo ">>>> see difference.out"
             mv result difference.out
          endif 
        else
          echo ">>>> features failed"
       endif
       $RM temp.msg temp.exp
      breaksw

    endsw

end #foreach

echo TEST OF XYLEM COMPLETE
