#!/bin/csh

#Version 11/23/02
# Run D3HOM, D4HOM as a command
#Synopsis: dxhom.csh seqfilex strand seqfiley outfile startx finishx starty\
#            finishy range minper compfact kind linesize program

#Convert arguments to variables
set SEQFILEX = $1
set STRAND  = $2
set SEQFILEY = $3
set OUTFILE = $4
set STARTX  = $5
set FINISHX = $6
set STARTY  = $7
set FINISHY = $8
set RANGE   = $9
set MINPER  = $10
set COMPFACT = $11
set KIND     = $12
set LINESIZE = $13
set PROGRAM = $14

#!/bin/csh

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

#-----------   error checking ---------------------
if ((-e $SEQFILEX) && (! -z $SEQFILEX) && (-e $SEQFILEY) && (! -z $SEQFILEY)) then
   set OKAY = true
   set PID = $$ #process id

   # Correct for bug in GDE2.2, which writes 'Circular' in wrong columns
   # in GenBank LOCUS line. This is done by determining whether seq is
   # linear or circular, converting to Pearson format, and appending the
   # '1' if linear or '2' if circular. We also have to work around a bug
   # in readseq that causes it to write 'Circular' on the DEFINITION line
   # ARRRGH!
   foreach file ($SEQFILEX $SEQFILEY)
      head -2 $file |grep 'Circular' > circ.$PID
      if (-z circ.$PID) then
         set TOPOLOGY = 1
      else 
         set TOPOLOGY = 2
      endif
#      readseq  -i1 -f8 $file >$file.$PID
      readseq  -i1 -ffasta -o=$file.$PID $file

      echo $TOPOLOGY >> $file.$PID
   end
   
   # Determine the size of seq. in seqfilex  and seqfiley and make sure that
   # start and finish values aren't greater than seq. length.
   tail +2 $SEQFILEX.$PID > TEMP.$PID
   @ SEQLENGTHX = (`wc -c < TEMP.$PID` - `wc -l < TEMP.$PID`) - 1
   if ($STARTX > $SEQLENGTHX) set STARTX = $SEQLENGTHX
   if ($FINISHX > $SEQLENGTHX) set FINISHX = $SEQLENGTHX
   
   tail +2 $SEQFILEY.$PID > TEMP.$PID
   @ SEQLENGTHY = (`wc -c < TEMP.$PID` - `wc -l < TEMP.$PID`) - 1
   if ($STARTY > $SEQLENGTHY) set STARTY = $SEQLENGTHY
   if ($FINISHY > $SEQLENGTHY) set FINISHY = $SEQLENGTHY
   
   
   # Abort if parameters aren't set up properly.
   if ($OKAY == false) then
      echo '>>> Aborting program.'
      exit
   
   else #----------------- generate keyboard input to send to program -----

   # Open files
   echo $SEQFILEX.$PID > $PID.COMFILE
   echo b              >> $PID.COMFILE
   echo $STRAND        >> $PID.COMFILE
   echo $SEQFILEY.$PID >> $PID.COMFILE
   echo b              >> $PID.COMFILE
   echo $OUTFILE       >> $PID.COMFILE
   
   # Set parameters
   echo 4        >> $PID.COMFILE #Choose parameter menu
   
   echo 1        >> $PID.COMFILE #choose startx
   echo $STARTX  >> $PID.COMFILE 
   
   echo 2        >> $PID.COMFILE  #choose finishx
   echo $FINISHX >> $PID.COMFILE 
   
   echo 3        >> $PID.COMFILE  #choose startx
   echo $STARTY  >> $PID.COMFILE 
   
   echo 4        >> $PID.COMFILE  #choose finishx
   echo $FINISHY >> $PID.COMFILE 
   
   echo 5        >> $PID.COMFILE  #choose range 
   echo $RANGE   >> $PID.COMFILE 
     
   echo 7        >> $PID.COMFILE  #choose minper 
   echo $MINPER  >> $PID.COMFILE 
   
   echo 8        >> $PID.COMFILE  #choose compfact 
   echo $COMPFACT >> $PID.COMFILE 
   
   echo 9        >> $PID.COMFILE  #choose kind
   echo $KIND    >> $PID.COMFILE 
   
   echo 10       >> $PID.COMFILE  #choose linesize
   echo $LINESIZE >> $PID.COMFILE 
   
   echo 0        >> $PID.COMFILE   #exit parameter menu
   
   echo 6        >> $PID.COMFILE  #write output to file
   
   echo 0        >> $PID.COMFILE   #exit program
   
   endif

   # run d3hom or d4hom and then clean up temporary files
   $PROGRAM < $PID.COMFILE
   $RM_CMD TEMP.$PID circ.$PID $PID.COMFILE
endif
