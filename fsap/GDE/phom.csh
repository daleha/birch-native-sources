#!/bin/csh

#Version 11/23/02
# Run P1HOM, P2HOM as a command
#Synopsis: phom.csh seqfilex seqfiley outfile startx finishx starty finishy\
#                    range minper compfact linesize program

#Convert arguments to variables
set SEQFILEX = $1
set SEQFILEY = $2
set OUTFILE = $3
set STARTX  = $4
set FINISHX = $5
set STARTY  = $6
set FINISHY = $7
set RANGE   = $8
set MINPER  = $9
set COMPFACT = $10
set LINESIZE = $11
set PROGRAM = $12

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

   # Determine the size of seq. in seqfilex  and seqfiley and make sure that
   # start and finish values aren't greater than seq. length.
   tail +2 $SEQFILEX > TEMP.$PID
   @ SEQLENGTHX = (`wc -c < TEMP.$PID` - `wc -l < TEMP.$PID`) - 1
   if ($STARTX > $SEQLENGTHX) set STARTX = $SEQLENGTHX
   if ($FINISHX > $SEQLENGTHX) set FINISHX = $SEQLENGTHX
   
   tail +2 $SEQFILEY > TEMP.$PID
   @ SEQLENGTHY = (`wc -c < TEMP.$PID` - `wc -l < TEMP.$PID`) - 1
   if ($STARTY > $SEQLENGTHY) set STARTY = $SEQLENGTHY
   if ($FINISHY > $SEQLENGTHY) set FINISHY = $SEQLENGTHY
   
   
   # Abort if parameters aren't set up properly.
   if ($OKAY == false) then
      echo '>>> Aborting program.'
      exit
   
   else #----------------- generate keyboard input to send to program -----
   
   # Open files
   echo $SEQFILEX  > $PID.COMFILE
   echo n          >> $PID.COMFILE
   echo $SEQFILEY  >> $PID.COMFILE 
   echo n          >> $PID.COMFILE 
   echo $OUTFILE   >> $PID.COMFILE 
   
   # Set parameters
   echo 4          >> $PID.COMFILE #Choose parameter menu
   
   echo 1          >> $PID.COMFILE #choose startx
   echo $STARTX    >> $PID.COMFILE 
   
   echo 2          >> $PID.COMFILE  #choose finishx
   echo $FINISHX   >> $PID.COMFILE 
   
   echo 3          >> $PID.COMFILE #choose starty
   echo $STARTY    >> $PID.COMFILE 
   
   echo 4          >> $PID.COMFILE  #choose finishy
   echo $FINISHY   >> $PID.COMFILE 
   
   echo 5          >> $PID.COMFILE #choose range 
   echo $RANGE     >> $PID.COMFILE 
     
   echo 7          >> $PID.COMFILE  #choose minper 
   echo $MINPER    >> $PID.COMFILE 
   
   echo 8          >> $PID.COMFILE  #choose compfact 
   echo $COMPFACT  >> $PID.COMFILE 
   
   echo 9          >> $PID.COMFILE  #choose linesize
   echo $LINESIZE  >> $PID.COMFILE 
   
   echo 0          >> $PID.COMFILE  #exit parameter menu
   
   echo 6          >> $PID.COMFILE #write output to file
   
   echo 0          >> $PID.COMFILE  #exit program
   
   endif

   $PROGRAM < $PID.COMFILE
   $RM_CMD *.$PID $PID.COMFILE
endif
