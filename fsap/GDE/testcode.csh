#!/bin/csh

#Version 11/23/02
# Run TESTCODE as a command
#Synopsis: testcode.csh infile outfile start finish which format window\
#                          skip

#Convert arguments to variables
set INFILE  = $1
set OUTFILE = $2
set START   = $3
set FINISH  = $4
set WHICH   = $5
set FORMAT  = $6
set WINDOW  = $7
set SKIP    = $8

set PID = $$ #process id

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

# Abort if INFILE does not exist or is of zero length
if ((-e $INFILE) && (! -z $INFILE)) then

   # Correct for bug in GDE2.2, which writes 'Circular' in wrong columns
   # in GenBank LOCUS line. This is done by determining whether seq is
   # linear or circular, converting to Pearson format, and appending the
   # '1' if linear or '2' if circular.
   head -1 $INFILE |grep 'Circular' > circ.$PID
   if (-z circ.$PID) then
      set TOPOLOGY = 1
   else 
      set TOPOLOGY = 2
   endif
   readseq  -i1 -f8  -o=$INFILE.$PID $INFILE
   echo $TOPOLOGY >> $INFILE.$PID
   
   # Determine the size of seq. ininfile and make sure that
   # start and finish values aren't greater than seq. length.
   tail +2 $INFILE.$PID > TEMP.$PID
   @ SEQLENGTH = (`wc -c < TEMP.$PID` - `wc -l < TEMP.$PID`) - 1
   if ($START > $SEQLENGTH) set START = $SEQLENGTH
   if ($FINISH > $SEQLENGTH) set FINISH = $SEQLENGTH
   @ MAXCODONS = ($SEQLENGTH / 3)
   if ($WINDOW > $MAXCODONS) set WINDOW = $MAXCODONS
   if ($SKIP > $MAXCODONS) set SKIP = $MAXCODONS
   
   #initial filenames 
   echo $INFILE.$PID > $PID.COMFILE #input filename
   echo b            >> $PID.COMFILE #Pearson format
   echo $OUTFILE     >> $PID.COMFILE #outfile
   
   # Set parameters
   echo 4       >> $PID.COMFILE #Choose parameter menu
   
   echo 1       >> $PID.COMFILE #choose start
   echo $START
   
   echo 2       >> $PID.COMFILE #choose finish
   echo $FINISH >> $PID.COMFILE 
   
   echo 3       >> $PID.COMFILE #choose which 
   echo $WHICH  >> $PID.COMFILE 
   
   echo 4       >> $PID.COMFILE #choose format 
   echo $FORMAT >> $PID.COMFILE 
   
   echo 5       >> $PID.COMFILE #choose window 
   echo $WINDOW >> $PID.COMFILE 
   
   echo 6       >> $PID.COMFILE #choose skip
   echo $SKIP   >> $PID.COMFILE 
   
   echo 0       >> $PID.COMFILE #exit parameter menu
   
   echo 6       >> $PID.COMFILE #Print output to file 
   echo ""      >> $PID.COMFILE #dummy prompt line
   echo 0       >> $PID.COMFILE #exit program
   
   testcode < $PID.COMFILE
   $RM_CMD *.$PID $PID.COMFILE

endif
