#!/bin/csh

#Version 11/ 5/00
# Convert free format file to pseudo GenBank format
# to be read by GDE.
#Synopsis: free2gb.csh infile outfile
#Convert arguments to variables
set INFILE  = $1
set OUTFILE = $2

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

   set PID = $$ #process id

   
   #run funnel to delete non-sequence characters
   echo $INFILE          > $PID.COMFILE #input filename
   echo $PID.raw         >> $PID.COMFILE #outfile
   echo 50                >> $PID.COMFILE # print 50 nt per line
   funnel < $PID.COMFILE

   # Create a Fasta format file for input to readseq.
   echo '>'$INFILE > $PID.wrp
   # copy any non-comment lines to end of fasta file
   grep -v ';' < $PID.raw >> $PID.wrp
   # add a blank line to end of file. A bug in the old
   # readseq loses some characters from the end of file.
   # This is not a problem in the new Java readseq.
   echo ' ' >> $PID.wrp
   # convert to pseudo GenBank format
   readseq -a -fGenBank -o=$OUTFILE $PID.wrp

   # delete temporary files
   $RM_CMD $PID.*

endif
