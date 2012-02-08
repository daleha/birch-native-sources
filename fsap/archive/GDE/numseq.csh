#!/bin/csh

#Version  11/23/02
# Run NUMSEQ as a command
#Synopsis: numseq.csh infile outfile gcfile start finish nuccase startno group gpl which strands kind numbers nucs peptides frames form

#Convert arguments to variables
set INFILE  = $1
set OUTFILE = $2
set GCFILE  = $3
set START   = $4
set FINISH  = $5
set NUCCASE = $6
set STARTNO = $7
set GROUP   = $8
set GPL     = $9
set WHICH   = $10
set STRANDS = $11
set KIND    = $12
set NUMBERS = $13
set NUCS    = $14
set PEPTIDES = $15
set FRAMES  = $16
set FORM    = $17

set PID = $$ #process id

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

# Abort if INFILE does not exist or is of zero length
if ((-e $INFILE) && (! -z $INFILE)) then

   # Correct for bug in GDE2.2, which writes 'Circular' in wrong columns
   # in GenBank LOCUS line. This is done by determining whether seq is
   #   linear or circular, converting to Pearson format, and appending the
   # '1' if linear or '2' if circular.
   head -1 $INFILE |grep 'Circular' > circ.$PID
   if (-z circ.$PID) then
      set TOPOLOGY = 1
   else 
      set TOPOLOGY = 2
   endif
   
   readseq -i1 -ffasta -o=$INFILE.$PID $INFILE 
   #readseq  -i1 -ffasta $INFILE> $INFILE.$PID  
   

   echo $TOPOLOGY >> $INFILE.$PID


   # Determine the size of seq. in infile and make sure that
   # start and finish values aren't greater than seq. length.
   tail +2 $INFILE.$PID > TEMP.$PID
   @ SEQLENGTH = (`wc -c < TEMP.$PID` - `wc -l < TEMP.$PID`) - 1
   if ($START > $SEQLENGTH) set START = $SEQLENGTH
   if ($FINISH > $SEQLENGTH) set FINISH = $SEQLENGTH
   
   #initial filenames 
   echo $INFILE.$PID > $PID.COMFILE #input filename
   echo b            >> $PID.COMFILE #Pearson format
   echo $OUTFILE     >> $PID.COMFILE #outfile
   echo 3            >> $PID.COMFILE #Genetic Code file
   echo $GCFILE >> $PID.COMFILE 
   

   # Set parameters
   echo 4     >> $PID.COMFILE #Choose parameter menu
   
   echo 1     >> $PID.COMFILE #choose start
   echo $START >> $PID.COMFILE 
   
   echo 2     >> $PID.COMFILE #choose finish
   echo $FINISH >> $PID.COMFILE 

   echo 3     >> $PID.COMFILE #choose nuccase 
   echo $NUCCASE >> $PID.COMFILE 
   
   # STARTNO = 0 means do not set STARTNO (default = START)
   if ($STARTNO != 0) then
     echo 4     >> $PID.COMFILE #choose startno 
     echo $STARTNO >> $PID.COMFILE 
   endif

   #if PEPTIDES = Yes, round GROUP to nearest multiple of three.
   if ($PEPTIDES == Y) then 
      @ remainder = $GROUP % 3 
      @ GROUP = $GROUP - $remainder
   endif
      
   echo 5      >> $PID.COMFILE #choose group 
   echo $GROUP >> $PID.COMFILE 
   
   echo 6     >> $PID.COMFILE #choose gpl 
   echo $GPL >> $PID.COMFILE 
   
   echo 7     >> $PID.COMFILE #choose which
   echo $WHICH >> $PID.COMFILE 
   
   echo 8     >> $PID.COMFILE #choose strands
   echo $STRANDS >> $PID.COMFILE 
   
   echo 9     >> $PID.COMFILE #choose kind
   echo $KIND >> $PID.COMFILE 
   
   echo 10    >> $PID.COMFILE #choose numbers 
   echo $NUMBERS >> $PID.COMFILE 
   
   echo 11    >> $PID.COMFILE #choose nucs    
   echo $NUCS >> $PID.COMFILE 
   
   echo 12    >> $PID.COMFILE #choose peptides    
   echo $PEPTIDES >> $PID.COMFILE 
   
   echo 13    >> $PID.COMFILE #choose frames    
   echo $FRAMES >> $PID.COMFILE 
   
   echo 14    >> $PID.COMFILE #choose form    
   echo $FORM >> $PID.COMFILE 
   
   echo 0     >> $PID.COMFILE #exit parameter menu

   echo 6     >> $PID.COMFILE #Print numbered sequence 
   echo ""    >> $PID.COMFILE #dummy prompt line
   echo 0     >> $PID.COMFILE #exit program
   
   #run NUMSEQ
   numseq < $PID.COMFILE
   
   $RM_CMD circ.$PID TEMP.$PID $PID.COMFILE
  endif
