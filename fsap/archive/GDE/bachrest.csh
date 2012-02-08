#!/bin/csh

#VERSION  3/29/06

# Run Bachrest as a command
#Synopsis: bachrest.csh infile refile outfile source proto prot3 blunt prot5 symm \
#            minsite maxsite fragleast fragmost fragprint
 
#Convert arguments to variables
set INFILE  = $1
set REFILE  = $2
set OUTFILE = $3
set SOURCE  = $4
set PROTO   = $5
set PROT3   = $6
set BLUNT   = $7
set PROT5   = $8
set SYMM    = $9
set MINSITE = $10
set MAXSITE = $11
set FRAGLEAST = $12
set FRAGMOST = $13
set FRAGPRINT = $14

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
set PID = $$ #process id


if ((-e $INFILE) && (! -z $INFILE)) then
   # Correct for bug in GDE2.2, which writes 'Circular' in wrong columns
   # in GenBank LOCUS line. This is done by determining whether seq is
   # linear or circular, converting to Pearson format, and appending the
   # '1' if linear or '2' if circular. We also have to work around a bug
   # in readseq that causes it to write 'Circular' on the DEFINITION line
   # ARRRGH!
   head -1 $INFILE |egrep -e '[Cc]ircular' > circ.$PID
   if (-z circ.$PID) then
      set TOPOLOGY = 1
   else 
      set TOPOLOGY = 2
   endif

   #readseq -f1 -o$INFILE.$PID $INFILE
   # Consistent with Java Readseq and original C readseq.
   readseq   -i1 -ffasta -o=$INFILE.$PID $INFILE 

   echo $TOPOLOGY >> $INFILE.$PID

   if ($FRAGLEAST > $FRAGMOST) then 
     set FRAGLEAST = $FRAGMOST
     endif
   if ($FRAGPRINT > $FRAGMOST) then
     set FRAGPRINT = $FRAGMOST
     endif
   if ($FRAGPRINT < 1) then
     set FRAGPRINT = 1
     endif
   if ($MAXSITE > 23) then
     set MAXSITE = 23
     endif
   if ($MINSITE > $MAXSITE) then
     set MINSITE = $MAXSITE  
     endif

   #initialize 
   echo $INFILE.$PID > $PID.COMFILE   #input filename
   echo b            >> $PID.COMFILE  #Intelligenetics format
   echo $REFILE      >> $PID.COMFILE  #Restriction site file
   echo $OUTFILE     >> $PID.COMFILE  #output file

   # Set parameters
   echo 4         >> $PID.COMFILE  #Choose parameter menu
   
   echo 1         >> $PID.COMFILE  #choose source 
   echo $SOURCE   >> $PID.COMFILE  
   
   echo 2         >> $PID.COMFILE  #choose prototype/all
   echo $PROTO    >> $PID.COMFILE 
 
   echo 3         >> $PID.COMFILE  #set 3' protruding ends
   echo $PROT3    >> $PID.COMFILE 
    
   echo 4         >> $PID.COMFILE  #set blunt ends
   echo $BLUNT    >> $PID.COMFILE  
   
   echo 5         >> $PID.COMFILE  #choose 5' protruding ends
   echo $PROT5    >> $PID.COMFILE  
 
   echo 6         >> $PID.COMFILE  #choose symmetric, asymmetric or both
   echo $SYMM     >> $PID.COMFILE   
    
   echo 7         >> $PID.COMFILE  #choose minsite
   echo $MINSITE  >> $PID.COMFILE  
   
   echo 8         >> $PID.COMFILE  #choose maxsite
   echo $MAXSITE  >> $PID.COMFILE  
 
   echo 9         >> $PID.COMFILE  #set FRAGLEAST
   echo $FRAGLEAST >> $PID.COMFILE   
    
   echo 10         >> $PID.COMFILE  #set FRAGMOST
   echo $FRAGMOST  >> $PID.COMFILE  
   
   echo 11          >> $PID.COMFILE  #set FRAGPRINT
   echo $FRAGPRINT  >> $PID.COMFILE  
      
   echo 0         >> $PID.COMFILE  #exit parameter menu

   echo 5         >> $PID.COMFILE  #Search for sites & write output to file

   echo 0         >> $PID.COMFILE  #exit program

   nice bachrest < $PID.COMFILE
   $RM_CMD circ.$PID $PID.COMFILE

endif
