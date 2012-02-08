#!/bin/csh
# mrtrans.csh                  Feb 20  2002
# This is a front end for mrtrans. It makes sure that the names of
# the sequences in PROTFILE and DNAFILE are the same, and re-orders
# the sequences in DNAFILE, if necessary, to be in the same order
# as in PROTFILE.
# This script assumes that sequence names in PROTFILE are IDENTICAL to
# the corresponding names in DNAFILE.

set PROTFILE = $1
set DNAFILE = $2
set JOBID = $$

#Make a list of sequence names from each file.
grep '>' $PROTFILE | cut -c2-80 | cut -f1 -d" " > $JOBID.pro.nam
grep '>' $DNAFILE  | cut -c2-80 | cut -f1 -d" " > $JOBID.dna.nam

echo Making sure all names in protein file have a counterpart in dna file...

sort < $JOBID.pro.nam > $JOBID.pro.nam.sorted
sort < $JOBID.dna.nam > $JOBID.dna.nam.sorted
comm -23 $JOBID.pro.nam.sorted $JOBID.dna.nam.sorted > $JOBID.missing

if (-z $JOBID.missing) then
  # All protein sequences have a DNA counterpart. 
  # Re-order the DNA file into the same order as the protein file.
  # A $JOBID.dna.num is the same ase $JOBID.dna.nam, except that
  # the former begins with a line number for each name in the 
  # latter. For each protein sequence, the appropriate line
  # number is chosen using grep, and readseq pulls out that
  # sequence from $DNAFILE.
echo Re-ordering the DNA file to correspond to order of sequences
echo in protein file.
  cat -n $JOBID.dna.nam |tr -s " " " " > $JOBID.dna.num
  foreach name (`cat $JOBID.pro.nam`)
     set LINENUM = `grep -w "$name" $JOBID.dna.num | cut -f1 -d"	" |tr -d ' ' `
     readseq -pipe -i$LINENUM -fPearson $DNAFILE >> $DNAFILE.reordered
  end  
else
  echo 'The following sequences are present in the protein file'
  echo 'but missing from the DNA file:'
  cat $JOBID.missing
endif

/usr/bin/rm $JOBID.*

 
