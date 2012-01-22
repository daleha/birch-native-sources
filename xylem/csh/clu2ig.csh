#!/bin/csh
# CLU2IG - Convert interleaved Clustal V .aln output to sequential
# .ig format suitable for use with MASE
# Brian Fristensky        August 3, 1992

# Usage: cul2ig clustalfile > masefile

# Clean up the file by removing title lines, '*' and blank lines
grep -v "CLUSTAL V" $1 | grep -v '*' | egrep -v -e ^$ > TEMP.CLU2MASE

# Create a namefile from the first 16 columns of the cleaned-up file.
cut -c1-16 TEMP.CLU2MASE > RAWNAMES.CLU2MASE

# Find out how many unique names there are, and store in the environment
# variable numnames. Write the first $NUMNAMES lines of RAWNAMES.CLU2MASE
# to UNIQUE.CLU2MASE
set numnames = `cat RAWNAMES.CLU2MASE | sort | uniq | wc -l`
head -$numnames  RAWNAMES.CLU2MASE> UNIQUE.CLU2MASE

# For each unique name in the namefile, grep out the corresponding lines
# and write to the output file

foreach name (`cat UNIQUE.CLU2MASE`)
    #create a dummy comment
    echo ";" 

    #write the name
    echo $name 

    #copy the sequence
    grep -w $name TEMP.CLU2MASE |cut -c17-80 
  end

/usr/bin/rm *.CLU2MASE
