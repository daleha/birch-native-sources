#! /usr/local/bin/python

# update May 13, 2002

# Create .fil files for use by FASTA
# Each division of GenBank is broken up into one or more sequence files.
# Where only one file exists for a division, the name is the sequence
# files is something like gbuna.wrp.
# Where the division is broken up into several files, the names contain
# numbers eg. gbest1.wrp, gbest2.wrp, gbest3.wrp, gbest4.wrp ...
# Input file, filnum is in the form
#
# div num
#
# where div is the 3 letter GenBank division code, and num
# is the number of files in the current release, for that
# division
#
# eg. est 169
#     gss 51
#
# 
# This script also sets the permissions so that the files are
# world-readable.

import sys
import os

INFILE = open('filnum','r')

# Get the path for the GenBank directory, as specified
# in the environment variable $GB.
GBPATH = os.environ['GB'] 

# Set file extension for sequence files
EXTENSION = 'wrp'
if len(sys.argv) > 2:
   if sys.argv[1] == '-e':
      EXTENSION = sys.argv[2]


# In addition to writing a .fil file for each division,
# all filenames are echoed to genbank.fil. genbank.fil is
# used when you want to search all divisions of genbank
# in a single search.
GBFILE = open('genbank.fil','w')
print "Creating genbank.fil"
GBFILE.write('<' + GBPATH + '\n') 

LINE = INFILE.readline()

# For each GenBank division (div) listed in INFILE, create
# a file called gbdiv.fil. 
while LINE != '':
      TOKEN = LINE.split()
      OUTFILENAME =  'gb' + TOKEN[0] + '.fil'
      print "Creating ", OUTFILENAME
      OUTFILE = open(OUTFILENAME,'w')
      OUTFILE.write('<' + GBPATH + '\n') 

      #For each .wrp file in a GenBank division, write the
      # name of that file to gbdiv.fil, and also to genbank.fil
      # For the special case in which the entire division is in
      # a single file, the filename does not contain a number.
      NUM = int(TOKEN[1])
      for i in range(1,NUM+1):
          if NUM == 1:
             NAME = 'gb' + TOKEN[0] + '.' + EXTENSION + ' 0' 
          else:
             NAME = 'gb' + TOKEN[0] + str(i) + '.' + EXTENSION + ' 0'
          OUTFILE.write(NAME + '\n') 
          GBFILE.write(NAME + '\n') 
      OUTFILE.close()
      os.chmod(OUTFILENAME,0644)
      LINE = INFILE.readline()

       
INFILE.close()
GBFILE.close()
os.chmod('genbank.fil',0644)

