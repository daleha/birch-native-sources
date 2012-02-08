#!/bin/csh
# UDS - Update a database subset - Version 8/26/95

# The union of olddatabase and newdatabase is written to $NEW.$DBEXT.UPDATED
# Entries present in olddatabase are replaced by their equivalents from
# newdatabase.      

if ($2 == "") then
   echo "Usage: uds olddatabase newdatabase"
   exit
   endif

set OLD = $1:r
set NEW = $2:r
set DBEXT = $2:e

# Make sure databases are of the same type
if ($DBEXT != $1:e) then
   echo ">>> Databases not of same type"
   exit
   endif

if ($DBEXT == gen | $DBEXT == pir) then

   # set database option for use with getloc
   switch ($DBEXT)
     case gen
       set DBFLAG = -g
       breaksw
     case pir
        set DBFLAG = -p
        breaksw
     endsw

   # carry out operations in a temporary working directory
   mkdir UDS.tmp
   cd UDS.tmp
   splitdb $DBFLAG ../$1 $OLD.ano $OLD.wrp $OLD.ind
   splitdb $DBFLAG ../$2 $NEW.ano $NEW.wrp $NEW.ind


  # Find accession numbers that are in the old database subset but
  # not in the new. First, lists of accession numbers from the old and
  # new indexes are compared, and any accession numbers unique to the old
  # index are written to notfound.tmp. Next, each accession number from 
  # notfound.tmp is searched for among the ACCESSION or #Accession lines
  # in the new database subset. If they are found, then it is assumed that
  # the old entry was merged into another entry. Any accession numbers not 
  # found at this point are written to $OLD.unique
  tr -s ' ' ' ' < $OLD.ind |cut -f2 -d" " | sort | uniq > $OLD.acc
  tr -s ' ' ' ' < $NEW.ind |cut -f2 -d" " | sort | uniq > $NEW.acc
  fgrep -v -f $NEW.acc $OLD.acc > notfound.tmp
  grep  -i 'ACCESSION' $NEW.ano > $NEW.acl
  foreach ACNO (`cat notfound.tmp`)
     egrep -e $ACNO $NEW.acl > has_been_merged.tmp 
     if (-z has_been_merged.tmp) then
      echo $ACNO >> $OLD.unique
      endif
    end
 
  if (! -z $OLD.unique) then

     # Extract these unique entries as individual files in UDS.tmp
     getloc $DBFLAG -c -f $OLD.unique $OLD.ano $OLD.wrp $OLD.ind  
    
     # Extract all of the entries from $NEW as individual files in UDS.tmp
     getloc $DBFLAG -f $NEW.ind $NEW.ano $NEW.wrp $NEW.ind

     #Rename files so that their names are LOCUS names, not ACCESSION 
     # numbers. This will ensure that files get written in alphabetical
     # order by LOCUS name.
     foreach file (*.$DBEXT)
       set LNAME = `grep LOCUS $file |tr -s  ' ' ' ' |cut -f2 -d" "`
       if ($file != $LNAME.$DBEXT) mv $file $LNAME.$DBEXT
       end
     # Write all of the entries collected in UDS.tmp into one file in
     # the original working directory
    echo Writing updated database to $NEW.$DBEXT.UPDATED
     if (-e ../$NEW.$DBEXT.UPDATED) /usr/bin/rm ../$NEW.$DBEXT.UPDATED
     cat *.$DBEXT >> ../$NEW.$DBEXT.UPDATED

   else
     echo "All accession numbers in $1 are also found in $2"
     echo "No action taken."
   endif

  # Cleanup - remove temporary directory
  cd ..
/usr/bin/rm -r UDS.tmp

else
  echo ">>> Unknown database file extension ."$DBEXT
endif
