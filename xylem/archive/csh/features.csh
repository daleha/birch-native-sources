#!/bin/csh

# FEATURES - a c-shell script (version  20 Nov 01)
# FEATURES is a front end for GETOB. FEATURES provides an interactive
# menu interface, as well as automating the process of resolving
# nested database expressions

# This script automates, as much as is possible, two tasks. First, getob 
# expects all of its input to come from files, as described in getob.doc.
# The parameter menu lets the user input some of the data directly at the 
# keyboard, for cases in which only one data element (eg. accession number)
# is needed. Otherwise, FEATURES prompts the user for a filename.
# The second thing that FEATURES does is to run getob and check to see
# whether there are unresolved expressions in the output. If there are,
# getob is called as many times as is necessary to resolve all expressions.

#####################################################################
# set parameter defaults
####################################################################
set mode = "interactive"
set NICELEVEL = 10
set jobid = $$
set featype = ""
set features = ""
set entrytype = "n"
set entries = ""
set dbtype = "g"
set database = ""
set where = "a"
set outfile  = ""
set resolve = ""
set acno = "OFF"
set needtosplit = "OFF"
set options   = ""

# Platform-specific syntax is chosen based on XYLEM_PLATFORM.
# (default = sun)
if !(${?XYLEM_PLATFORM}) set XYLEM_PLATFORM = sun #Sun4, Sparcstations

switch ($XYLEM_PLATFORM)
  case "sun":
    set RM = "rm -f"
    breaksw
  case "HP":
    set RM = "rm -f"
    breaksw
  default:
    set RM = "rm -f"
    breaksw
endsw

#####################################################################
# Read options from the command line and set parameters,
# or prompt user for parameters at terminal.
####################################################################
  set numargs = $#argv

if ($numargs != 0) then
  #---------------------------------------------------
  #parameters given in command line
  set mode = command
  set files = ""

  #One argument form: a FEATURES expression is given as the input, returning ONLY
  #the resultant sequence to the standard output.
  if ($numargs == 1) then
     # -h : print usage message
     if ($1 == "-h") then
        echo "Usage:    features"
        echo "          features expression"
        echo "          features [-f featurekey | -F keyfile]"
        echo "                   [-n name     |-a accession    | -e expression |"
        echo "                    -N namefile |-A accfile      | -E expfile]"
        echo "                   [-u dataset  | -U dataset     | -g ]"
        exit
      else
	set resolve = "r"
	set features = ""
	set entries = $1
        #Make sure the expression begins with '@'
        if (`echo $entries |cut -c1-1` == @) then
        else
          set entries = "@"$entries
        endif

	echo "OBJECTS" > FEA.$jobid.inf #make dummy infile
	echo "SITES" >> FEA.$jobid.inf  #for getob
	echo $entries > FEA.$jobid.nam
        #accession number will be used for output filenames
	set base = `cut -c2-80 FEA.$jobid.nam | cut -f1 -d':'` 
	set acno = "ON"
     endif
  else
  #Multiple argument form: feature keys, entries, expressions, and databases are
  #specified on the command line. Message, sequence and expression files are written.
      set index = 1
      while ($index < $numargs)
	switch ($argv[$index])
	  case "-f":
	    @ index++
	    set features = $argv[$index]
	    set resolve = ""
	    echo "OBJECTS" > FEA.$jobid.inf
	    echo $features >> FEA.$jobid.inf 
	    echo "SITES" >> FEA.$jobid.inf
	    breaksw
	  case "-F":
	    @ index++
	    set features = $argv[$index]
	    set resolve = ""
	    set features = $argv[$index]
	    echo "OBJECTS" > FEA.$jobid.inf
	    cat $features >> FEA.$jobid.inf 
	    echo "SITES" >> FEA.$jobid.inf
	    breaksw
	  case "-n":
	    @ index++
	    set entries = $argv[$index]
	    set base = $entries  #entry name used for output file names
	    echo $entries > FEA.$jobid.nam
	    set acno = "OFF"
	    breaksw
	  case "-a":
	    @ index++
	    set entries = $argv[$index]
	    set base = $entries #accession number used for output file names
	    echo $entries > FEA.$jobid.nam
	    set acno = "ON"
	    breaksw
	  case "-e":
	    set resolve = "r"
	    set features = ""
	    @ index++
	    set entries = $argv[$index]

            #Make sure the expression begins with '@'
            if (`echo $entries |cut -c1-1` == @) then
            else
              set entries = "@"$entries
            endif

	    echo "OBJECTS" > FEA.$jobid.inf #make dummy infile
	    echo "SITES" >> FEA.$jobid.inf  #for getob
	    echo $entries > FEA.$jobid.nam
            #accession number will be used for output filenames
	    set base = `cut -c2-80 FEA.$jobid.nam |cut -f1 -d':'` 
	    set acno = "ON"
	    breaksw
	  case "-N":
	    @ index++
	    set entries = $argv[$index]
	    set base = $entries:t
	    set base = $base:r    #raw filename used for output filenames
	    cat $entries > FEA.$jobid.nam
	    set acno = "OFF"
	    breaksw
	  case "-A":
	    @ index++
	    set entries = $argv[$index]
	    set base = $entries:t
	    set base = $base:r    #raw filename used for output filenames
	    cat $entries > FEA.$jobid.nam
	    set acno = "ON"
	    breaksw
	  case "-E":
	    set resolve = "r"
	    set features = ""
	    @ index++
	    set entries = $argv[$index]
	    set base = $entries:t
	    set base = $base:r    #raw filename used for output filenames
	    cat $entries > FEA.$jobid.nam
	    echo "OBJECTS" > FEA.$jobid.inf #make dummy infile
	    echo "SITES" >> FEA.$jobid.inf  #for getob
	    set acno = "ON"
	    breaksw
	  case "-u":
	    set dbtype = u
	    @ index++
	    set database = $argv[$index]
	    set dbname = $database:r
	    set extension = $database:e
	    if  ($extension == gen) then 
		 #If user-defined database is one or more
		 #GenBank entries, it will have to be 
		 #split by splitdb for use by programs.
		set needtosplit = ON
	    else 
		set needtosplit = OFF
	    endif
	    set anofile = $dbname.ano
	    set seqfile = $dbname.wrp
	    set indfile = $dbname.ind
	    breaksw
	  case "-U":
	    set dbtype = U
	    @ index++
	    set database = $argv[$index]
	    set dbname = $database:r
            set base = $dbname
	    set extension = $database:e
	    if  ($extension == gen) then 
		 #If user-defined database is one or more
		 #GenBank entries, it will have to be 
		 #split by splitdb for use by programs.
		set needtosplit = ON
	    else 
		set needtosplit = OFF
	    endif
	    set anofile = $dbname.ano
	    set seqfile = $dbname.wrp
	    set indfile = $dbname.ind
	    breaksw
	  endsw
	  @ index++
      end # while

      #Make sure all necessary data has been entered
      if ($resolve != r && $features == "") then
           echo '>>> Features must be specified' 
           exit
         endif
        # csh subtlety: $entries must be quoted to allow for blanks in
        # the string
        if ("$entries" == "" && $dbtype != "U") then 
           echo '>>> Entries must be specified'
           exit
        endif

  endif
else 
  #---------------------------------------------------------------
  # Interactive parameter input
  set mode = interactive
  set complete = 0
  while ($complete == 0)
      #Display current parameter settings
      echo '___________________________________________________________________'
      echo '                     FEATURES - Version   20 Nov 01                '
      echo '    Please cite: Fristensky (1993) Nucl. Acids Res. 21:5997-6003'
      echo '___________________________________________________________________'
      echo 'Features: ' $features
      echo 'Entries:  ' $entries
      echo 'Dataset:  ' $database
      echo '___________________________________________________________________'
      echo '   Parameter              Description                      Value'
      echo '-------------------------------------------------------------------'
      echo '1).................... FEATURES TO EXTRACT ....................> '$featype
      echo '  f:Type a feature at the keyboard '
      echo '  F:Read a list of features from a file'
      echo '2)....................ENTRIES TO BE PROCESSED (choose one).....> '$entrytype
      echo '  Keyboard input - n:name     a:accession #     e:expression'
      echo '  File input     - N:name(s)  A:accession #(s)  E:expression(s)' 
      echo '3)....................WHERE TO GET IT .........................> '$dbtype
      echo '  u:Genbank dataset   g:complete GenBank database'
      echo '  U: same as u, but all entries'
      echo '4)....................WHERE TO SEND IT ........................> '$where
      echo '  s:Each feature to a separate file  a:All output to same file'
      echo '   ---------------------------------------------------------------'
      echo '   Type number of your choice or 0 to continue:'

      #Prompt for parameter to change
      set paramnum = $<
      switch ($paramnum)
        case 0:
          set complete = 1
          breaksw
        case 1:
          if ($resolve == r) then
             echo '>>> No features needed to resolve expressions.'
          else
            echo 'Enter one of {fF}:'
            set temp = $<
            if ($temp == f | $temp == F) then 
               set featype = $temp
               switch ($featype)
                 case f:
                   set resolve = ""
                   echo 'Enter feature:'
                   set features = $<
                   echo "OBJECTS" > FEA.$jobid.inf
                   echo $features >> FEA.$jobid.inf 
                   echo "SITES" >> FEA.$jobid.inf               
                   breaksw
                 case F:
                   set resolve = ""
                   echo 'Enter filename:'
                   set features = $<
                   cat $features > FEA.$jobid.inf
                   breaksw
                 endsw
            else
              echo '>>> Invalid choice'
            endif # f,F
          endif # resolve == r
          breaksw # 1
        case 2:
          echo 'Enter one of {naeNAE}:'
          set temp = $<
          if ($temp == n | $temp == a | $temp == e |\
              $temp == N | $temp == A | $temp == E) then 
             set entrytype = $temp
             switch ($entrytype)
               case n:
                 echo 'Enter name:'
                 set entries = $<
                 set base = $entries  #entry name used for output file names
                 echo $entries > FEA.$jobid.nam
                 set acno = "OFF"
                 breaksw
               case a:
                 echo 'Enter accession number:'
                 set entries = $<
                 set base = $entries #accession number used for
                                     # output file names
                 echo $entries > FEA.$jobid.nam
                 set acno = "ON"
                 breaksw
               case e:
                 set resolve = "r"
                 set features = ""
                 echo 'Enter expression:'
                 set entries = $<

                 #Make sure the expression begins with '@'
                 if (`echo $entries |cut -c1-1` == @) then
                 else
                   set entries = "@""$entries"
                 endif

                 echo "OBJECTS" > FEA.$jobid.inf #make dummy infile
                 echo "SITES" >> FEA.$jobid.inf  #for getob
                 echo $entries > FEA.$jobid.nam
                 #accession number will be used for output filenames
                 set base = `cut -c2-80 FEA.$jobid.nam |cut -f1 -d':'` 
                 set acno = "ON"
                 breaksw
               case N:
                 echo 'Enter filename:'
                 set entries = $<
                 set base = $entries:t
                 set base = $base:r    #raw filename used for output filenames
                 cat $entries > FEA.$jobid.nam
                 set acno = "OFF"
                 breaksw
               case A:
                 echo 'Enter filename:'
                 set entries = $<
                 set base = $entries:t
                 set base = $base:r    #raw filename used for output filenames
                 cat $entries > FEA.$jobid.nam
                 set acno = "ON"
                 breaksw
               case E:
                 set resolve = "r"
                 set features = ""
                 echo 'Enter filename:'
                 set entries = $<
                 set base = $entries:t
                 set base = $base:r    #raw filename used for output filenames
                 cat $entries > FEA.$jobid.nam
                 echo "OBJECTS" > FEA.$jobid.inf #make dummy infile
                 echo "SITES" >> FEA.$jobid.inf  #for getob
                 set acno = "ON"
                 breaksw
               endsw
          else 
            echo '>>> Invalid choice'
          endif #naeNAE
          breaksw
        case 3:
          echo 'Enter one of {uUg}:'
          set temp = $<
          switch ($temp)
             case u:
               set dbtype = u
               echo 'Enter filename:'
               set database = $<
               set dbname = $database:r
               set extension = $database:e
               if  ($extension == gen) then
                   #If user-defined database is one or more
                   #GenBank entries, it will have to be
                   #split by splitdb for use by programs.
                   set needtosplit = ON
               else
                   set needtosplit = OFF
               endif
               set anofile = $dbname.ano
               set seqfile = $dbname.wrp
               set indfile = $dbname.ind
               breaksw
             case U:
               set dbtype = U
               echo 'Enter filename:'
               set database = $<
               set dbname = $database:r
               set base = $dbname
               set extension = $database:e
               if  ($extension == gen) then
                   #If user-defined database is one or more
                   #GenBank entries, it will have to be
                   #split by splitdb for use by programs.
                   set needtosplit = ON
               else
                   set needtosplit = OFF
               endif
               set anofile = $dbname.ano
               set seqfile = $dbname.wrp
               set indfile = $dbname.ind
               breaksw
             case g:
               set dbtype = g
               breaksw
             default:
               echo '>>> Invalid choice'
               breaksw
             endsw
           breaksw
        case 4:
          echo 'Enter one of {sa}:'
          set temp = $<
          if ($temp == s | $temp == a) then 
             set where = $temp
          else
            echo '>>> Invalid choice'
          endif
          breaksw
        endsw

      #Make sure all necessary data has been entered
      if ($complete == 1) then
         if ($resolve != r && $features == "") then
           echo '>>> Features must be specified (option 1)' 
           set complete = 0
         endif
        # csh subtlety: $entries must be quoted to allow for blanks in
        # the string
        if ("$entries" == "" && $dbtype != "U") then 
           echo '>>> Entries must be specified (option 2)'
           set complete = 0
        endif
      endif
    end #while 
  endif # argv 


  #set options
  if ($resolve == r) then 
     set options = "$options -r"
  endif
  if ($acno == ON) then
     set options = "$options -c"
  endif 
  if ($where == s) then 
     set options = "$options -f"
  endif

#####################################################################
#  For a given database, extract the specified features.
####################################################################

  if ($mode == interactive) then
    set termout = /dev/tty
    echo 'Messages will be written to '$base.msg
    if ($where == a) then
        echo 'Final sequence output will be written to ' $base.out
     else
        echo 'Output for each entry will be written to individual files.'
        echo 'Filenames will be in the form:    name.obj'
     endif
     if ($resolve != r) then
        echo 'Expressions will be written to '$base.exp
     endif
  else
    set termout = /dev/null
  endif

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # The first pass will usually resolve all feature expressions. 
  #Retrieve GenBank entries, if necessary
  if ($dbtype == g) then
     if ($resolve == r) then
        #Unresolved expressions will begin with '@' in first column.
        #Find lines with unresolved expressions and extract the
        #accession numbers for use by fetch.
        egrep -e ^@ FEA.$jobid.nam > UNRESOLVED.$jobid.fea
        cut -c2-80 UNRESOLVED.$jobid.fea |cut -f1 -d':' > UNRESOLVED.$jobid.nam
        set fetchfile = UNRESOLVED.$jobid.nam
     else
        set fetchfile = FEA.$jobid.nam
     endif
     set anofile = $base.ano
     set seqfile = $base.wrp
     set indfile = $base.ind
     fetch $fetchfile FEA.$jobid.gen > $termout
     if (-e FEA.$jobid.gen) then
        nice -$NICELEVEL splitdb FEA.$jobid.gen $anofile $seqfile $indfile
        $RM FEA.$jobid.gen
        else
          exit
        endif
     endif

  #Extract features from the entries
if ($mode == interactive)  echo 'Extracting features...'
  if ($needtosplit == ON) then
     nice -$NICELEVEL splitdb $database $anofile $seqfile $indfile
     endif
  if ($dbtype == 'U') cp $indfile FEA.$jobid.nam
  if ($where == a) then
     if ($numargs == 1) then
        nice -$NICELEVEL getob $options FEA.$jobid.inf FEA.$jobid.nam $anofile $seqfile\
               $indfile /dev/null FEA.$jobid.out /dev/null  > $base.$jobid.term

     else
       nice -$NICELEVEL getob $options FEA.$jobid.inf FEA.$jobid.nam $anofile $seqfile $indfile\
               $base.msg FEA.$jobid.out $base.exp  > $base.$jobid.term
     endif
     else
       nice -$NICELEVEL getob $options FEA.$jobid.inf FEA.$jobid.nam $anofile $seqfile $indfile\
               $base.msg  $base.exp  > $base.$jobid.term
  endif

  # Clean up a bit. Database files could be huge!
  if ($dbtype == g | $needtosplit == ON) then 
     $RM $anofile $seqfile $indfile
     endif

  #Pull out all lines containing indirect references
  if ($where == a) then
     egrep -e ^@ FEA.$jobid.out > UNRESOLVED.$jobid.fea
  else
     egrep -e ^@  *.obj > UNRESOLVED.$jobid.fea
  endif


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
  # This part resolves any unresolved expressions that can be resolved.
  # Only in very rare cases will this loop execute more than once.

  while !(-z UNRESOLVED.$jobid.fea) #is not empty
    if ($where == s) then
       echo ">>> Can't resolve references sent to separate files."
       echo '>>> Re-run FEATURES separately for each unresolved reference'
       echo ">>> See $base.msg to find unresolved references."
       exit 
    else
      if ($mode == interactive) echo 'Resolving Feature expression(s)...'
      #extract accession numbers to be retrieved
      cut -c2-80 UNRESOLVED.$jobid.fea | cut -f1 -d':' > FEA.$jobid.nam
 
      #retrieve the sequences into a new file, and create
      #a database subset to be used by getob. If $GB is set, then 
      #entries will be fetched from the online GenBank database.
      #Otherwise, features will attempt to get entries from the
      #user-generated database subset.
      if (${?GB} | ${?XYLEM_RHOST}) then
         set fetchoptions = ""
      else 
         set fetchoptions = "-G $database"
      endif
      fetch $fetchoptions FEA.$jobid.nam FEA.$jobid.gen >$termout
      splitdb FEA.$jobid.gen FEA.$jobid.ano FEA.$jobid.wrp FEA.$jobid.ind
                 
      #run getob again to resolve indirect references
      echo 'Extracting features...' > $termout
      mv FEA.$jobid.out UNRESOLVED.$jobid.out
      if ($numargs == 1) then
        nice -$NICELEVEL getob -r -c FEA.$jobid.inf UNRESOLVED.$jobid.out FEA.$jobid.ano\
                  FEA.$jobid.wrp FEA.$jobid.ind /dev/null FEA.$jobid.out  >>\
                    $base.$jobid.term

      else
        nice -$NICELEVEL getob -r -c FEA.$jobid.inf UNRESOLVED.$jobid.out FEA.$jobid.ano\
                    FEA.$jobid.wrp FEA.$jobid.ind FEA.$jobid.msg FEA.$jobid.out\
                    >> $base.$jobid.term
        cat FEA.$jobid.msg >> $base.msg 
      endif

      #Pull out all lines containing indirect references 
      #ie. lines beginning with '@' 
      egrep -e ^@ FEA.$jobid.out > UNRESOLVED.$jobid.fea
     endif
     end

 if ($where == a) then
    #One argument form, command mode
    #Sequence output only goes to standard output.
    if ($numargs == 1) then
      cat FEA.$jobid.out
    #Multi-argument form or interactive mode
    else
      mv FEA.$jobid.out $base.out
    endif
  endif

#####################################################################
#  Clean up.
####################################################################
$RM FEA.$jobid.* UNRESOLVED.$jobid.* *.$jobid.term
