#!/bin/csh

# FINDKEY - a c-shell script (version 15 Feb 02)
# FINDKEY is a front end for IDENTIFY. It provides an interactive
# menu interface, as well as overseeing the process of locating database
# entries by keywords.
#
# Revisions:

#    15 Feb 02  Run identify using nice to minimize impact of
#               runaway jobs. More careful checking is needed 
#               in the future.
#    15 Feb 02  Support for VecBase removed
#    13 Mar 97  Added GSS and HTG divisions
#
#####################################################################
# set default parameters
####################################################################
# It is assumed that the following environment variables have been set
# by .cshrc or .login:
#   GB  - GenBank directory
#   PIR - PIR/NBRF directory
#   VEC - VecBase directory

# Running findkey on a remote host
# If your databases are on a different host, you can still run findkey
# locally, and it will use rsh to run a findkey job on the remote host.
# To facilitate remote execution, you must set the environment variable
# XYLEM_RHOST to the name of the remote host, and XYLEM_USERID to your
# userid on that remote host (should be done in .cshrc)

set keyword = ""
set keyfile = ""
set dbfile = ""
set namefile = ""
set findfile = ""
set wheretolook = "p"
set mode      = interactive
set options   = "" #string to hold command line options for getloc
set jobid = $$
set runlocal = ""  # runlocal = y prevents remote execution as
                   # specified by XYLEM_RHOST.
set NICE = 'nice +8' 

# Platform-specific syntax is chosen based on XYLEM_PLATFORM.
# (default = sun)
if !(${?XYLEM_PLATFORM}) set XYLEM_PLATFORM = sun #Sun4, Sparcstations

switch ($XYLEM_PLATFORM)
  case "sun":
    set RM = "/usr/bin/rm -f"
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
  set mode  = command
  set index = 1
  while ($index <= $numargs) 
   set a = $argv[$index]
   switch ($a)
     case "-p":
       set wheretolook = p
       breaksw
     case "-G":
       set wheretolook = G
       @ index++
       set dbfile = $argv[$index] 
       breaksw
     case "-P":
       set wheretolook = P
       @ index++
       set dbfile = $argv[$index] 
       breaksw
     case "-b":
       set wheretolook = b
       breaksw
     case "-m":
       set wheretolook = m
       breaksw
     case "-g":
       set wheretolook = g
       breaksw
     case "-r":
       set wheretolook = r
       breaksw
     case "-d":
       set wheretolook = d
       breaksw
     case "-u":
       set wheretolook = u
       breaksw
     case "-t":
       set wheretolook = t
       breaksw
     case "-z":
       set wheretolook = z
       breaksw
     case "-i":
       set wheretolook = i
       breaksw
     case "-l":
       set wheretolook = l
       breaksw
     case "-n":
       set wheretolook = n
       breaksw
     case "-s":
       set wheretolook = s
       breaksw
     case "-a":
       set wheretolook = a
       breaksw
     case "-x":
       set wheretolook = x
       breaksw
     case "-e":
       set wheretolook = e
       breaksw
     case "-S":
       set wheretolook = S
       breaksw
     case "-h":
       set wheretolook = h
       breaksw
     case "-L":
       set runlocal = y
       breaksw
     case "-h":
       echo 'Usage: findkey [-pbmgrdutielnsaxShzL] keywordfile [namefile findfile]'
       echo '       findkey [-P PIR_dataset] keywordfile [namefile findfile]'
       echo '       findkey [-G GenBank_dataset] keywordfile [namefile findfile]'
       breaksw
     default:
       if ($keyfile == "") then
         set keyfile = $a
       else 
         if ($namefile == "") then
            set namefile = $a
         else
           if ($findfile == "") set findfile = $a
         endif
       endif
       breaksw
     endsw
     @ index++
  end #while

  if ($keyfile == "") then
     echo 'No keyfile specified'
     exit 
  endif
   
else 
  #---------------------------------------------------------------
  # Interactive parameter input
  set complete = 0
  while ($complete == 0)
      #Display current parameter settings
      echo '___________________________________________________________________'
      echo '                     FINDKEY - Version 15 Feb 02'
      echo '    Please cite: Fristensky (1993) Nucl. Acids Res. 21:5997-6003'
      echo '___________________________________________________________________'
      echo 'Keyfile:' $keyfile
      echo 'Dataset:' $dbfile
      echo '-------------------------------------------------------------------'
      echo '   Parameter              Description                      Value'
      echo '-------------------------------------------------------------------'
      echo '1) Keyword     Keyword to find                               '$keyword
      echo '2) Keyfile     Get list of keywords from Keyfile'
      echo '3) WhereToLook p:PIR                                         '$wheretolook
      echo '               GenBank - b:bacterial     i:invertebrate'
      echo '                         m:mamalian      e:expressed seq. tag'
      echo '                         g:phage         l:plant'
      echo '                         r:primate       n:rna'
      echo '                         d:rodent        s:synthetic'
      echo '                         u:unannotated   a:viral'
      echo '                         t:vertebrate    x:patented'
      echo '                         z:STS           S:Genome Survey Seq.'
      echo '                         h:HTG'
      echo '               G: GenBank dataset        P: PIR dataset'
      echo '   -------------------------------------------------------------'
      echo '   Type number of your choice or 0 to continue:'

      #Prompt for parameter to change
      set paramnum = $<
      switch ($paramnum)
        case 0:
          if ($keyword != "" | $keyfile != "") then
             set complete = 1
          else echo '>>> Must specify value for Keyword or Keyfile'
          endif
          breaksw
        case 1:
          echo 'Type keyword to search for:'
          set keyword = $<
          set keyfile = ""
          breaksw
        case 2:
          echo 'Name of file containing keywords to search for:'
          set keyfile = $<
          set keyword = ""
          breaksw
        case 3:
          echo 'Choose one of {pevbimoglrndsuatGP}'
          set temp = $<
          if ($temp == p | $temp == b | $temp == i | $temp == m | $temp == e | \
              $temp == g | $temp == l | $temp == r | $temp == n | $temp == d | $temp == s | \
              $temp == u | $temp == a | $temp == t | $temp == z | $temp == G | $temp == P  | $temp == S | $temp == h ) then
             if ($temp == G | $temp == P) then
                echo 'Name of file containing user-defined dataset:'
                set dbfile = $<
                endif
             set wheretolook = $temp
          else echo '>>> Invalid choice'
          endif
          breaksw
        default:
        endsw
       #If parameter chosen is 0, and a minimal set of parameters are
      #set, terminate the loop (complete=1)
    end #while
endif

# findkey -L forces local execution, overriding XYLEM_RHOST. This 
# is necessary to prevent an infinite chain of calls to findkey
# on different hosts.
if (${?XYLEM_RHOST} && $wheretolook != G && $wheretolook != P && $runlocal != y) then 
#####################################################################
#  Run FINDKEY remotely, if XYLEM_RHOST and XYLEM_USERID are set
####################################################################
   # Remote hosts can be chosen by having a script called choosehost
   # in your bin directory. choosehost returns the name of a remote
   # host. While one possible implementation of choosehost is provided
   # with XYLEM, choosehost can be tailored to your particular
   # configuration of servers.
   if ($XYLEM_RHOST == choosehost) set XYLEM_RHOST = `choosehost`

   set tempname = "TMP_$jobid"
   set remotefn = $XYLEM_USERID@$XYLEM_RHOST\:$tempname
   if ($mode == interactive) echo "Copying keyword file to $remotefn.kw ..."
   if ($keyfile == "") then
    echo $keyword > $tempname.kw
    rcp $tempname.kw $remotefn.kw
    $RM $tempname.kw
   else
     rcp $keyfile $remotefn.kw
   endif
   if ($mode == interactive) echo "Running FINDKEY remotely on $XYLEM_RHOST as"\
                                  "user $XYLEM_USERID ..."
   rsh $XYLEM_RHOST -l $XYLEM_USERID findkey -L -$wheretolook $tempname.kw \
                                      $tempname.nam $tempname.fnd

   if $keyfile == "" then 
      set keyname = $keyword
   else
      set keyname = $keyfile:r
   endif

   if ($namefile == "") set namefile = $keyname.nam
   if ($findfile == "") set findfile = $keyname.fnd

   if ($mode == interactive) echo Copying $remotefn.nam to $namefile ...
   rcp $remotefn.nam $namefile
   if ($mode == interactive) echo Copying $remotefn.fnd to $findfile ...
   rcp $remotefn.fnd $findfile

   if ($mode == interactive) echo Removing temporary files...
   rsh $XYLEM_RHOST -l $XYLEM_USERID $RM  "$tempname.*"

else 
#####################################################################
#  For a given database, search annotation files for all entries
#  containing the specified keyword(s).
####################################################################

  # If only one keyword is present in keyfile, copy it to $keyword
  # This will cause FINDKEY to use egrep, which is faster than fgrep.
  if ($keyfile != "") then
     if (`wc -l <$keyfile` == 1) then
      set keyword =  `cat $keyfile` 
      set keyfile = ""
     endif
  endif

if ($wheretolook == p) then
    #------------------------ PIR/NBRF --------------------------
    set divisions = (pir1 pir2 pir3 pir4)
    foreach div ($divisions)
      set base = $PIR/$div
      if ($mode == interactive) echo "Searching $base.ano..."
      # egrep through the .ano file
      if ($keyfile == "") then
        set key = $keyword
        nice egrep -i -n -e $keyword $base.ano > $jobid.grep
      else
        set key = $keyfile:r
        nice fgrep -i -n -f $keyfile $base.ano > $jobid.grep
      endif

      if ($namefile == "") set namefile = $key~pir.nam
      if ($findfile == "") set findfile = $key~pir.fnd

      # identify each sequence found
      # temporarily store in $jobid..nam & $jobid..fnd, and then append
      # to $namefile and $findfile
      if (-z $jobid.grep) then
        if ($mode == interactive) echo "No matches found in $base.ano" 
      else
        if ($mode == interactive) then
           echo "Sequence names will be written to $namefile"
           echo "Lines containing keyword(s) will be written to $findfile"
        endif
        $NICE identify $jobid.grep $base.ind $jobid.nam $jobid.fnd
        cat $jobid.nam >> $namefile
        cat $jobid.fnd >> $findfile
        $RM $jobid.* 
      endif
     end

else

if ($wheretolook == G | $wheretolook == P) then
    #---------------- User-defined dataset ------------------
      # If dataset is not split, split it 
      set base = $dbfile:r
      set dbextension = $dbfile:e
      if ($dbextension == "gen") then #GenBank
         set needtosplit = true
         set base = TMP$jobid
         splitdb $dbfile $base.ano $base.wrp $base.ind
      else
        if ($dbextension == "pir") then #PIR
           set needtosplit = true
           set base = TMP$jobid
           splitdb -p $dbfile $base.ano $base.wrp $base.ind
        else
          set needtosplit = false
        endif
      endif

      if ($mode == interactive) echo "Searching $base.ano..."
      # egrep through the .ano file
      if ($keyfile == "") then
        set key = $keyword
        nice egrep -i -n -e $keyword $base.ano > $jobid.grep
      else
        set key = $keyfile:r
        nice fgrep -i -n -f $keyfile $base.ano > $jobid.grep
      endif

      if ($namefile == "") set namefile = $key~$base.nam
      if ($findfile == "") set findfile = $key~$base.fnd

      if (-z $jobid.grep) then
        echo "No matches found in $base.ano" 
      else
        if ($mode == interactive) then
           echo "Sequence names will be written to $namefile"
           echo "Lines containing keyword(s) will be written to $findfile"
        endif
        $NICE identify $jobid.grep $base.ind $namefile $findfile
      endif    
      if ($needtosplit == true) $RM $base.ano $base.wrp $base.ind
else
#------------------------  GenBank  ------------------
  # Set $div to the name of the database division to search
  switch ($wheretolook)
     case b:
       set div = bct
       breaksw
     case m:
       set div = mam
       breaksw
     case g:
       set div = phg
       breaksw
     case r:
       set div = pri
       breaksw
     case d:
       set div = rod
       breaksw
     case u:
       set div = una
       breaksw
     case t:
       set div = vrt
       breaksw
     case z:
       set div = sts
       breaksw
     case i:
       set div = inv
       breaksw
     case l:
       set div = pln
       breaksw
     case n:
       set div = rna
       breaksw
     case s:
       set div = syn
       breaksw
     case a:
       set div = vrl
       breaksw
     case x:
       set div = pat
       breaksw
     case e:
       set div = est
       breaksw
     case S:
       set div = gss
       breaksw
     case h:
       set div = htg
       breaksw
     default:
    endsw

  # $base is the name of the database file, without the extension
  # Most GenBank divisions are present in one file eg. gbrna.
  # Large GenBank divisions such as EST and Primate are split
  # eg. gbest1, gbest2, gbest3...
  # Regardless of how many divisions there are, BASESET creates
  # the list of all files  for that division.
    if ($div == vec) then
       set BASESET = $VEC/vecbase
    else
        if (-e $GB/gb$div.ind) then
           set BASESET = $GB/gb$div
        else
           set index = 1
           set BASESET = ()
           while (-e $GB/gb$div$index.ind)
             set BASESET = ($BASESET $GB/gb$div$index)
             @ index++
             end # while     
        endif
    endif 
   
  #use grep to find the keyword(s), and then use identify to determine which 
  #entries correspond to the occurrences of the keyword(s) found

  foreach division ($BASESET)
    if (-e $division.ind) then
       set base = $division
       if ($mode == interactive) echo "Searching $base.ano..."
       if ($keyfile == "") then
          set key = $keyword
          nice egrep -i -n -e $keyword $base.ano > $jobid.grep
       else
          set key = $keyfile:r
          nice fgrep -i -n -f $keyfile $base.ano > $jobid.grep
       endif
       if ($namefile == "") set namefile = $key~$div.nam
       if ($findfile == "") set findfile = $key~$div.fnd

       if (-z $jobid.grep) then
                if ($mode == interactive) echo "No matches found in $base.ano" 
       else
           if ($mode == interactive) then
              echo "Sequence names will be written to $namefile"
              echo "Lines containing keyword(s) will be written to $findfile"
           endif
           $NICE identify $jobid.grep $base.ind $jobid.nam $jobid.fnd
       endif
       cat $jobid.nam >> $namefile
       cat $jobid.fnd >> $findfile
     endif
   end    
endif
endif
endif

#####################################################################
#  Clean up.
####################################################################
$RM $jobid.* 
