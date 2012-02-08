#PRINTDOC - a c-shell script (version  30 Apr 93)
#PRINTDOC  determines from the file extension what sort of file is
#          to be printed. It then generates the appropriate print
#          command for the PostScript printer
# 
# Synopsis:  printdoc <filename>

set filename = $1
set extension = $filename:e

switch ($extension)
  #PostScript
  case "ps"
    lpr -Pccp $filename
    breaksw

  #nroff
  case "nroff"
    nroff -ms < $filename | xlp -f BDJ1 -d xerox    
    breaksw

  #man (same as nroff)
  case "man"
    nroff -ms < $filename | xlp -f BDJ1 -d xerox    
    breaksw

  #nroff.me
  case "me"
    nroff -me < $filename | xlp -f BDJ1 -d xerox    
    breaksw

  # .doc files print at 12cpi with a 1in. left margin 
  case "doc"
    pr -f $filename | xlp -f BDJ1 -d xerox
    breaksw

  # .tex files contain page breaks to begin each page 
  # header and footer supplied by pr are suppressed to
  # accommodate this type of 'preformatting'
  case "tex"
    /usr/5bin/pr -t -o5 $filename | xlp -f BDJ1 -d xerox
    breaksw

  case "txt"
    pr -f $filename | xlp -f BDJ1 -d xerox
    breaksw

  default
    pr -f $filename | xlp -f BDJ1 -d xerox
    breaksw
  endsw
