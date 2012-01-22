  (* UPDATE HISTORY
      2 Sep 2001 GenBank has shortened LOCUS names to 16 char.
     15 Aug 2001 Revised to allow LOCUS names up to 18 char., to
     conform with upcomming changes in GenBank LOCUS line. Also
     removed support for the ancient LIMB database.
     28 Mar 98 added -t option to include genus & species in .ind file 
     26 Aug 95 Writes ACCESSION number to INDFILE using 8 characters.
     30 May 95 Changed CRUNCHOFFSET to ASCII 33 = "!". VT100 emulators
     didn't like ASCII 29.
     21 Dec 94 -c option compresses leading blanks in .ano file. 
     This version of splitdb now looks first for an ACCESSIONS line, and if
     none is found, it looks for the first #accession line.
  *)
 
  (***********************************************************)
  (*                                                         *)
  (*  SPLITDB   VERSION    9/ 2/2001  Standard Pascal        *)
  (*            Brian Fristensky                             *)
  (*            Dept. of Plant Science                       *)
  (*            University of Manitoba                       *)
  (*            Winnipeg, MB Canada R3T 2N2                  *)
  (*                                                         *)
  (*  Splits GENBANK file into three files:                  *)
  (*    ANOFILE - contains annotation parts of entries       *)
  (*    SEQFILE - contains sequence parts of entries         *)
  (*    INDFILE - contains linenumbers for beginning of each *)
  (*              entry in ANOFILE and SEQFILE               *)
  (*                                                         *)
  (*           -c compress leading blanks in .ano file       *)
  (*           -t append the contents of the first ORGANISM  *)
  (*              line to each line in INDFILE               *)
  (*                                                         *)
  (*  Copyright (c) 1988 - 2001      by Brian Fristensky     *)
  (*  !!! in comment indicates feature which may need change *)
  (***********************************************************)
  program SPLITDB(DBFILE,ANOFILE,SEQFILE,INDFILE);
(*!!!  Some Pascals require file parameters in program heading *)

  const MAXLINE = 132;
        MAXWORD = 25;
	SEQLINELEN=75; (* length of sequence line in SEQFILE *)
        CRUNCHOFFSET = 33;
        
(* BEGIN MODULE STARTARGNUM *)
	STARTARGNUM=1;    (* SUN Pascal: ARG(1) is 1st command line argument*)
      (*STARTARGNUM=2;       ATT Pascal: ARG(2) is 1st command line argument*)
(* END MODULE STARTARGNUM         VERSION= 'SUNMODS     Version  8/ 9/94'; *)

  type  DBTYPE = (GB,PIR,EMBL,VECTOR);
       
(* BEGIN MODULE TYPE.WORD *)
       (*   <word>::= <non-blank char>[<non-blank char>] *)
       WORD    = record
                 LEN:integer;
                 STR:array[1..MAXWORD] of char
                 end;
(* END MODULE TYPE.WORD         VERSION= 'SUNMODS     Version  8/ 9/94'; *)
 
(* BEGIN MODULE TYPE.LINE *)
       CHARARRAY = packed array[1..MAXLINE] of char;
       LINE = record
                STR:CHARARRAY;
                LEN:integer
                end;
(* END MODULE TYPE.LINE         VERSION= 'SUNMODS     Version  8/ 9/94'; *)
 
  var DBFILE,ANOFILE,SEQFILE,INDFILE:text;
      NAMEID,DEFID,ACID,SEQID,SOURCEID,ORGID, (* keyword identifiers *)
      OLDACID,REFID,  (*can be eliminated after PIR changes *)
      NAME,ACCESSION:WORD;     (* name and accession numbers *)
      ORGANISM:LINE;
      CRUNCHFLAG,        (* =chr(CRUNCHOFFSET), indicates -c compression *)
      CH: char;
      DATABASE:DBTYPE;
      TAXONOMY,            (* =true, include genus & sp. in index file *)
      COMPRESSION:boolean; (* =true, then compress leading blanks *)
      CLINE,TITLE:LINE;  (* current line, title line for SEQFILE *)
      POSITION,          (* next unprocessed position in CLINE *)
      CANO,CSEQ,         (* current line pointers for ANOFILE & SEQFILE *)
      FIRSTANO,          (* first annotation line for an entry in ANOFILE*)
      FIRSTSEQ,          (* title line for a sequence in SEQFILE *)
      ARGNUM,            (* number of command line argument *) 
      ANOARGNUM          (* arg # of annotation file, used for rewrite *)
      :integer;

      ANOFILENAME:packed array [1..132] of char;

  (***************************************************************)
  (* Read options from command line and set DATABASE             *)
  (***************************************************************)
   procedure READOPTIONS(var ARGNUM:integer; var DATABASE:DBTYPE;
                         var COMPRESSION,TAXONOMY:boolean);
     var ARGUMENT:packed array[1..132] of char;
         OPTIONSDONE:boolean;
     begin (* READOPTIONS *)
       ARGNUM:=STARTARGNUM;
       DATABASE:=GB; TAXONOMY:= false; COMPRESSION:=false; OPTIONSDONE:=false;
       repeat
               argv(ARGNUM,ARGUMENT);
	       if ARGUMENT[1]='-' then begin
		 if ARGUMENT[2] in ['g','p','e','v','t','c'] then 
		   case ARGUMENT[2] of 
		    'g':DATABASE:=GB;
		    'p':DATABASE:=PIR;
		    'e':DATABASE:=EMBL;
		    'v':DATABASE:=VECTOR;
		    't':TAXONOMY:=true;
		    'c':COMPRESSION:=true
		     end; (*case*)
		 ARGNUM:=ARGNUM+1
		 end
	       else OPTIONSDONE:=true
        until OPTIONSDONE
     end; (* READOPTIONS *)      
 
  (***************************************************************)
  (* Initialize keywords. Keywords identify the line on which a  *)
  (* data item occurs.  NAMEID is the keyword that identifies the*)
  (* line on which the name of the entry occurs (eg. 'LOCUS' for *)
  (* GenBank, 'ENTRY' for PIR and so forth.)  DEFID is the TITLE,*)
  (* or DEFINITION line, to be written to SEQFILE.  ACID is the  *)
  (* ACCESSION number line. SEQID denotes where the sequence     *)
  (* begins. Not all databases have all types of keywords, and   *)
  (* the order of occurrence varies with database. SOURCEID and  *)
  (* ORGID are used for GenBank SOURCE and ORGANISM lines.       *)
  (***************************************************************)
   procedure INITKEYS(var NAMEID,DEFID,ACID,SEQID,SOURCEID,ORGID:WORD);
     begin
        with NAMEID do 
          case DATABASE of
             GB,VECTOR:begin
                  STR[1]:='L';STR[2]:='O';STR[3]:='C';STR[4]:='U';STR[5]:='S';
                  LEN:=5
                  end; (* GB,VECTOR *)
             PIR:begin
                  STR[1]:='E';STR[2]:='N';STR[3]:='T';STR[4]:='R';STR[5]:='Y';
                  LEN:=5
                  end; (* PIR *)
             EMBL:begin
                  STR[1]:='I';STR[2]:='D';
                  LEN:=2
                  end (* EMBL *)
             end; (*case*)

        with DEFID do 
          case DATABASE of
             GB,VECTOR:begin
                  STR[1]:='D';STR[2]:='E';STR[3]:='F';STR[4]:='I';STR[5]:='N';
                  STR[6]:='I';STR[7]:='T';STR[8]:='I';STR[9]:='O';STR[10]:='N';
                  LEN:=10
                  end; (* GB,VECTOR *)
             PIR:begin
                  STR[1]:='T';STR[2]:='I';STR[3]:='T';STR[4]:='L';STR[5]:='E';
                  LEN:=5
                  end; (* PIR *)
             EMBL:begin
                  STR[1]:='D';STR[2]:='E';
                  LEN:=2
                  end (* EMBL *)
             end; (*case*)

        with ACID do 
          case DATABASE of
             GB,VECTOR:begin
                  STR[1]:='A';STR[2]:='C';STR[3]:='C';STR[4]:='E';STR[5]:='S';
                  STR[6]:='S';STR[7]:='I';STR[8]:='O';STR[9]:='N';
                  LEN:=9
                  end; (* GB,VECTOR *)
             PIR:begin
                  STR[1]:='#';STR[2]:='a';STR[3]:='c';STR[4]:='c';STR[5]:='e';
                  STR[6]:='s';STR[7]:='s';STR[8]:='i';STR[9]:='o';STR[10]:='n';
                  LEN:=10
                  end; (* PIR *)
             EMBL:begin
                  STR[1]:='A';STR[2]:='C';
                  LEN:=2
                  end (* EMBL *)
             end; (*case*)

        with SEQID do 
          case DATABASE of
             GB,VECTOR:begin
                  STR[1]:='O';STR[2]:='R';STR[3]:='I';STR[4]:='G';STR[5]:='I';
                  STR[6]:='N';
                  LEN:=6
                  end; (* GB,VECTOR *)
             PIR:begin
                  STR[1]:='S';STR[2]:='E';STR[3]:='Q';STR[4]:='U';STR[5]:='E';
                  STR[6]:='N';STR[7]:='C';STR[8]:='E';
                  LEN:=8
                  end; (* PIR *)
             EMBL:begin
                  STR[1]:='S';STR[2]:='Q';
                  LEN:=2
                  end (* EMBL *)
             end; (*case*)

       (* Used for GenBank SOURCE and ORGANISM lines. *)
       with SOURCEID do begin
          STR[1]:='S';STR[2]:='O';STR[3]:='U';STR[4]:='R'; STR[5]:='C';
          STR[6]:='E';
          LEN:=6
          end;
       with ORGID do begin
          STR[1]:='O';STR[2]:='R';STR[3]:='G';STR[4]:='A'; STR[5]:='N';
          STR[6]:='I';STR[7]:='S';STR[8]:='M';
          LEN:=8
          end;


       (* These keywords can be eliminated after PIR makes the switch over
          to the #accession protocol *)
       with OLDACID do begin
          STR[1]:='A';STR[2]:='C';STR[3]:='C';STR[4]:='E'; STR[5]:='S';
          STR[6]:='S';STR[7]:='I';STR[8]:='O';STR[9]:='N'; STR[10]:='S';
          LEN:=10
          end;
       with REFID do begin
          STR[1]:='R';STR[2]:='E';STR[3]:='F';STR[4]:='E';STR[5]:='R';
          STR[6]:='E';STR[7]:='N';STR[8]:='C';STR[9]:='E';
          LEN:=9
         end
     end; (* INITKEYS *)

      
  (***************************************************************)
  (* I/O procedures                                              *)
  (***************************************************************)
   procedure OPENFILES(var DBFILE,ANOFILE,SEQFILE,INDFILE:text;
                       var ANOARGNUM:integer);

(* BEGIN MODULE FILEARGS *)
  (* This procedure overcomes one of the stupidest aspects of UNIX Pascal,
     namely the fact that filenames in the program statement are supposed to
     be actual UNIX filenames!  To overcome this, the 2-argument version of
     reset and rewrite must be used with string variables.  This module
     need only contain the reset and rewrite statements in any normal 
     implementation of Pascal. *)
  procedure FILEARGS(var F:text; FTYPE:char; var ARGNUM:integer);
    var  ARGUMENT : packed array[1..132] of char;
    begin
      argv(ARGNUM,ARGUMENT);
      if FTYPE='I' then reset(F,ARGUMENT)
      else rewrite(F,ARGUMENT);
      ARGNUM:=ARGNUM+1
    end; (* FILEARGS *)
(* END MODULE FILEARGS         VERSION= 'SUNMODS     Version  8/ 9/94'; *)

     begin
       FILEARGS(DBFILE,'I',ARGNUM);

       (* We need to rewrite ANOFILE later, so save the ARGNUM as ANOARGNUM *)
       ANOARGNUM:=ARGNUM;
       FILEARGS(ANOFILE,'O',ARGNUM);
      
       if DATABASE in [GB,PIR,EMBL,VECTOR] then FILEARGS(SEQFILE,'O',ARGNUM);
       FILEARGS(INDFILE,'O',ARGNUM)
     end; (* OPENFILES *)


(* BEGIN MODULE READLINE *)
  (* Read a line from a file, omitting trailing blanks *)  
  procedure READLINE(var F:text; var L:LINE);
    var LASTNONBLANK,I:integer;
        CH:char;
    begin
      with L do begin
        LEN:=0; LASTNONBLANK:=0;
        while not eoln(F) do begin
          read(F,CH);
          if LEN < MAXLINE then begin 
            LEN:=LEN+1; STR[LEN]:=CH;
            if CH <> ' ' then LASTNONBLANK:=LEN
            end
          end;
        if not eof(F) then readln(F);
        LEN:=LASTNONBLANK;
        for I:= LEN+1 to MAXLINE do STR[I]:=' '
        end; (* with L*) 
    end; (* READLINE *)
(* END MODULE READLINE         VERSION= 'SUNMODS     Version  8/ 9/94'; *)

(* BEGIN MODULE WRITELINE *)
   (*  Write a line to a file using L char, left-justified.  *)
   procedure WRITELINE(var F:text;W:LINE;L:integer);
     var I :integer;
     begin
       with W do
         for I := 1 to L do
           if I <= LEN then write(F,STR[I])
           else write(F,' ')
     end; (* WRITELINE *)
(* END MODULE WRITELINE         VERSION= 'SUNMODS     Version  8/ 9/94'; *)
    
(* BEGIN MODULE GETWORD *)
  (*  Read a word from LINE L *)
  procedure GETWORD(CLINE:LINE;var W:WORD;var P:integer);
    var I : integer;
    begin
      with W do begin
        LEN:=0;
        while CLINE.STR[P]=' ' do P:=P+1; 
        while CLINE.STR[P] <> ' ' do begin 
          if LEN < MAXWORD then begin
            LEN := LEN + 1;
            STR[LEN]:= CLINE.STR[P]
            end;
          P:= P+1
          end;
        for I := LEN+1 to MAXWORD do STR[I]:= ' '
        end
    end; (* GETWORD *)
(* END MODULE GETWORD *)

(* BEGIN MODULE WRITEWORD *)
   (*  Write a word to a file using L char, left-justified.  *)
   procedure WRITEWORD(var F:text;W:WORD;L:integer);
     var I :integer;
     begin
       with W do
         for I := 1 to L do
           if I <= LEN then write(F,STR[I])
           else write(F,' ')
     end;
(* END MODULE WRITEWORD         VERSION= 'SUNMODS     Version  8/ 9/94'; *)

(* BEGIN MODULE SAMEWORD *)
      (* Compare two WORDS for equality *)
      function SAMEWORD(var W1,W2:WORD):boolean;
        var I:integer;
            T:boolean;
        begin
          if W1.LEN = W2.LEN then begin
            T:=true;I:=1;
            while (I <= W1.LEN) and T do
              if W1.STR[I] = W2.STR[I] then I:=I+1
              else T:=false;
            SAMEWORD:=T
            end 
          else SAMEWORD:=false
        end; (* SAMEWORD *)
(* END MODULE SAMEWORD         VERSION= 'SUNMODS     Version  8/ 9/94'; *)


  (***************************************************************)
  (* Read in each line and write to output file, until KEYWORD   *)
  (* is found starting in the column specified by POSITION       *)
  (* When done, POSITION returns the column containing the next  *)
  (* unread character in CLINE, after KEYWORD is read.           *)
  (***************************************************************)
  procedure ADVANCE(var CLINE:LINE; KEYWORD:WORD; var POSITION,LINENUM:integer);
    var FLAG:char;
        FLAGPOSN:integer;
        ID:WORD;

      (* If the number of leading blanks is greater than 3, represent*)
      (* the leading blanks by <CRUNCHFLAG><CRUNCHOFFSET+NUMBLANKS>  *)
      (* where CRUNCHFLAG is the character used to signal that the   *)
      (* line contains compressed blanks, and CRUNCHOFFSET+NUMBLANKS *)
      (* is a character whose ASCII value is equal to number of      *)
      (* leading blanks that the line should have plus CRUNCHOFFSET. *)
      (* CRUNCHOFFSET is added to prevent the use of the control     *)
      (* characters with low ASCII values. Since none of the data-   *)
      (* bases commonly used with XYLEM should have more then 80     *)
      (* leading blanks, ASCII 127 should never be exceeded.         *)
      procedure WRITECRUNCH(var F:text; var CLINE:LINE);
	const MINBLANKS = 3;
	var NUMBLANKS,I:integer;
	begin
	  with CLINE do
	    if LEN >= MINBLANKS then
	      if (STR[1]=' ') and (STR[2]=' ') and (STR[3]=' ') then begin
		 write(F,CRUNCHFLAG);
		 NUMBLANKS:=3; I:=4;
		 (* Set NUMBLANKS = the number of leading blanks. *)
		 (* Write the corresponding ASCII character. *)
		 while (STR[I]=' ') and (I<=LEN) do begin
		    NUMBLANKS:= NUMBLANKS +1;
		    I:= I + 1
		    end;
		 write(F,chr(CRUNCHOFFSET+NUMBLANKS));
    
		 (* Write the non-blank portion of the line. *)
		 for I:= NUMBLANKS+1 to LEN do begin
		     write(F,STR[I]);
		     end
		 end
	     else WRITELINE(F,CLINE,CLINE.LEN)
	   end; (* WRITECRUNCH *)

    begin (* ADVANCE *)
      FLAGPOSN:=POSITION;
      FLAG:=KEYWORD.STR[1];     
      repeat
        repeat 
          READLINE(DBFILE,CLINE);
          if COMPRESSION then WRITECRUNCH(ANOFILE,CLINE)
          else WRITELINE(ANOFILE,CLINE,CLINE.LEN);
          writeln(ANOFILE);
          LINENUM:=LINENUM+1
        until eof(DBFILE) or (CLINE.STR[FLAGPOSN]=FLAG);
        POSITION:=FLAGPOSN;
        GETWORD(CLINE,ID,POSITION)
      until eof(DBFILE) or SAMEWORD(KEYWORD,ID)
    end; (* ADVANCE *)

  (******************************************************************)
  (* Transfer lines until the DEFINITION (or equivalent) is reached.*)
  (******************************************************************)
  procedure READDEF(var CLINE:LINE; var CANO,CSEQ:integer);
    var I,POSITION:integer;
    begin
         (* Get DEFINITION line *)
         POSITION:=1;
         ADVANCE(CLINE,DEFID,POSITION,CANO);

        (* Read in partial title from DEFINITION line *)
         with CLINE do 
           while (STR[POSITION]=' ') and (POSITION<LEN) do POSITION:=POSITION+1;
         I:=0;
         while (I<50) and (POSITION<=CLINE.LEN) do begin
           I:=I+1;
           TITLE.STR[I]:=CLINE.STR[POSITION];
           POSITION:=POSITION+1
           end;
         TITLE.LEN:= I;

         (* Title line:
            >name - parital title from DEFINITION line *)
         write(SEQFILE,'>');
         WRITEWORD(SEQFILE,NAME,NAME.LEN);
         write(SEQFILE,' - ');
         WRITELINE(SEQFILE,TITLE,TITLE.LEN);
         writeln(SEQFILE);
         CSEQ:=CSEQ+1;
         FIRSTSEQ:=CSEQ;
     end; (* READDEF *)

  (***************************************************************)
  (* Transfer lines until the ACCESSION line is found. Read      *)
  (* ACCESSION from CLINE.                                       *)
  (***************************************************************)
  procedure READAC(var CLINE:LINE; var ACCESSION:WORD; var CANO:integer);
    var POSITION:integer;
    (* Search for ACCESSION line *)
    procedure FINDAC(var CLINE:LINE; var POSITION,LINENUM:integer);
      var GOTIT:boolean;
          FLAGPOSN:integer;
          ID:WORD;
      begin
        FLAGPOSN:=1;
        (* Keep searching until either an ACCESSIONS or REFERENCE line is
           found. *)    
        GOTIT:= false; 
        repeat
          repeat 
            READLINE(DBFILE,CLINE);
            WRITELINE(ANOFILE,CLINE,CLINE.LEN);
            writeln(ANOFILE);
            LINENUM:=LINENUM+1
          until eof(DBFILE) or (CLINE.STR[FLAGPOSN] in ['A','R']);
          POSITION:=FLAGPOSN;
          GETWORD(CLINE,ID,POSITION);
          if SAMEWORD(OLDACID,ID) then GOTIT:= true
          else if SAMEWORD(REFID,ID) then GOTIT:= true
        until eof(DBFILE) or GOTIT;

        if SAMEWORD(REFID,ID) then begin
        (* We're on the REFERENCE line. No ACCESSIONS line was found. Move to
           the first #accession line. *)
           POSITION:=4;
           ADVANCE(CLINE,ACID,POSITION,CANO)
           end
      end; (* FINDAC *)


    begin
         (* Get ACCESSION line *)
         if DATABASE = PIR then begin
            (* Since the primary accession number could be on an ACCESSIONS
               line (old format) or an #accession line (new format) we
               have be able to find either. In either case, FINDAC will
               set POSITION to the first character of the accession num. *)
           FINDAC(CLINE,POSITION,CANO);
           GETWORD(CLINE,ACCESSION,POSITION)
           end (* PIR *)
            
         else begin
           POSITION:=1;
           ADVANCE(CLINE,ACID,POSITION,CANO);
           GETWORD(CLINE,ACCESSION,POSITION)
           end;
         (* PIR - If there were secondary accession numbers, the
            primary accession number will be terminated by a semicolon (;)
            which must be removed *)
         if DATABASE = PIR then
            with ACCESSION do 
               if STR[LEN] = ';' then LEN:= LEN-1
     end; (* READAC *)

  (******************************************************************)
  (* Transfer lines until the SOURCE line is reached. Source may span*)
  (* more than 1 line.                                               *)
  (* ORGANISM line may take the form                                 *)
  (* [genome] genus species subspecies                               *)
  (* so we just read in the rest of the line.                        *)
  (******************************************************************)
  procedure READSOURCE(var CLINE:LINE; var CANO:integer;
                       var ORGANISM:LINE);
    var I,J:integer;
    begin
         (* Get SOURCE line *)
(*         I:=1;
         ADVANCE(CLINE,SOURCEID,I,CANO);
*)
         (* Read in ORGANISM line and write to .ano file. *)
         I:=3;
         ADVANCE(CLINE,ORGID,I,CANO);

         (* Copy remaining part of line to ORGANISM *)
         with CLINE do begin
            while ((STR[I] = ' ') and (I <= LEN)) do  I:= I + 1;
            J:=1;
            while I <= LEN do begin
              ORGANISM.STR[J] := STR[I];
              I:= I+1; J:= J + 1
              end;
            ORGANISM.LEN:= J-1
            end (* with CLINE *)

     end; (* READSOURCE *)

  (***************************************************************)
  (* Transfer lines until sequence is found. Write sequence to   *)
  (* SEQFILE, stripping off blanks and numbers.                  *)
  (***************************************************************)
  procedure READSEQ(var CLINE:LINE; var CSEQ:integer);
    var WIDTH,POSITION:integer;
    begin
         (* Advance to ORIGIN line. Sequence begins on next line *)
         POSITION:=1;
         ADVANCE(CLINE,SEQID,POSITION,CANO);
         if DATABASE=PIR then writeln(ANOFILE,'///')
         else writeln(ANOFILE,'//');
         CANO:=CANO+1;

         (* Read in sequence and write to output. *)
         (* Omit extra blanks and numbers. *)
         WIDTH:=0;
         while not (DBFILE^='/') do begin
           while not eoln(DBFILE) do begin
             read(DBFILE,CH);
             if not (CH in [' ','0'..'9']) then begin
                write(SEQFILE,CH);
                WIDTH:=WIDTH+1;
                if WIDTH=SEQLINELEN then begin writeln(SEQFILE); 
                                               CSEQ:=CSEQ+1;  
                                               WIDTH:=0
                                               end
                end (* CH in [' ','0'..'9'] *)
             end; (* not eoln *)
           readln(DBFILE)
           end; (* not (DBFILE^='/') *)
	 if WIDTH>0 then begin
            writeln(SEQFILE);
            CSEQ:=CSEQ+1
            end;
         if not eof(DBFILE) then readln(DBFILE);

         (* PIR files terminate with '\\\'. However, splitdb does not
            demand that this terminator line be there. It is not re-
            generated by getloc. *)
         if DBFILE^='\' then readln(DBFILE);

         (* Read past trailing blank lines, if any *)
         while not(eof(DBFILE)) and (DBFILE^=' ') do readln(DBFILE); 
     end; (* READSEQ *)

       
  (* ----------------------------------------------------------  *)
  (* ----------------- MAIN  PROCEDURE  -----------------------  *)
    begin
       (* Determine which database is being processed. *)
       READOPTIONS(ARGNUM,DATABASE,COMPRESSION,TAXONOMY);

       (* Initialize keywords according to database type. *)
       INITKEYS(NAMEID,DEFID,ACID,SEQID,SOURCEID,ORGID);    

       (* Open files. *)
       OPENFILES(DBFILE,ANOFILE,SEQFILE,INDFILE,ANOARGNUM); 
       argv(ANOARGNUM,ANOFILENAME); (* save ANOFILENAME for rewrite *)

       CANO:=0; CSEQ:=0;

       (* Advance to first entry and throw away header information
          by rewriting ANOFILE.*)
       POSITION:=1;
       CRUNCHFLAG:= chr(CRUNCHOFFSET);
       ADVANCE(CLINE,NAMEID,POSITION,CANO);
      
       rewrite(ANOFILE,ANOFILENAME);
    
       WRITELINE(ANOFILE,CLINE,CLINE.LEN);writeln(ANOFILE);
       CANO:=1; 
       GETWORD(CLINE,NAME,POSITION);
       WRITEWORD(INDFILE,NAME,16);
       FIRSTANO:=CANO;
       (* Note on integer output to INDFILE: *)
       (* ATT Pascal writes integers using a field width of 10 unless a field
          width is set.  If the field width is too small to print the integer,
          the minimum required field width is used.  Thus, I have to set a 
          field width to insure that the minimum field width is used in 
          all cases, if I want to minimize the size of the output file.*)

       while not eof(DBFILE) do begin 

         case DATABASE of
           GB,PIR,VECTOR:begin
                  READDEF(CLINE,CANO,CSEQ);
                  READAC(CLINE,ACCESSION,CANO);
                  if (TAXONOMY and (DATABASE=GB))
                     then READSOURCE(CLINE,CANO,ORGANISM);
                  READSEQ(CLINE,CSEQ);
                  write(INDFILE,' ');
                  WRITEWORD(INDFILE,ACCESSION,8);
                  write(INDFILE,' ',FIRSTANO:8,' ',FIRSTSEQ:8);
                  if (TAXONOMY and (DATABASE=GB))
                     then begin
                          write(INDFILE,' '); 
                          WRITELINE(INDFILE,ORGANISM,ORGANISM.LEN)
                          end;
                  writeln(INDFILE) 
                  end; (* GB,PIR,VECTOR *)
           EMBL:begin
                  READAC(CLINE,ACCESSION,CANO);
                  READDEF(CLINE,CANO,CSEQ);
                  READSEQ(CLINE,CSEQ);
                  write(INDFILE,' ');
                  WRITEWORD(INDFILE,ACCESSION,8);
                  writeln(INDFILE,' ',FIRSTANO:8,' ',FIRSTSEQ:8)
                end (* EMBL *)

           end; (* case *)

         (* Advance to next LOCUS line. Read in NAME. *)
         (* Write NAME to INDFILE and set FIRSTANO to current line number*)
         if not eof(DBFILE) then begin
           POSITION:=1;
           ADVANCE(CLINE,NAMEID,POSITION,CANO);
           GETWORD(CLINE,NAME,POSITION);
           WRITEWORD(INDFILE,NAME,16);
           FIRSTANO:=CANO
           end
         end; (* while not eof *)
 
    (*!!!CLOSE(ANOFILE); CLOSE(SEQFILE); CLOSE(INDFILE)*)
    end. (* SPLITDB  *)
