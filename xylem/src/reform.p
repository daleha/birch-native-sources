  (***********************************************************)
  (*                                                         *)
  (*  REFORM    VERSION   4/ 1/2000  Standard Pascal           *)
  (*            Brian Fristensky                             *)
  (*            Dept. of Plant Science                       *)
  (*            University of Manitoba                       *)
  (*            Winnipeg, MB R3T 2N2  CANADA                 *)
  (*                                                         *)
  (*  SYNOPSIS                                               *)
  (*  reform [-gpcnm] [-ln] [-sn] file                       *)
  (*                           or                            *)
  (*  ralign file window wordsize mm in del |reform [options]*)
  (*                                                         *)
  (*  DESCRIPTION                                            *)
  (*  Reformats multiple alignment output                    *)
  (*                                                         *)
  (*      g    Gaps are to be represented by dashes (-).     *)
  (*      p    Bases which agree with the consensus are      *)
  (*           represented by periods (.).                   *)
  (*      c    Positions at which all sequences agree are    *)
  (*           capitalized in the consensus.                 *)
  (*      n    Sequence data is nucleic acid. Defult protein.*)
  (*      fx   Specify input file format, where x is:        *)
  (*           r:RALIGN  p:PEARSON  i:INTELLIGENETICS (MASE) *)
  (*      m    each input sequence may run over many lines.  *)
  (*           Sequences begin with >name, after Pearson.    *)
  (*           Obsolete, equivalent to -fp                   *)
  (*      ln   The output linelength is set to n.            *)
  (*           Default is 70.                                *)
  (*      sn   numbering starts with n (0 default).          *)
  (*                                                         *)
  (*    file   Sequence file as described in ralign docu-    *)
  (*           mentation.  reform needs to re-read the       *)
  (*           sequence file read by ralign to get the       *)
  (*           names of the sequences, which ralign ignores. *)
  (*                                                         *)
  (*  Copyright (c) l988-2000  by Brian Fristensky           *)
  (*  !!! in comment indicates feature which may need change *)
  (***********************************************************)
  program REFORM (input, output (*,SEQFILE*));
(*!!!  Some Pascals require file parameters in program heading *)

  const MAXLINE = 80;   (* Maximum length of input line *)
        MAXWORD = 25;   (* Maximum length of word eg. seq. name *)
        MAXSEQ  = 250;  (* Maximum number of sequences in alignment *)
        MAXLEN  = 25000; (* Maximum length of sequences in alignment *)
(* BEGIN MODULE STARTARGNUM *)
	STARTARGNUM=1;    (* SUN Pascal: ARG(1) is 1st command line argument*)
      (*STARTARGNUM=2;       ATT Pascal: ARG(2) is 1st command line argument*)
(* END MODULE STARTARGNUM         VERSION= 'SUNMODS     Version  8/ 9/94'; *)
 
  type 
        ALIGNTYPE = (RALIGN, IG, PEARSON); (* alignment type *)
   
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

 
  var (* Names, sequences, and related parameters *)
      NAMES         :array[0..MAXSEQ] of WORD;  
      ALIGNMENT     :array[0..MAXSEQ,1..MAXLEN] of char;
      SEQLEN        :array[0..MAXSEQ] of integer;
      NUMSEQ:integer;

      (* Variables associated with command line options *)
      SEQFILE:text; (* sequence file, also read by ralign *)
      ARGUMENT:packed array[1..132] of char; (* command line argument *)
      FILEARGUMENT,        (*=false if argument preceeded by '-' *)
      CAPS,PERIODS:boolean;
      FORMAT:ALIGNTYPE;
      GAPS,UNKNOWN:char;
      I,
      ARGNUM,TOTALARGS,
      STARTNUM,SIGN,
      LINELENGTH:integer;
      DUMMY:LINE;     

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
    
(* BEGIN MODULE READWORD *)
  (*  Read a word from a textfile           *)
  procedure READWORD(var F:text;var W:WORD);
    var I : integer;
        CH: char;
    begin
      with W do begin
        LEN:=0;
        while F^ = ' ' do
          if not eoln(F) then read(F,CH)
          else if not eof(F) then readln(F);
        while (F^ <> ' ') do 
          if LEN < MAXWORD then begin
            LEN := LEN + 1;
            read(F,STR[LEN])
            end
          else read(F,CH);
        for I := LEN+1 to MAXWORD do STR[I]:= ' '
        end
    end; (* READWORD *)
(* END MODULE READWORD         VERSION= 'SUNMODS     Version  8/ 9/94'; *)

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

(* BEGIN MODULE TOLOWER *)
    (* Change a character from upper to lowercase *)
      function TOLOWER(CH:char):char;
       begin
       if not(CH in ['A'..'Z']) then TOLOWER:=CH
       else case CH of
        'A':TOLOWER:='a'; 'B':TOLOWER:='b'; 'C':TOLOWER:='c'; 'D':TOLOWER:='d';
        'E':TOLOWER:='e'; 'F':TOLOWER:='f'; 'G':TOLOWER:='g'; 'H':TOLOWER:='h';
        'I':TOLOWER:='i'; 'J':TOLOWER:='j'; 'K':TOLOWER:='k'; 'L':TOLOWER:='l';
        'M':TOLOWER:='m'; 'N':TOLOWER:='n'; 'O':TOLOWER:='o'; 'P':TOLOWER:='p';
        'Q':TOLOWER:='q'; 'R':TOLOWER:='r'; 'S':TOLOWER:='s'; 'T':TOLOWER:='t';
        'U':TOLOWER:='u'; 'V':TOLOWER:='v'; 'W':TOLOWER:='w'; 'X':TOLOWER:='x';
        'Y':TOLOWER:='y'; 'Z':TOLOWER:='z'
       end
      end; (* TOLOWER *)
(* END MODULE TOLOWER         VERSION= 'SUNMODS     Version  8/ 9/94'; *)

  (***************************************************************)
  (* Read aligned sequences into arrays.                         *)
  (***************************************************************)
  procedure READALIGNMENT(var INFILE,OUTFILE:text);
    var CLINE:LINE;
        LONGEST,      (* length of longest sequence *)
        I:integer;
        ENDOFSEQUENCE:boolean;
        CH:char;

    begin
      if FORMAT = RALIGN then begin (*ie. RALIGN output *)
         (* Copy consensus word lines of RALIGN output to OUTFILE. *)
      	if not eof(INFILE) then readln(INFILE); (*ignore 1st line of pralign
						 output *)
          while INFILE^ in ['0'..'9'] do begin
            READLINE(INFILE,CLINE);
            WRITELINE(OUTFILE,CLINE,CLINE.LEN);writeln(OUTFILE)
            end;
          writeln(OUTFILE);
          readln(INFILE) (* ignore consensus line *)
          end; (* RALIGN *)
        
      (* Read in each sequence. *)
      NUMSEQ:=0; LONGEST:=0;
      for I:= 0 to MAXSEQ do SEQLEN[I]:=0;

      while not eof(INFILE) do begin
        (* IG,PEARSON: move past comments to name *)
        if FORMAT in [IG,PEARSON] then begin
           while INFILE^=';' do readln(INFILE); (* ignore comment lines *)
           if not eof(INFILE) then begin
              if INFILE^='>' then read(INFILE,CH); (*read past '>' *)
              NUMSEQ:=NUMSEQ+1;              
              READWORD(INFILE,NAMES[NUMSEQ]); (* read name *)
              if not eof(INFILE) then readln(INFILE)
              end
           end (* FORMAT in [IG,PEARSON] *)
        else NUMSEQ:=NUMSEQ+1;

        (* Read sequence *)
        ENDOFSEQUENCE:= false;
        while not ENDOFSEQUENCE do begin
          while not eoln(INFILE) do begin
            read(INFILE,CH);
            (* Truncate alignment if > MAXLEN *)
            if SEQLEN[NUMSEQ] < MAXLEN then begin
               SEQLEN[NUMSEQ]:=SEQLEN[NUMSEQ]+1;
               if CH in [' ','-'] then CH:=GAPS;
               ALIGNMENT[NUMSEQ,SEQLEN[NUMSEQ]]:=TOLOWER(CH)
               end; (* if *)
            end;
          if not eof(INFILE) then readln(INFILE);
          if (eof(INFILE)) or (FORMAT=RALIGN) then ENDOFSEQUENCE:=true
          else if INFILE^ in [';','>'] then ENDOFSEQUENCE:=true;
          if ENDOFSEQUENCE then
             if SEQLEN[NUMSEQ] > LONGEST then LONGEST:= SEQLEN[NUMSEQ]
          end; (* not ENDOFSEQUENCE *)
        end; (* not eof(INFILE) *)
      SEQLEN[0]:=LONGEST
    end; (* READALIGNMENT *)

(* BEGIN MODULE TOUPPER *)
    (* Change a character from lower to uppercase *)
      function TOUPPER(CH:char):char;
       begin
       if not(CH in ['a'..'z']) then TOUPPER:=CH
       else case CH of
        'a':TOUPPER:='A'; 'b':TOUPPER:='B'; 'c':TOUPPER:='C'; 'd':TOUPPER:='D';
        'e':TOUPPER:='E'; 'f':TOUPPER:='F'; 'g':TOUPPER:='G'; 'h':TOUPPER:='H';
        'i':TOUPPER:='I'; 'j':TOUPPER:='J'; 'k':TOUPPER:='K'; 'l':TOUPPER:='L';
        'm':TOUPPER:='M'; 'n':TOUPPER:='N'; 'o':TOUPPER:='O'; 'p':TOUPPER:='P';
        'q':TOUPPER:='Q'; 'r':TOUPPER:='R'; 's':TOUPPER:='S'; 't':TOUPPER:='T';
        'u':TOUPPER:='U'; 'v':TOUPPER:='V'; 'w':TOUPPER:='W'; 'x':TOUPPER:='X';
        'y':TOUPPER:='Y'; 'z':TOUPPER:='Z'
       end
      end; (* TOUPPER *)
(* END MODULE TOUPPER         VERSION= 'SUNMODS     Version  8/ 9/94'; *)

  (***************************************************************)
  (* Modify the alignment according to the command line options. *)
  (***************************************************************)
  procedure MODIFY;
    var I,J:integer;
        CONCHAR:char;
        AGREE:boolean;
       
    (* Return the majority character position J, or if there is no majority,
       "n" for DNA and "x" for proteins *)
    function CONSENSUS(J:integer):char;
      var I,K,SETSIZE,KNOWNCOUNT,MAXNUMBER:integer;
          LETTER:array[1..MAXSEQ] of char;
          NUMBER:array[1..MAXSEQ] of integer;
      begin
        (* Count the number of each different type of character at position
           J. SETSIZE=number of different characters examined so far. *)
        SETSIZE:=1;
        LETTER[1]:=ALIGNMENT[1,J]; NUMBER[1]:=1;
        for I:= 2 to NUMSEQ do begin
          K:=1;
          while (ALIGNMENT[I,J]<>LETTER[K]) and (K<=SETSIZE) do K:=K+1;
          if K > SETSIZE then begin
               LETTER[K]:=ALIGNMENT[I,J]; NUMBER[K]:=1; SETSIZE:=SETSIZE+1 end
          else NUMBER[K]:=NUMBER[K]+1
          end; (* for I *)
        
        (* Set CONSENSUS equal to the majority character, if there is one *)
        KNOWNCOUNT:=0;
        if LETTER[1]=UNKNOWN then I:=2 else I:=1; 
        MAXNUMBER:=I;
        for K:= I to SETSIZE do
          if LETTER[K]<>UNKNOWN then begin
            if (NUMBER[K] > NUMBER[MAXNUMBER]) then MAXNUMBER:=K; 
            KNOWNCOUNT:=KNOWNCOUNT+NUMBER[K] 
            end;
        if NUMBER[MAXNUMBER] > KNOWNCOUNT div 2 
           then CONSENSUS:=LETTER[MAXNUMBER]
           else CONSENSUS:=UNKNOWN        
      end; (* CONSENSUS *) 
 
    begin (*MODIFY*)
     (* The consensus will be the length of the longest sequence. Any 
        sequence shorter than the consensus will be padded with gaps. *)
     for I:= 1 to NUMSEQ do
         for J:= SEQLEN[I]+1 to SEQLEN[0] do ALIGNMENT[I,J]:= GAPS;

     (* For each position, calculate the consensus *)
     for J:= 1 to SEQLEN[0] do begin
       CONCHAR:=CONSENSUS(J);
       ALIGNMENT[0,J]:=CONCHAR;
       AGREE:=true;
       for I:= 1 to NUMSEQ do begin
         if ALIGNMENT[I,J]=CONCHAR then
           if PERIODS and (CONCHAR<>UNKNOWN) then ALIGNMENT[I,J]:='.'
           else
         else begin
              AGREE:=false;
              if ALIGNMENT[I,J] in [' ','-'] then ALIGNMENT[I,J]:=GAPS
              end
         end; (* for I *)
       if AGREE and CAPS then 
          ALIGNMENT[0,J]:= TOUPPER(ALIGNMENT[0,J])
       end (* for J *)
   end; (* MODIFY *)
    
  (***************************************************************)
  (* Write the alignment to OUTFILE.                             *)
  (***************************************************************)
  procedure WRITEALIGNMENT(var OUTFILE:text);
    var NUCSPRINTED:array[0..MAXSEQ] of integer;
        DONE:boolean; (* true when all of each sequence has been printed*)
        I,J,
        THISLINE, (* index of last nucleotide printed on a line *)
        NUMBER:integer;
    begin
      DONE:=false;
      for I:= 0 to NUMSEQ do NUCSPRINTED[I]:=0;
      NUMBER:=STARTNUM-1; (* number used for printing *)
      while not DONE do begin
        DONE:=true;
        
        (* Write a line of numbers *)
        write(OUTFILE,' ':10);
        THISLINE:=NUCSPRINTED[0] + LINELENGTH;
        if THISLINE > SEQLEN[0] then THISLINE:=SEQLEN[0];
        for J:= 1 to (THISLINE-NUCSPRINTED[0]) div 10 do begin
          (* The coordinate 0 is not skipped. *)
          if (NUMBER>=-10) and (NUMBER<0) then NUMBER:=NUMBER+11 
          else NUMBER:=NUMBER+10;
          write(OUTFILE,NUMBER:10)
          end;
        writeln(OUTFILE);

        (* Write LINELENGTH nucleotides for each position *)
        for I:= 0 to NUMSEQ do begin
          WRITEWORD(OUTFILE,NAMES[I],10);
          THISLINE:=NUCSPRINTED[I] + LINELENGTH;
          if THISLINE > SEQLEN[I] then THISLINE:=SEQLEN[I];
          for J:= NUCSPRINTED[I]+1 to THISLINE do
            write(OUTFILE,ALIGNMENT[I,J]);
          NUCSPRINTED[I]:=THISLINE;
          writeln(OUTFILE);
          if NUCSPRINTED[I] < SEQLEN[I] then DONE:=false
          end; (* for I *)
        writeln(OUTFILE)
        end (* not DONE *)
    end; (* WRITEALIGNMENT *)
  
  
  (* ----------------------------------------------------------  *)
  (* ----------------- MAIN  PROCEDURE  -----------------------  *)
    begin
      (* Read options from command line *)
      GAPS:=' '; CAPS:=false; PERIODS:=false; UNKNOWN:='x'; FORMAT:=RALIGN;
      LINELENGTH:=70; STARTNUM:=1;
      ARGNUM:=STARTARGNUM; FILEARGUMENT:=false;
      TOTALARGS:=argc-1;

      while (ARGNUM <= TOTALARGS) and (not  FILEARGUMENT)do begin
        argv(ARGNUM,ARGUMENT);
        if ARGUMENT[1]='-' then begin
          if ARGUMENT[2] in ['g','c','n','f','p','m','l','s'] then
            case ARGUMENT[2] of
              'g':GAPS:='-';
              'p':PERIODS:=true;
              'c':CAPS:=true;
              'n':UNKNOWN:='n';
              'f':begin
                    if ARGUMENT[3] in ['r','i','p'] then
                       case ARGUMENT[3] of
                         'r':FORMAT:=RALIGN;
                         'i':FORMAT:=IG;
                         'p':FORMAT:=PEARSON
                         end (* case *)
                  end;
              'm': FORMAT:=PEARSON;
              's':begin
                STARTNUM:=0; I:=3;
                if ARGUMENT[I]='-' then begin SIGN:=-1; I:=I+1 end else SIGN:=1;
                while ARGUMENT[I] in ['0'..'9'] do begin
                    STARTNUM:= (STARTNUM*10) + ord(ARGUMENT[I])-ord('0');
                    I:=I+1
                    end;
                STARTNUM:=STARTNUM*SIGN
                end;(* s *)
              'l':begin
                  LINELENGTH:=0; I:=3;
                  while ARGUMENT[I] in ['0'..'9'] do begin
                    LINELENGTH:= (LINELENGTH*10) + ord(ARGUMENT[I])-ord('0');
                    I:=I+1
                    end;
                  if LINELENGTH <1 then LINELENGTH:=70
                  end (* l *)
              end; (* case *)
          ARGNUM:=ARGNUM+1
          end (* ARGUMENT[1]='-' *)
        else FILEARGUMENT:=true
        end; (* while *)

      (* Read sequence names from SEQFILE. (RALIGN output only) *)
      if FORMAT = RALIGN then begin
        FILEARGS(SEQFILE,'I',ARGNUM);
        I:=0; NAMES[0].LEN:=0; (* consensus name is blank *)
        while not eof(SEQFILE) do begin
          (* skip comment lines *)
          while SEQFILE^=';' do readln(SEQFILE);
          if not eof(SEQFILE) then begin
            (* read the name *)     
            I:=I+1;
            READWORD(SEQFILE,NAMES[I]);
            if not eof(SEQFILE) then readln(SEQFILE);
            (* skip the sequence *)
            repeat
              READLINE(SEQFILE,DUMMY)
            until (DUMMY.LEN=0) or (eof(SEQFILE))
            end (* if not eof(SEQFILE) *)
          end; (* while *)
        end; (* RALIGN *)

      (* Read in the aligned sequences produced by ralign.*)        
      READALIGNMENT(input,output);
      
      (* Modify the output according to command line options.*)
      MODIFY;
      
      (* Write the re-formatted alignment with numbering.*)
      WRITEALIGNMENT(output)
    end. (* REFORM  *)
 
