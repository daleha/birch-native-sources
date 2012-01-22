  (***********************************************************)
  (*                                                         *)
  (*  IDENTIFY  VERSION   8/19/90  UNIX Pascal               *)
  (*            Brian Fristensky                             *)
  (*            Dept. of Plant Science                       *)
  (*            University of Manitoba                       *)
  (*            Winnipeg, MB  R3T 2N2  CANADA                *)
  (*                                                         *)
  (*  Uses GREP output to identify locus names in .IND file  *)
  (*  that correspond to lines found in .ANO file.           *)
  (*                                                         *)
  (*    GREPFILE- contains lines found in .ANO file by GREP  *)
  (*    INDFILE - contains linenumbers for beginning of each *)
  (*              entry in ANOFILE and SEQFILE               *)
  (*    NAMEFILE- output file listing all loci found         *)
  (*    FINDFILE- output file with locus names and lines     *)
  (*              found for each locus                       *)
  (*                                                         *)
  (*  Copyright (c) l988, 1990 by Brian Fristensky           *)
  (*  !!! in comment indicates feature which may need change *)
  (***********************************************************)
  program IDENTIFY(GREPFILE,INDFILE,NAMEFILE,FINDFILE);
(*!!!  Some Pascals require file parameters in program heading *)

  const MAXWORD = 25;
        MAXLINE = 132;
(*BEGIN MODULE STARTARGNUM *)
        STARTARGNUM=1;
(*END MODULE STARTARGNUM *)

  type        
(* BEGIN MODULE TYPE.WORD *)
       (*   <word>::= <non-blank char>[<non-blank char>] *)
       WORD    = record
                 LEN:integer;
                 STR:array[1..MAXWORD] of char
                 end;
(* END MODULE TYPE.WORD         VERSION= 'SUNMODS     Version  9/12/91'; *)

(* BEGIN MODULE TYPE.LINE *)
       CHARARRAY = packed array[1..MAXLINE] of char;
       LINE = record
                STR:CHARARRAY;
                LEN:integer
                end;
(* END MODULE TYPE.LINE         VERSION= 'SUNMODS     Version  9/12/91'; *)

  var GREPFILE,INDFILE,NAMEFILE,FINDFILE:text;
      GREPLINE:LINE;
      CURNAME,NEXTNAME,DUMMY:WORD;
      ARGNUM,GREPNUM,CURNUM,NEXTNUM:integer;
            
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
(* END MODULE FILEARGS         VERSION= 'SUNMODS     Version  9/12/91'; *)

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
(* END MODULE READLINE         VERSION= 'SUNMODS     Version  9/12/91'; *)

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
(* END MODULE WRITELINE         VERSION= 'SUNMODS     Version  9/12/91'; *)

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
(* END MODULE READWORD         VERSION= 'SUNMODS     Version  9/12/91'; *)

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
(* END MODULE WRITEWORD         VERSION= 'SUNMODS     Version  9/12/91'; *)

  (*********************************************************)
  (* Read a line number and line from GREPFILE.            *)
  (*********************************************************)
  procedure GLINE(var GREPFILE:text; var GREPLINE:LINE; var GREPNUM:integer);
    var CH:char;
        ORDZERO:integer;
    begin
      ORDZERO:=ord('0');
      GREPNUM:=0;   
      repeat
        read(GREPFILE,CH);
        if CH <> ':' then GREPNUM:=(GREPNUM*10) + (ord(CH)-ORDZERO)          
      until CH=':';
      READLINE(GREPFILE,GREPLINE)
    end; (* GLINE *)
    
  (*********************************************************)
  (* Read INDFILE until CURNUM <= GREPNUM < NEXTNUM.       *)
  (*********************************************************)
  procedure SCAN(var INDFILE:text; var CURNAME,NEXTNAME:WORD;
                 var CURNUM,NEXTNUM:integer);
    var DUMMY:WORD;
    begin
      if not eof(INDFILE) then
        while (NEXTNUM <= GREPNUM) and (not eof(INDFILE)) do begin
          CURNAME:=NEXTNAME;
          CURNUM:=NEXTNUM;
          READWORD(INDFILE,NEXTNAME);
          READWORD(INDFILE,DUMMY); (* ignore ACCESSION number *)
          readln(INDFILE,NEXTNUM)
          end
      else begin (* Special case for last entry in INDFILE *)
        CURNAME:=NEXTNAME;
        CURNUM:=NEXTNUM;
        NEXTNUM:=10000000 (* GREPNUM can never be this large *)
        end
    end; (* SCAN *)
    
                     
  (* ----------------------------------------------------------  *)
  (* ----------------- MAIN  PROCEDURE  -----------------------  *)
    begin
       (* Open files.*)
       ARGNUM:=STARTARGNUM; 
       FILEARGS(GREPFILE,'I',ARGNUM);
       FILEARGS(INDFILE,'I',ARGNUM);
       FILEARGS(NAMEFILE,'O',ARGNUM);
       FILEARGS(FINDFILE,'O',ARGNUM);

       (* Read past comment lines at beginning of INDFILE, if any.*)
       (* Read the first two loci in INDFILE *)
       while (INDFILE^=';') and (not eof(INDFILE)) do readln(INDFILE); 
       READWORD(INDFILE,CURNAME);
       READWORD(INDFILE,DUMMY); (* ignore ACCESSION number *)
       readln(INDFILE,CURNUM);
       READWORD(INDFILE,NEXTNAME);
       READWORD(INDFILE,DUMMY); (* ignore ACCESSION number *)
       readln(INDFILE,NEXTNUM);
       
       (* Read the first line in GREPFILE *)
       GLINE(GREPFILE,GREPLINE,GREPNUM);

       (* Find the LOCUS corresponding to GREPNUM *)
       SCAN(INDFILE,CURNAME,NEXTNAME,CURNUM,NEXTNUM);
       
       (* Write the first LOCUS name to NAMEFILE and FINDFILE and the
          GREPLINE to FINDFILE *)
       WRITEWORD(FINDFILE,CURNAME,CURNAME.LEN);writeln(FINDFILE);
       WRITEWORD(NAMEFILE,CURNAME,CURNAME.LEN);writeln(NAMEFILE);
       WRITELINE(FINDFILE,GREPLINE,GREPLINE.LEN);
       writeln(FINDFILE);            


       (* MAIN LOOP *)             
       (* Invariant: CURNUM <= GREPNUM < NEXTNUM *)
       while not eof(GREPFILE) do begin

         (* Read a line from GREPFILE.  The first characters on the line tell
            the line number, the rest is from the .ANO file. *)
         GLINE(GREPFILE,GREPLINE,GREPNUM);
         
         (* Search INDFILE until the position of the next entry is greater
            than the line number from the GREPFILE.*)
         if GREPNUM >= NEXTNUM then begin
           writeln(FINDFILE,'//');
           SCAN(INDFILE,CURNAME,NEXTNAME,CURNUM,NEXTNUM);
           WRITEWORD(FINDFILE,CURNAME,CURNAME.LEN);writeln(FINDFILE);
           WRITEWORD(NAMEFILE,CURNAME,CURNAME.LEN);writeln(NAMEFILE)
           end;
            
         (* Write the line into FINDFILE *)
         WRITELINE(FINDFILE,GREPLINE,GREPLINE.LEN);
         writeln(FINDFILE)            
         end;

       writeln(FINDFILE,'//');
(*!!!  CLOSE(NAMEFILE); *)
(*!!!  CLOSE(FINDFILE)  *)
    end. (* IDENTIFY  *) 
