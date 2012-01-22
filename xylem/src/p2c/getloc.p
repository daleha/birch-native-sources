  (* 30 May 95 CRUNCHOFFSET changed to 33 *)
  (* 29 May 95 Copied FINDDATA from GETOB. Now, NAMEFILE does not have
               be sorted, although sorted files will be searched faster. *)
  (* 21 Dec 94 restores leading blanks to lines for which they have been
               crunched by splitdb -c *)
  (**************************************************************)
  (*                                                            *)
  (*  GETLOC    VERSION   5/30/95  Standard Pascal              *)
  (*            Brian Fristensky                                *)   
  (*            Dept. of Plant Science                          *)
  (*            University of Manitoba                          *) 
  (*            Winnipeg, MB R3T 2N2  CANADA                    *)
  (*                                                            *)
  (* SYNOPSIS                                                   *)
  (* getloc -a -s -f -c namefile anofile seqfile indfile outfile*)
  (*                                                            *)
  (* DESCRIPTION                                                *)
  (*  Gets GENBANK files from database and writes to OUTFILE    *)
  (*                                                            *)
  (*    -a        write annotation data only                    *)
  (*    -s        write sequence data only                      *)
  (*    -f        write each entry to a separate file.          *)
  (*              the filename consists of the locus name,      *)
  (*              followed by a 3-letter file extesion:         *)
  (*                   .gen (default, complete genbank entry    *)
  (*                   .ano (annotation only)                   *)
  (*                  .wrp (sequence only)                      *)
  (*    -c        namefile contains accession numbers           *)
  (*    -database    one of:                                    *)
  (*                 g  - GenBank (default)                     *)
  (*                 p  - PIR (NBRF)                            *)
  (*                 e  - EMBL                                  *)
  (*                 l  - LiMB                                  *)
  (*                 v  - Vecbase (VECTOR)                      *)
  (*                                                            *)
  (*    NAMEFILE- contains names of entries to get              *)
  (*    ANOFILE - contains annotation parts of entries          *)
  (*    SEQFILE - contains sequence parts of entries            *)
  (*    INDFILE - contains linenumbers for beginning of each    *)
  (*              entry in ANOFILE and SEQFILE                  *)
  (*                                                            *)
  (*  Copyright (c)  1990 by Brian Fristensky                   *)
  (*  !!! in comment indicates feature which may need change    *)
  (**************************************************************)
  program GETLOC(NAMEFILE,ANOFILE,SEQFILE,INDFILE,OUTFILE);
(*!!!  Some Pascals require file parameters in program heading *)

  const MAXWORD = 25;
        MAXLINE = 80;
        CRUNCHOFFSET = 33; (* ASCII value of CRUNCHFLAG *)
(* BEGIN MODULE STARTARGNUM *)
	STARTARGNUM=1;    (* SUN Pascal: ARG(1) is 1st command line argument*)
      (*STARTARGNUM=2;       ATT Pascal: ARG(2) is 1st command line argument*)
(* END MODULE STARTARGNUM         VERSION= 'SUNMODS     Version  8/ 9/94'; *)

  type        
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

       DBTYPE = (GB,PIR,EMBL,VECTOR,LIMB);

       LOC = record
             NAME:WORD;
             ANOLINE,
             SEQLINE:integer 
             end;
 
  var NAMEFILE,ANOFILE,SEQFILE,INDFILE,OUTFILE:text;
      SEQNAME,
      DUMMY:WORD;
      ARGUMENT:packed array[1..132] of char;
      LOCATIONS:LOC; (* locations of a locus in .ANO and .SEQ files *)
      DATABASE:DBTYPE;
      PRINTSEQ,PRINTANO,FILES,ACNO,
      FOUND:boolean;
      ARGNUM,  
      CA,CS,     (* current lines in annotation and sequence files*)
      I:integer;
      FILENAME:LINE;
      CRUNCHFLAG:char;

  (***************************************************************)
  (* Read options from command line.                             *)
  (***************************************************************)
   procedure READOPTIONS(var ARGNUM:integer;
                         var PRINTANO,PRINTSEQ,FILES,ACNO:boolean;
			 var DATABASE:DBTYPE);
     var ARGUMENT:packed array[1..132] of char;
     begin (* READOPTIONS *)
       (* Set defaults *)
       PRINTSEQ:=true; PRINTANO:=true; FILES:=false; ACNO:=false;DATABASE:=GB;

       (* Read options.*)
       ARGNUM:=STARTARGNUM;
       repeat
         argv(ARGNUM,ARGUMENT);
         if ARGUMENT[1]='-' then begin
           if ARGUMENT[2] in ['a','s','f','c','g','p','e','l','v'] then 
             case ARGUMENT[2] of 
              'a':PRINTSEQ:=false;
              's':PRINTANO:=false;
              'f':FILES:=true;
              'c':ACNO:=true;
              'g':DATABASE:=GB;
              'p':DATABASE:=PIR;
              'e':DATABASE:=EMBL;
              'l':DATABASE:=LIMB;
              'v':DATABASE:=VECTOR
              end; (*case*)
           ARGNUM:=ARGNUM+1
           end
        until ARGUMENT[1] <> '-'
     end; (* READOPTIONS *)      
 
  (***************************************************************)
  (* Open files.                                                 *)
  (***************************************************************)
   procedure OPENFILES(var NAMEFILE,ANOFILE,SEQFILE,INDFILE,OUTFILE:text);

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
       FILEARGS(NAMEFILE,'I',ARGNUM);
       if PRINTANO then FILEARGS(ANOFILE,'I',ARGNUM);
       if PRINTSEQ then FILEARGS(SEQFILE,'I',ARGNUM);
       FILEARGS(INDFILE,'I',ARGNUM);
       if not FILES then FILEARGS(OUTFILE,'O',ARGNUM);
     end; (* OPENFILES *)

  (***************************************************************)
  (* I/O ROUTINES                                                *)
  (***************************************************************)

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

  (**************************************************************)
  (* Read a LOCUS name for next sequence to pull from database. *)
  (**************************************************************)
  procedure READNAME(var NAMEFILE:text; var SEQNAME:WORD);
    begin
      (* Skip comment lines *)
      while (NAMEFILE^=';') and (not eof(NAMEFILE)) do readln(NAMEFILE);

      (* Read the first token in the next non-comment line *)
      SEQNAME.LEN:=0;
      if not eof(NAMEFILE) then READWORD(NAMEFILE,SEQNAME);
      if not eof(NAMEFILE) then readln(NAMEFILE)
    end; (* READNAME *)

  (**********************************************************************)
  (* Look for a given LOCUS in the index file.  Since names in NAMEFILE *)
  (* are supposed to be sorted, reaching the end of file means that the *)
  (* locus in question is not in the database.                          *)
  (**********************************************************************)
  procedure FINDDATA(var INDFILE:text; ID:WORD; var LOCATIONS:LOC;
                     var FOUND:boolean);
    var FLAG:char;
        DUMMY:WORD;

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

   (* Advance to first line with first char of ID *)
   procedure ZOOM;
     begin
       while (not eof(INDFILE)) and (INDFILE^<>FLAG) do readln(INDFILE)
     end; (* ZOOM *)

   (* Search for LOCUS corresponding to ID *) 
   procedure LCRAWL;
     begin
       while (INDFILE^=FLAG) and (not FOUND) do
         with LOCATIONS do begin
           READWORD(INDFILE,NAME);
           if SAMEWORD(ID,NAME) then begin
             FOUND:=true;
             READWORD(INDFILE,DUMMY); (* ignore ACCESSION number *)
             read(INDFILE,ANOLINE,SEQLINE)
             end;
           if not eof(INDFILE) then readln(INDFILE)
           end (* with *)
      end; (* LCRAWL *)

   (* Search for ACCESSION number corresponding to ID *) 
   procedure ACRAWL;
     begin
       while (not FOUND) and (not eof(INDFILE)) do
         with LOCATIONS do begin
           READWORD(INDFILE,DUMMY); (* skip over LOCUS name *)
           READWORD(INDFILE,NAME);
           if SAMEWORD(ID,NAME) then begin
             FOUND:=true;
             read(INDFILE,ANOLINE,SEQLINE)
             end;
           if not eof(INDFILE) then readln(INDFILE)
           end (* with *)
      end; (* ACRAWL *)

   begin (* FINDDATA *)
     if eof(INDFILE) then reset(INDFILE);
     FOUND:=false;
     (* Make sure that the name or accession number is capitalized. *)
     with ID do for I:= 1 to LEN do STR[I]:= TOUPPER(STR[I]);
     FLAG:=ID.STR[1];

     (* Search to end of file or until FOUND *)
     if ACNO then ACRAWL
     else while (not FOUND) and (not eof(INDFILE)) do begin ZOOM; LCRAWL end;

     (* If list isn't sorted, go back to the beginning and search from top *)
     if not FOUND then begin 
       reset(INDFILE);
       if ACNO then ACRAWL
       else while (not FOUND) and (not eof(INDFILE)) do begin ZOOM; LCRAWL end;
       end
  end; (* FINDDATA *)

  (*****************************************************)
  (* Create a filename using the locus.                *)
  (*****************************************************)
  procedure MAKEFN(NAME:WORD; var FILENAME:LINE);
    var I:integer;
    begin
      with FILENAME do begin
        I:=1;
        while I <= NAME.LEN do begin
          STR[I]:= NAME.STR[I];
          I:=I+1
          end;
        STR[I]:='.';
        if not PRINTSEQ and PRINTANO then begin
          STR[I+1]:='a'; STR[I+2]:='n'; STR[I+3]:='o' end
        else if not PRINTANO and PRINTSEQ then begin
          STR[I+1]:='w'; STR[I+2]:='r'; STR[I+3]:='p' end
        else begin
          STR[I+1]:='g'; STR[I+2]:='e'; STR[I+3]:='n' end;
        LEN:=I+3;
        for I:= LEN+1 to MAXLINE do STR[I]:=' '
        end
    end; (* MAKEFN *)


  (*********************************************************)
  (* Read annoation data from ANOFILE and write to OUTFILE.*)
  (*********************************************************)
  procedure GETANO(ANOLINE:integer; var CA:integer);
    var CH:char;
    
    begin
      (* Advance to ANOLINE, the beginning of the entry.*)
      if CA > ANOLINE then begin reset(ANOFILE); CA:=1 end;
      while CA < ANOLINE do begin
        readln(ANOFILE);
        CA:=CA+1
        end;

      (* Transfer data until the end of entry marker is reached.*)
      while not(ANOFILE^='/') do begin

        (* If the second character holds a compression value. Expand
           the leading blanks so that there are 
           ord(CH)-CRUNCHOFFSET leading blanks. *)
        if ANOFILE^=CRUNCHFLAG then begin
           read(ANOFILE,CH);
           read(ANOFILE,CH);
           write(OUTFILE,' ':(ord(CH)-CRUNCHOFFSET))
           end; 

        while not eoln(ANOFILE) do begin
          read(ANOFILE,CH);
          write(OUTFILE,CH)
          end; (* not eoln(ANOFILE) *)
        readln(ANOFILE);
        CA:=CA+1;
        writeln(OUTFILE)
        end; (* not(ANOFILE^='/') *)

      (* if -a option is specified, divide entries with a line for easy reading. *)
      if not PRINTSEQ then writeln(OUTFILE,
	  '//----------------------------------------------------------------------');
      readln(ANOFILE);
      CA:=CA+1
    end; (* GETANO *)

  (********************************************************)
  (* Read sequence data from SEQFILE and write to OUTFILE.*)
  (********************************************************)
  procedure GETSEQ(SEQLINE:integer; var CS:integer);
    var CH:char;
        MORESEQ:boolean; (* =true if more sequence remains to be read. *)
        SEQLEN,      (* length of sequence printed so far *)
        NUMWIDTH,    (* width of numbers printed to left of sequence *)
        SPACERWIDTH, (* width of white space between number and sequence *)
        GROUP,LINESIZE, (* print <= LINESIZE nucleotides per line in groups 
                           of GROUP *)
        I,              (* width of current group *)
	WIDTH:integer; (* width of current line *)
        NAMELINE:LINE;
    begin
      (* Advance to SEQLINE, the beginning of the entry.*)
      if CS > SEQLINE then begin reset(SEQFILE); CS:=1 end;
      while CS < SEQLINE do begin
        readln(SEQFILE);
        CS:=CS+1
        end;
      (* Read past name line *)
      READLINE(SEQFILE,NAMELINE);
      CS:=CS+1;

      (* Transfer data until the next sequence or end of file is reached.*)
      if PRINTANO then (* write in database-specific format *)
         case DATABASE of
           GB,EMBL,VECTOR:
              begin
              (* Set printing parameters according to DATABASE *)
              case DATABASE of 
                  GB: begin NUMWIDTH:=9; SPACERWIDTH:=1; LINESIZE:=60 end;
                  EMBL:begin SPACERWIDTH:=5; LINESIZE:=60 end;
                  VECTOR: begin NUMWIDTH:=8; SPACERWIDTH:=2; LINESIZE:=50 end
                  end; (* case *)
              GROUP:=10;SEQLEN:=1;

              if (not eof(SEQFILE)) and (SEQFILE^<>'>') then begin
                 MORESEQ:=true;
                 if DATABASE<>EMBL then write(OUTFILE,SEQLEN:NUMWIDTH);
                 write(OUTFILE,' ':SPACERWIDTH)
                 end
              else MORESEQ:=false;

              (* For each character in the sequence: *)
              WIDTH:=0; I:=0;
              while MORESEQ do begin
                (* Read a character, write to OUTFILE, and determine if there
                   if more sequence left to read. *)
                read(SEQFILE,CH);
                write(OUTFILE,CH);
                if eoln(SEQFILE) and (not eof(SEQFILE)) then begin
                  readln(SEQFILE);
                  CS:=CS+1
                  end;
                if eof(SEQFILE) or (SEQFILE^='>') then MORESEQ:=false;

                if MORESEQ then begin (*write spaces or /newline, if necessary*)
                  I:=I+1;
                  WIDTH:=WIDTH+1;
                  if I=GROUP then begin
                     I:=0;
                     if WIDTH<LINESIZE then write(OUTFILE,' ')
                     else begin
                       writeln(OUTFILE);
                       SEQLEN:=SEQLEN+WIDTH;
                       WIDTH:=0;
                       if DATABASE<>EMBL then write(OUTFILE,SEQLEN:NUMWIDTH);
                       write(OUTFILE,' ':SPACERWIDTH)
                       end (* WIDTH=LINESIZE *)
                     end (* I=GROUP*)
                   end (* if MORESEQ *)
                 else writeln(OUTFILE)
                end; (* while MORESEQ *)
              writeln(OUTFILE,'//')
              end; (* GB,EMBL,VECTOR *)
           PIR:begin 
              (* write number header *)
              write(OUTFILE,' ':7);
              I:=5;
              while I<=30 do begin write(OUTFILE,I:10); I:=I+5 end;
              writeln(OUTFILE);

              (* write sequence in PIR format *)
              (* punctuation characters occupy odd numbered columnes 1..61
                 and amino acids occupy even numbered columns 2..60. This
                 algorithm assumes that punctuation characters may not occur
                 adjacent to one another. *)
              NUMWIDTH:=7; LINESIZE:=60;SPACERWIDTH:=1;
              SEQLEN:=1; WIDTH:=0;

              if (not eof(SEQFILE)) and (SEQFILE^<>'>') then begin
                 MORESEQ:=true;
                 write(OUTFILE,SEQLEN:NUMWIDTH)
                 end
              else MORESEQ:=false;

              while MORESEQ do begin
                read(SEQFILE,CH);
                if eoln(SEQFILE) and (not eof(SEQFILE)) then begin
                  readln(SEQFILE);
                  CS:=CS+1;
                  end;
               if eof(SEQFILE) or (SEQFILE^='>') then MORESEQ:=false;
               if CH in ['(',')','=','/','.',','] then begin
                   write(OUTFILE,CH);
                   WIDTH:=WIDTH+1
                   end
                else begin
                   if not odd(WIDTH) then begin
                      write(OUTFILE,' '); WIDTH:=WIDTH+1 end;
                   write(OUTFILE,CH);
                   WIDTH:=WIDTH+1;
                   end;
                if MORESEQ then 
                  if WIDTH>=LINESIZE then begin
                     writeln(OUTFILE);
                     SEQLEN:=SEQLEN+30;
                     WIDTH:=0;
                     write(OUTFILE,SEQLEN:NUMWIDTH);
                     end (* WIDTH>=LINESIZE *)
                   else
                 else writeln(OUTFILE)
                 end; (* while MORESEQ *)
               writeln(OUTFILE,'///')           
               end; (*PIR*)
           LIMB:;
           end (* case DATABASE *)

      else begin (* write in Pearson format *)
        WRITELINE(OUTFILE,NAMELINE,NAMELINE.LEN);writeln(OUTFILE);
        while not((SEQFILE^='>') or (eof(SEQFILE))) do begin
          while not eoln(SEQFILE) do begin
            read(SEQFILE,CH);
            write(OUTFILE,CH)
            end;
          readln(SEQFILE);
          CS:=CS+1;
          writeln(OUTFILE)
          end;
        end  
    end; (* GETSEQ *)

  (* ----------------------------------------------------------  *)
  (* ----------------- MAIN  PROCEDURE  -----------------------  *)
    begin
      
       (* Read options from the command line *)
       READOPTIONS(ARGNUM,PRINTANO,PRINTSEQ,FILES,ACNO,DATABASE);

       (* Open files.*)
       OPENFILES(NAMEFILE,ANOFILE,SEQFILE,INDFILE,OUTFILE);
     
       
       (* For each locus in NAMEFILE, find the locations of the annotation and
          sequence parts in ANOFILE and SEQFILE, and then get the data
          and write to OUTFILE in GenBank tape format. *)
       CRUNCHFLAG:=chr(CRUNCHOFFSET);
       CA:=1; CS:=1;
       while not eof(NAMEFILE) do begin
         READNAME(NAMEFILE,SEQNAME);
         if SEQNAME.LEN > 0 then FINDDATA(INDFILE,SEQNAME,LOCATIONS,FOUND);
         if FOUND then begin
           if FILES then begin
              MAKEFN(SEQNAME,FILENAME);         
              rewrite(OUTFILE,FILENAME.STR)
              end;
           if PRINTANO then GETANO(LOCATIONS.ANOLINE,CA);
           if PRINTSEQ then GETSEQ(LOCATIONS.SEQLINE,CS);
(*!!!      if FILES then CLOSE(OUTFILE) *)
           end
         end; (* while not eof(NAMEFILE) *)
         
(*!!! if not FILES then CLOSE(OUTFILE); *)
    end. (* GETLOC  *) 
 
