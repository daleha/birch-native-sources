  (* ********************************************************  *)
  (*                                                           *)
  (*  FLAT2PHYL  VERSION  5/ 2/2000, Standard Pascal           *)
  (*             Brian Fristensky                              *)
  (*             Dept. of Plant Science                        *)
  (*             University of Manitoba                        *)
  (*             Winnipeg, MB  R3T 2N2 CANADA                  *)
  (*                                                           *)
  (* SYNOPSIS                                                  *)
  (*    flat2phyl [-i] < GDE_flatfile > Phylip_file            *)
  (*                                                           *)
  (* DESCRIPTION                                               *)
  (*    Convert a GDE flatfile into a                          *)
  (*    PHYLIP discrete character file.                        *)
  (*                                                           *)
  (*             -s  write in Phylip sequential format.        *)
  (*                 rather than the default, which is         *)
  (*                 Phylip interleaved format.                *)
  (*             -i  invert characters, so that 0 -> 1 and     *)
  (*                 1 -> 0.                                   *)
  (*                                                           *)
  (* Copyright (c) 1996        by Brian Fristensky.            *)
  (* !!! in comment indicates feature which may need change.   *)
  (*  *******************************************************  *)
 
(*!!!*)program PHYL2FLAT(input, output);
 (*!!! Some Pascals require file parameters in heading *)
 
  const MAXCHAR  =  10000;
        MAXSEQ =  500;
        MAXWORD = 10;
        STARTARGNUM=1;

        VERSION = 'FLAT2PHYL     Version 5/ 2/2000';
 
  type 
(* BEGIN MODULE TYPE.WORD *)
       (*   <word>::= <non-blank char>[<non-blank char>] *)
       WORD    = record
                 LEN:integer;
                 STR:array[1..MAXWORD] of char
                 end;
(* END MODULE TYPE.WORD         VERSION= 'SUNMODS     Version  8/ 9/94'; *)

       NAMEARRAY   = array[1..MAXSEQ] of WORD;
       SEQARRAY    = array[1..MAXSEQ,1..MAXCHAR] of char;
       INTARRAY    = array[1..MAXSEQ] of integer;
      
  var SEQ          : SEQARRAY;
      NAME         : NAMEARRAY;
      LENGTH       : INTARRAY;
      NUMSEQ,
      SEQLEN       : integer;
      (* global processing parameters *)
      ARGNUM,
      POS:integer;
      SEQUENTIAL,
      INVERT,
      UNREADOPTIONS:boolean;
      ARGUMENT:packed array[1..132] of char; (* command line argument *)

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

  (***************************************************************)
  (* Read aligned sequences into arrays.                         *)
  (***************************************************************)
  procedure READFLAT(var F:text; var NUMSEQ:integer; var NAME:NAMEARRAY;
                     var SEQ:SEQARRAY; var SEQLEN:integer);
    var I,MINREAD,MAXREAD:integer;
        ENDOFSEQUENCE:boolean;
        CH:char;


    (* --------------------------------------------*)
    (* Calculate the minimum and maximum sequence lengths read so far. *)
    procedure MINMAXREAD(var LENGTH:INTARRAY; var MINREAD,MAXREAD:integer); 
      var I:integer;
      begin
        MINREAD:=LENGTH[1]; MAXREAD:=MINREAD;
        for I:= 2 to NUMSEQ do begin
            if LENGTH[I] < MINREAD then MINREAD:= LENGTH[I];
            if LENGTH[I] > MAXREAD then MAXREAD:= LENGTH[I]
            end
      end; (* MINMAXREAD *)

    begin (* READFLAT ----------------------------------------------- *)     
      (* Read in each sequence. *)
      NUMSEQ:=0;
      for I:= 1 to MAXSEQ do LENGTH[I]:=0;

      while not eof(F) do begin
        (* Read a name *)
        if F^='"' then begin
           read(F,CH); (*read past '"' *)
           NUMSEQ:=NUMSEQ+1;              
           READWORD(F,NAME[NUMSEQ]); (* read name *)
           if not eof(F) then readln(F)
           end;

        (* Read sequence *)
        ENDOFSEQUENCE:= false;
        while not ENDOFSEQUENCE do begin
          while not eoln(F) do begin
            read(F,CH);
            if CH in ['0','1','+','-','?','P','B'] then begin
              LENGTH[NUMSEQ]:=LENGTH[NUMSEQ]+1;
              SEQ[NUMSEQ,LENGTH[NUMSEQ]]:=CH
              end
            end; (* while not eoln(F) *)
          if not eof(F) then readln(F);
          if (eof(F)) then ENDOFSEQUENCE:=true
          else if F^ = '"' then ENDOFSEQUENCE:=true
          end; (* not ENDOFSEQUENCE *)
        end; (* not eof(F) *)

        MINMAXREAD(LENGTH,MINREAD,MAXREAD);
        SEQLEN:=MINREAD
    end; (* READFLAT *)


  (* ******************************************************************* *)
  (* Invert the sequence ie. 1->0 and 0->1                               *)
  (* ******************************************************************* *)
  procedure INVERTSEQ(var SEQ:SEQARRAY);

    var I,K:integer;
    
      begin 
        for I:= 1 to NUMSEQ do
            for K:= 1 to SEQLEN do 
                if SEQ[I,K] in ['0','1','+','-'] then 
                   case SEQ[I,K] of
                     '0':SEQ[I,K]:= '1';
                     '1':SEQ[I,K]:= '0';
                     '-':SEQ[I,K]:= '+';
                     '+':SEQ[I,K]:= '-'
                     end (* case *)
      end; (* INVERTSEQ *)

  (* ******************************************************************* *)
  (* Write names and sequences in Phylip interleaved format.             *)
  (* ******************************************************************* *)
   procedure WRITEPHYL(var F:text; NUMSEQ,SEQLEN:integer; NAME:NAMEARRAY; 
                       SEQ:SEQARRAY; SEQUENTIAL:boolean);
     const LINELEN = 50;
     var I,K,THISLINE,SEQPRINTED:integer;
         FIRSTLINE:boolean;
     begin
       writeln(F,NUMSEQ:10,SEQLEN:10);
       
       if SEQUENTIAL then begin (* Phylip Sequential Format*)
          for I:= 1 to NUMSEQ do begin
              (* some PHYLIP programs want 10 letter names, 
                 with at least 1 char of padding at the end. *)
              WRITEWORD(F,NAME[I],11);writeln(F);
              SEQPRINTED:=0;
              while SEQPRINTED < SEQLEN do begin
                 (* Write a line of sequence data *)
                 K:=1;
                 THISLINE:=SEQPRINTED+LINELEN;
                 if THISLINE > SEQLEN then THISLINE:=SEQLEN;
                 for K:= SEQPRINTED+1 to THISLINE do
                        write(F,SEQ[I,K]);
                 writeln(F);
                 SEQPRINTED:=THISLINE
                 end (*while*)
              end (* for *)

          end (* Phylip Sequential Format*)
          
       else begin (* Phylip Interleaved Format*)
          FIRSTLINE:= true; SEQPRINTED:=0;
          while SEQPRINTED < SEQLEN do begin
            for I:= 1 to NUMSEQ do begin

               (* Write out names on the first set of lines *)
               if FIRSTLINE then WRITEWORD(F,NAME[I],10);

               (* Write a line of sequence data *)
               K:=1;
               THISLINE:=SEQPRINTED+LINELEN;
               if THISLINE > SEQLEN then THISLINE:=SEQLEN;
               for K:= SEQPRINTED+1 to THISLINE do
                      write(F,SEQ[I,K]);
               writeln(F)
               end; (* for I:= 1 to NUMSEQ *)
             FIRSTLINE:=false; 
             SEQPRINTED:=THISLINE
           end (* while SEQPRINTED < SEQLEN *)
         end (* Phylip Interleaved Format*)
      end; (* WRITESEQ *)
 
  (* ----------------------------------------------------------  *)
  (* ----------------- MAIN  PROCEDURE  -----------------------  *)
    begin

      (* Read options from the command line. Note that -s is required *)
      INVERT:= false; SEQUENTIAL:= false;
      ARGNUM:= STARTARGNUM;
      UNREADOPTIONS:=true;
      while UNREADOPTIONS do
        if ARGNUM < argc then begin
          argv(ARGNUM,ARGUMENT);
          if ARGUMENT[1]='-' then begin
            if ARGUMENT[2] in ['i','s'] then begin
              POS:=3;
              case ARGUMENT[2] of
                'i':INVERT:=true;
                's':SEQUENTIAL:=true
                end (* case*)
              end (* i *)
            end (* '-' *)
          else UNREADOPTIONS:=false;
          ARGNUM:=ARGNUM+1
          end (* if ARGNUM <= argc *)
        else UNREADOPTIONS:=false;

      (* Read in names and sequence *)
      READFLAT(input,NUMSEQ,NAME,SEQ,SEQLEN);

      (* If -i, then invert ie. 1->0 and 0->1 *)
      if INVERT then INVERTSEQ(SEQ);

      (* Write names and sequence to output file *)
      WRITEPHYL(output,NUMSEQ,SEQLEN,NAME,SEQ,SEQUENTIAL);
   
     end. (* FLAT2PHYL *)
