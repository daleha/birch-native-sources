  (* ********************************************************  *)
  (*                                                           *)
  (*  PHYL2FLAT  VERSION  3/12/97, Standard Pascal             *)
  (*             Brian Fristensky                              *)
  (*             Dept. of Plant Science                        *)
  (*             University of Manitoba                        *)
  (*             Winnipeg, MB  R3T 2N2 CANADA                  *)
  (*                                                           *)
  (* SYNOPSIS                                                  *)
  (*    phyl2flat [-i] < Phylip_file > GDE_flatfile            *)
  (*                                                           *)
  (* DESCRIPTION                                               *)
  (*    Convert a phylip discrete character file into a        *)
  (*    GDE flatfile.                                          *)
  (*                                                           *)
  (*             -i  invert characters, so that 0 -> 1 and     *)
  (*                 1 -> 0.                                   *)
  (*             -n  read non-interleaved format               *)
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
(* BEGIN MODULE REALLIMITS *)
(*!!!*)MAXREAL = 1.7E38;
(* END MODULE REALLIMITS         VERSION= 'SUNMODS     Version  8/ 9/94'; *)

        VERSION = 'PHYL2FLAT     Version 3/12/97';
 
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
      OKAY:boolean;
      INTERLEAVED,
      INVERT,
      UNREADOPTIONS:boolean;
      ARGUMENT:packed array[1..132] of char; (* command line argument *)

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

  (* ******************************************************************* *)
  (* Read in the  names and sequences.                                   *)
  (* ******************************************************************* *)
  procedure NAMESANDSEQ(var F:text; var NAME:NAMEARRAY; var SEQ:SEQARRAY;
                        var SEQLEN:integer; INTERLEAVED:boolean;
                        var OKAY:boolean);

    var I,MINREAD,MAXREAD:integer;

    (* --------------------------------------------*)
    procedure READNAME(var F:text; var NAME:NAMEARRAY; I:integer; 
                       var OKAY:boolean); 
      var J:integer;
          CH:char;
      begin
        J:=1;
        while (J<=MAXWORD) and OKAY do
           if not eoln(F) then begin read(F,CH);
                                     if CH <> ' ' then begin
                                        NAME[I].STR[J]:=CH;
                                        NAME[I].LEN:=NAME[I].LEN+1
                                        end; (* CH <> ' ' *)
                                     J:= J+1
                                     end
            else begin
                 OKAY:=false;
                 writeln('PHYL2FLAT: >>>Truncated sequence line.')
                 end
      end; (* READNAME *)
    
    (* --------------------------------------------*)
    procedure READSEQLINE(var F:text; var SEQ:SEQARRAY; var LENGTH:INTARRAY;
                         I :integer); 
      var K:integer;
          CH:char;
      begin
        (* Read in the sequence *)
        K:=LENGTH[I];
        while (not eoln(F) and (K<MAXCHAR)) do 
              if F^ in ['0','1','+','-','?','P','B']
                     then begin K:=K+1; read(F,SEQ[I,K]) end
                     else read(F,CH);
        if not eof(F) then readln(F);
        LENGTH[I]:= K
      end; (* READSEQLINE *)

    (* --------------------------------------------*)
    (* Calculate the minimum and maximum sequence lengths read so far. *)
    procedure MINMAXREAD(LENGTH:INTARRAY; var MINREAD,MAXREAD:integer); 
      var I:integer;
      begin
        MINREAD:=LENGTH[1]; MAXREAD:=MINREAD;
        for I:= 2 to NUMSEQ do begin
            if LENGTH[I] < MINREAD then MINREAD:= LENGTH[I];
            if LENGTH[I] > MAXREAD then MAXREAD:= LENGTH[I]
            end
      end; (* MINMAXREAD *)

      begin (* NAMESANDSEQ -------------------------------------- *)
      reset(F);

      (* Read line 1, telling number of 'sequences' and number of
         characters per sequence *)
      readln(F,NUMSEQ,SEQLEN);
      for I:= 1 to NUMSEQ do begin
          NAME[I].LEN:= 0;
          LENGTH[I]:=0
          end;

      (* Read in the names and sequences *)
      OKAY:=true;
      MINREAD:=0; MAXREAD:=0;

      if INTERLEAVED then begin
         while (MAXREAD<SEQLEN) and OKAY do begin
           I:=1;
           while (I<= NUMSEQ) and OKAY do begin
             (* First set of seq. lines, read in the name.
                In interleaved format, as many additional sets of 
                sequence lines, in the same order as the first,
                will follow, without names. *)
             if NAME[I].LEN = 0 then READNAME(F,NAME,I,OKAY);
             READSEQLINE(F,SEQ,LENGTH,I);
             I:= I + 1
             end; (* while I <= NUMSEQ *)
           MINMAXREAD(LENGTH,MINREAD,MAXREAD)
           end; (* while MAXREAD < SEQLEN *)
   
         if MINREAD < MAXREAD then begin
            writeln('>>> PHYL2FLAT: sequences not all same length');
            writeln('>>> truncating to ',MINREAD);
            SEQLEN:=MINREAD
            end 
         end (* if INTERLEAVED *)
      else begin (* non-interleaved format *)
           I:=1;
           while ((I<= NUMSEQ) and OKAY) do begin
             (* First set of seq. lines, read in the name.
                In interleaved format, as many additional sets of 
                sequence lines, in the same order as the first,
                will follow, without names. *)
             if NAME[I].LEN = 0 then READNAME(F,NAME,I,OKAY);
             while (LENGTH[I] < SEQLEN)  and (not eof(F)) do
                   READSEQLINE(F,SEQ,LENGTH,I);
             I:= I + 1
             end; (* while (I <= NUMSEQ) and OKAY *)
         end (* non-interleaved format *)
      
    end; (* NAMESANDSEQ *)
 
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
  (* Write names and sequences in GDE flatfile format.                   *)
  (* ******************************************************************* *)
   procedure WRITEFLAT(var F:text; NUMSEQ,SEQLEN:integer; NAME:NAMEARRAY; 
                       SEQ:SEQARRAY);
     const LINELEN = 50;
     var I,K,THISLINE:integer;
     begin
       for I:=1 to NUMSEQ do begin
           write(F,'"'); (* indicates text in GDE flatfile *)
           WRITEWORD(F,NAME[I],NAME[I].LEN);writeln(F);
           K:=1; THISLINE:=0;
           while K<=SEQLEN do begin
             THISLINE:=THISLINE+LINELEN;
             if THISLINE > SEQLEN then THISLINE:=SEQLEN;
             while K <= THISLINE do begin
                   write(F,SEQ[I,K]);
                   K:=K+1
                   end; (* K <= THISLINE *)
             writeln(F) 
             end (* K <= SEQLEN*)        
        end (* for I *)
      end; (* WRITESEQ *)
 
  (* ----------------------------------------------------------  *)
  (* ----------------- MAIN  PROCEDURE  -----------------------  *)
    begin

      (* Read options from the command line. Note that -s is required *)
      INVERT:= false; INTERLEAVED:=true;
      ARGNUM:= STARTARGNUM;
      UNREADOPTIONS:=true;
      while UNREADOPTIONS do
        if ARGNUM < argc then begin
          argv(ARGNUM,ARGUMENT);
          if ARGUMENT[1]='-' then begin
            if ARGUMENT[2] in ['i','n'] then begin
              POS:=3;
              case ARGUMENT[2] of
                'i':INVERT:=true;
                'n':INTERLEAVED:=false
                end (* case*)
              end (* i,n *)
            end (* '-' *)
          else UNREADOPTIONS:=false;
          ARGNUM:=ARGNUM+1
          end (* if ARGNUM <= argc *)
        else UNREADOPTIONS:=false;

      (* Read in names and sequence *)
      NAMESANDSEQ(input,NAME,SEQ,SEQLEN,INTERLEAVED,OKAY);

      if OKAY then begin
          (* If -i, then invert ie. 1->0 and 0->1 *)
          if INVERT then INVERTSEQ(SEQ);

          (* Write names and sequence to output file *)
          WRITEFLAT(output,NUMSEQ,SEQLEN,NAME,SEQ);
         end (* if OKAY *)
      
 
    end. (* PHYL2FLAT *)
