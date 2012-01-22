  (* ********************************************************  *)
  (*                                                           *)
  (*   SHUFFLE   VERSION 10/23/94, Standard Pascal             *)
  (*             Brian Fristensky                              *)
  (*             Dept. of Plant Science                        *)
  (*             University of Manitoba                        *)
  (*             Winnipeg, MB  R3T 2N2 CANADA                  *)
  (*                                                           *)
  (* SYNOPSIS                                                  *)
  (*    shuffle -s<n> [-w<n> -o<n>]                            *)
  (*                                                           *)
  (* DESCRIPTION                                               *)
  (*    Shuffles nucleici acid or protein sequences in a       *)
  (*    Pearson (.wrp) format file. see Lipman et al           *)
  (*                                                           *)
  (*    -s<n>  n is a random integer between 0 and 32767.      *)
  (*           This number must be provided for each run.      *)
  (*                                                           *)
  (*    -w<n>  n is an integer, indicating the width of the    *)
  (*           window for random localization. It should       *)
  (*           never exceed the length of the shortest input   *)
  (*           sequence.                                       *)
  (*                                                           *)
  (*    -o<n>  n is an integer, indicating the number of nuc-  *)
  (*           leotides overlap between adjacent windows. It   *)
  (*           should never exceed the window size.  o def-    *)
  (*           aults to 0 of not specified.                    *)
  (*                                                           *)
  (*    If w and o are specified, overlapping windows of w     *)
  (*    nucleotides are shuffled, thus preserving the local    *)
  (*    characteristic base composition. Windows overlap by    *)
  (*    o nucleotides.                                         *)
  (*                                                           *)
  (*    If w and o are not specified, each sequence is shuf-   *)
  (*    fled globally, thus preserving the overall base        *)
  (*    composition, but not the local variations in comp.     *)
  (*                                                           *)
  (* Copyright (c) 1988, 1990  by Brian Fristensky.            *)
  (* !!! in comment indicates feature which may need change.   *)
  (*  *******************************************************  *)
 
(*!!!*)program SHUFFLE(input, output);
 (*!!! Some Pascals require file parameters in heading *)
 
  const MAXSEQ  =  750000;
        MAXWORD = 25;
        MAXLINE = 120;
        STARTARGNUM=1;
(* BEGIN MODULE REALLIMITS *)
(*!!!*)MAXREAL = 1.7E38;
(* END MODULE REALLIMITS         VERSION= 'SUNMODS     Version  8/ 9/94'; *)

        VERSION = 'SHUFFLE       Version 10/23/94';
 
  type SEQCHAR  = 
                (T,U,C,A,G,N,R,Y,M,W,S,K,D,H,V,B,L,Z,F,P,E,I,Q,ASTERISK,X,GAP);
       SEQUENCE    = array[1..MAXSEQ] of SEQCHAR;
      
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
 
  var SEQ          : SEQUENCE;
      SEQLEN       : integer;
      SCHAR       : array[SEQCHAR] of char;
      (* global processing parameters *)
      ARGNUM,
      POS:integer;
      UNREADOPTIONS:boolean;
      ARGUMENT:CHARARRAY; 
      RSEED,OVERLAP,WINDOW: integer;
      J:integer;

(* BEGIN MODULE NUMBER *)
  (* Extract an integer from a CHARARRAY (see TYPE.LINE), starting
     at the current position. *)
  function NUMBER(TARGET:CHARARRAY; var POS:integer):integer;
    var N:integer;
        ORDZERO:integer;
    begin
      ORDZERO:=ord('0');
      N:=0;
      (* evaluate characteristic *)
      while TARGET[POS] in ['0'..'9'] do begin
        N:= (N*10) + (ord(TARGET[POS])-ORDZERO);
        POS:=POS+1
        end; (* while *)
      NUMBER:=N
    end; (* NUMBER *)
(* END MODULE NUMBER         VERSION= 'SUNMODS     Version  8/ 9/94'; *)
            
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
 
(* BEGIN MODULE NUCPROT *)
  (*****************************************************************)
  (*  Convert a character to the appropriate sequence symbol.      *)
  (*****************************************************************)
  function NUCPROT(CH:char):SEQCHAR;
    begin
      case CH of
       'A','a': NUCPROT:= A;'C','c': NUCPROT:= C;'G','g': NUCPROT:= G; 
       'U','u': NUCPROT:=U; 'L','l':NUCPROT:=L;
       'T','t': NUCPROT:= T;'R','r': NUCPROT:= R;'M','m': NUCPROT:= M;
       'B','b': NUCPROT:= B;'N','n': NUCPROT:= N;'Y','y': NUCPROT:= Y;
       'K','k': NUCPROT:= K;'D','d': NUCPROT:= D;'S','s': NUCPROT:= S;
       'W','w': NUCPROT:= W;'H','h': NUCPROT:= H;'V','v': NUCPROT:= V;
       'Z','z':NUCPROT:=Z;'F','f':NUCPROT:=F;'P','p':NUCPROT:=P;
       'E','e':NUCPROT:=E;'I','i':NUCPROT:=I;'Q','q':NUCPROT:=Q;
       '*':NUCPROT:=ASTERISK;'X','x':NUCPROT:=X;'-':NUCPROT:=GAP;
       end
    end;
(* END MODULE NUCPROT *)

  (*****************************************************)
  (* Initialization procedures.                        *)
  (*****************************************************)
  procedure INITIALIZE;
    begin
      (*Set default values for parameters *)
      WINDOW:=0; OVERLAP:=0; RSEED:=0;
 
      (*   Initialize SCHAR array which holds the          *)
      (*   character values of the SEQCHARS.               *)
      SCHAR[A]:='A'; SCHAR[C]:='C';SCHAR[G]:='G'; SCHAR[T]:='T';
      SCHAR[U]:='U'; SCHAR[L]:='L';
      SCHAR[R]:='R'; SCHAR[Y]:='Y';SCHAR[M]:='M'; SCHAR[W]:='W';
      SCHAR[S]:='S'; SCHAR[K]:='K';SCHAR[D]:='D'; SCHAR[H]:='H';
      SCHAR[V]:='V'; SCHAR[B]:='B';SCHAR[N]:='N'; SCHAR[Z]:='Z';
      SCHAR[F]:='F'; SCHAR[P]:='P';SCHAR[E]:='E'; SCHAR[I]:='I';
      SCHAR[Q]:='Q'; SCHAR[X]:='X';SCHAR[ASTERISK]:='*'; SCHAR[GAP]:='-';
 
     end; (* INITIALIZE *)


(* BEGIN MODULE RANDSEED *)
{SKIP}
(* The c function srand seeds the random number generator.
   This declaration is tolerated by some C compilers, but considered
   redundant by others. Therefore, the SKIP comments cause p2c to
   omit it from the C-code. *)
(*!!!*) procedure srand(X:integer); external c;
{NOSKIP}
  procedure RANDSEED(X:integer);
    begin
      srand(X)
    end; (* RANDINT *)
(* END MODULE RANDSEED *)

(* BEGIN MODULE RANDREAL *)
{SKIP}
(* The c function rand generates a random integer.
   This declaration is tolerated by some C compilers, but considered
   redundant by others. Therefore, the SKIP comments cause p2c to
   omit it from the C-code. *)
(*!!!*) function rand:integer; external c;
{NOSKIP}

  (* It's real hard to make random number generation portable, because the 
     C rand() function generates random integers in an architecture dependent
     range, usually in a range between 0 and (2**16)-1 or (2**32)-1.*)
  function RANDREAL:real;
    var X:integer;
    begin
{SKIP}
      (* This section will only be read by Pascal. *)
      X:=rand;
{NOSKIP}

{EMBED
    X_=rand();
}
      (* if  random integer is a  31bit integer, perform an arithmetic shift
         right 16 bits to make it a 16 bit integer. *) 
      if X > 65536 then 
         X:= X div 65536;

      (* Convert random integer into random real between 0 and 1.*)
      RANDREAL:=X/32768
    end; (* RANDREAL *)
(* END MODULE RANDREAL *) 

  (* ******************************************************************* *)
  (*  Randomize each sequence by shuffling groups of WINDOW bases .      *)
  (*  Each adjecent window overlaps by OVERLAP positions.                *)
  (* ******************************************************************* *)
  procedure SCRAMBLE(var SEQFILE,OUTFILE:text);

    var COMMENTLINE:LINE;
        RNUM:real;
    (* - - - - - - - - - - - - - - - - - - - - - - - - - - -*)
    (* Read in the next sequence                            *)
     procedure READSEQUENCE(var SEQFILE:text; var SEQ:SEQUENCE);
        var CH:char;
            TOOLONG:boolean;
        begin
          SEQLEN := 0; TOOLONG:=false;
          while (not eof(SEQFILE)) and (not (SEQFILE^ in ['>',';'])) do begin
            while not eoln(SEQFILE) do begin
              read(SEQFILE,CH);
              if CH in ['A','a','C','c','G','g','T','t','U','u','N','n',
                        'R','r','Y','y','M','m','W','w','S','s','K','k',
                        'L','l',
                        'D','d','H','h','V','v','B','b','Z','z','F','f',
                        'P','p','E','e','I','i','Q','q','*','X','x','-'] then
                if SEQLEN < MAXSEQ-2 then begin
                  SEQLEN := SEQLEN + 1;
                  SEQ[SEQLEN]:=NUCPROT(CH);
                  end
                else TOOLONG:= true
              end; (* eoln *)
              readln(SEQFILE)
              end; (* eof *)
          if TOOLONG then writeln(
             '>>> WARNING! Sequence length exceeds MAXSEQ-2. Seq. truncated.')
        end; (* READSEQUENCE *)

    (* - - - - - - - - - - - - - - - - - - - - - - - - - - -*)
   procedure RANDOMIZE(var SEQ:SEQUENCE; WINDOW:integer);
 
     var LEFT,RIGHT,OFFSET:integer;

     (* Randomize nucleotides/aa's within a window by swapping. *)
     procedure SCRAMBLE_WINDOW(LEFT,RIGHT:integer);

       var I,J:integer;

       (* Swap nucleotides or amino acids between two positions. *)
       procedure SWAP(var N1,N2:SEQCHAR);
         var TEMP:SEQCHAR;
         begin
           TEMP:=N1;
           N1:=N2;
           N2:=TEMP
         end; (* SWAP *)

       begin (* SCRAMBLE_WINDOW *)
         for I:= LEFT to RIGHT do begin
           RNUM:= RANDREAL;
           J:=LEFT + round(RNUM * WINDOW);
           if J>RIGHT then J:=RIGHT; (* Make sure to stay within bounds *)
           SWAP(SEQ[I],SEQ[J])
           end
       end; (* SCRAMBLE_WINDOW *)

     begin (* RANDOMIZE *)
       LEFT:=1;
       if WINDOW < 1 then WINDOW:=SEQLEN
       else if WINDOW > SEQLEN then WINDOW:=SEQLEN;
       RIGHT:=WINDOW;
       OFFSET:=WINDOW-OVERLAP;

       while RIGHT <= SEQLEN do begin 
         SCRAMBLE_WINDOW(LEFT,RIGHT);
         LEFT:= LEFT+OFFSET;
         RIGHT:=RIGHT+OFFSET;
         if RIGHT > SEQLEN then begin (* end of sequence, a special case *)
           RIGHT:=SEQLEN;
           WINDOW:=RIGHT-LEFT+1;
           if WINDOW>1 then SCRAMBLE_WINDOW(LEFT,RIGHT);
           RIGHT:=RIGHT+OFFSET;
           end (* RIGHT > SEQLEN *)
         end (* RIGHT <= SEQLEN *)
     end; (* RANDOMIZE *)
 
   (* - - - - - - - - - - - - - - - - - - - - - - - - - - -*)
   procedure WRITESEQ(var F:text; var SEQ:SEQUENCE);
     var SEQPRINTED,THISLINE:integer;
     begin
       SEQPRINTED:=0;
       while SEQPRINTED < SEQLEN do begin
         THISLINE:=SEQPRINTED + 50;
         if THISLINE > SEQLEN then THISLINE:=SEQLEN;
         while SEQPRINTED < THISLINE do begin
           SEQPRINTED:=SEQPRINTED + 1;
           write(F,SCHAR[SEQ[SEQPRINTED]])
           end;
         writeln(F)
         end (* SEQPRINTED < SEQLEN *)
      end; (* WRITESEQ *)
     
      begin (* ----------------- SCRAMBLE ----------------------- *)
      reset(SEQFILE);
      (* The next two output lines appear as comments in the first two
         lines of the output.*)
      writeln(OUTFILE,'>',VERSION);
      writeln(OUTFILE,'> RANDOM SEED: ',RSEED:10,' ':10,'WINDOW: ',WINDOW:5,
                      ' ':10,'OVERLAP: ',OVERLAP:5);

      RANDSEED(RSEED);
      while not eof(SEQFILE) do begin
   
        (* Read in comment lines *)
        while SEQFILE^ in ['>',';'] do begin
           READLINE(SEQFILE,COMMENTLINE);
           WRITELINE(OUTFILE,COMMENTLINE,COMMENTLINE.LEN);
           writeln(OUTFILE)
           end; (* SEQFILE^ in ['>',';'] *)

        (* Process a sequence *)
        if not eof(SEQFILE) then begin
          READSEQUENCE(SEQFILE,SEQ);
          RANDOMIZE(SEQ,WINDOW);
          WRITESEQ(OUTFILE,SEQ) 
          end
        end (* not eof(SEQFILE) *)
     end; (* SCRAMBLE *)
 

  (* ----------------------------------------------------------  *)
  (* ----------------- MAIN  PROCEDURE  -----------------------  *)
    begin
      INITIALIZE;

      (* Read options from the command line. Note that -s is required *)
      ARGNUM:= STARTARGNUM;
      UNREADOPTIONS:=true;
      while UNREADOPTIONS do
        if ARGNUM < argc then begin
          argv(ARGNUM,ARGUMENT);
          if ARGUMENT[1]='-' then begin
            if ARGUMENT[2] in ['s','w','o'] then begin
              POS:=3;
              case ARGUMENT[2] of
                's':RSEED:=NUMBER(ARGUMENT,POS);
                'w':WINDOW:=NUMBER(ARGUMENT,POS);
                'o':OVERLAP:=NUMBER(ARGUMENT,POS)
                end (* case*)
              end (* s,w,o *)
            end (* '-' *)
          else UNREADOPTIONS:=false;
          ARGNUM:=ARGNUM+1
          end (* if ARGNUM <= argc *)
        else UNREADOPTIONS:=false;

      (* Scramble the sequence(s) as specified. *)      
      if RSEED > 0 then SCRAMBLE(input,output)
      else writeln(output,
                   'ERROR >>> Random seed (-s) must be set in command line.');
    end. (* SHUFFLE *)
