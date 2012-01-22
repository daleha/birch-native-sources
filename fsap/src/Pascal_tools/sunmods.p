  (* ********************************************************  *)
  (*                                                           *)
  (* FRISTENSKY SEQUENCE ANALYSIS PACKAGE SUBROUTINE MODULES   *)
  (*             VERSION  6/26/01  SUN    Pascal               *)
  (*             Brian Fristensky                              *)
  (*             Dept. of Plant Science                        *)
  (*             University of Manitoba                        *)
  (*             Winnipeg, MB R3T 2N2 CANADA                   *)
  (*                                                           *)
  (* Copyright (c) 1986-2001               by Brian Fristensky *)
  (* !!! in comment indicates feature which may need change.   *)
  (*  *******************************************************  *)
  
(* --------  REVISION HISTORY --------------------------------
26 Jun 01 - Readseq updated to accommodate changes to the
GenBank LOCUS line format, which should take effect with
GenBank 126.0. The new format allows names up to 18 char. long,
and sequences up to 9,999,999 bases. One side effect is that
the topology has moved from columns 43-52 to columns 56-63.
------------------------------------------------------------*)

  program CSAPMOD(input, output (*,INFILE*));
(*!!! Some Pascals require file parameters in program heading *)
 
  const MAXSEQ  = 750000;
(* BEGIN MODULE VERSION *)
        VERSION= 'SUNMODS     Version  6/26/01';
(* END MODULE VERSION *)
 
(* BEGIN MODULE REALLIMITS *)
(*!!!*)MAXREAL = 1.7E38;
(* END MODULE REALLIMITS *)
 
(* BEGIN MODULE INTLIMITS *)
(*!!!  MAXINT =  2147483647; *)
(*!!!  MININT = -2147483647; *)
(* END MODULE INTLIMITS *)
 
(* BEGIN MODULE STARTARGNUM *)
	STARTARGNUM=1;    (* SUN Pascal: ARG(1) is 1st command line argument*)
      (*STARTARGNUM=2;       ATT Pascal: ARG(2) is 1st command line argument*)
(* END MODULE STARTARGNUM *)

        MAXWORD = 25;
        MAXLINE = 70;
        MAXRANGE=  0;

  type NUCLEOTIDE  = (CANT,A,C,R,D,V,M,K,B,H,Y,G,T,W,S,N,Z,WONT);
       AMINOACID   = (GLY,ALA,VALINE,LEU,ILE,MET,PHE,PRO,SER,THR,CYS,ASN,GLN,
                      TYR,TRP,ASP,GLU,HIS,LYS,ARG,ASX,GLX,TERM,UNKX,UNKY);
       SEQUENCE    = array[1..MAXSEQ] of NUCLEOTIDE;
       PROTEIN     = array[1..3000] of AMINOACID;
       CHSET       = set of 'A'..'Z';
       LETTERS = packed array[1..10] of char;
       VECTOR  = array[1..MAXLINE] of integer;
 
(* BEGIN MODULE TYPE.WORD *)
       (*   <word>::= <non-blank char>[<non-blank char>] *)
       WORD    = record
                 LEN:integer;
                 STR:array[1..MAXWORD] of char
                 end;
(* END MODULE TYPE.WORD         VERSION= 'CSAPMODS     Version  4/13/88'; *)
 
(* BEGIN MODULE TYPE.LINE *)
       CHARARRAY = packed array[1..MAXLINE] of char;
       LINE = record
                STR:CHARARRAY;
                LEN:integer
                end;
(* END MODULE TYPE.LINE         VERSION= 'SUNMODS     Version  3/21/90'; *)
 
(* BEGIN MODULE TYPE.FRAG *)
       (* LIST OF RESTRICTION SITES FOUND *)
       FRAG = ^FRAGMENT;
       FRAGMENT = record
         START,FINISH,SIZE:integer;
         PREV,NEXT:FRAG
         end;
       FRAGSFOUND = record
         LNUM:integer;
         HEAD,TAIL:FRAG
         end;
(* END MODULE TYPE.FRAG *)
 
  var INFILE      : text;
      CH:char;
      WO : WORD;
      LI : LINE;
      RE:real;
      SEQ         : SEQUENCE;
      SEQLEN,INT  : integer;
      CIRCULAR    : boolean;
      FREEFRAG:FRAG;
      ORDER: array[1..10] of FRAG;
 
(************ DUMMY LINES FOR MODULE STRIPPING **********)
(* BEGIN MODULE DUMMY *)
(* BEGIN MODULE REALLIMITS *)
(* END MODULE REALLIMITS *)
 
(* BEGIN MODULE INTLIMITS *)
(* END MODULE INTLIMITS *)
 
(* BEGIN MODULE TYPE.WORD *)
(* END MODULE TYPE.WORD *)

(* BEGIN MODULE TYPE.LINE *)
(* END MODULE TYPE.LINE *)
  
(* BEGIN MODULE TYPE.FRAG *)
(* END MODULE TYPE.FRAG *)

(* BEGIN MODULE INPLINE *)
(* END MODULE INPLINE *)

(* BEGIN MODULE READLINE *)
(* END MODULE READLINE *)

(* BEGIN MODULE WRITELINE *)
(* END MODULE WRITELINE *)

(* BEGIN MODULE FILEARGS *)
(* END MODULE FILEARGS *)

(* BEGIN MODULE GETFILE *)
(* END MODULE GETFILE *)
 
(* BEGIN MODULE INPWORD *)
(* END MODULE INPWORD *)
 
(* BEGIN MODULE READWORD *)
(* END MODULE READWORD *)

(* BEGIN MODULE WRITEWORD *)
(* END MODULE WRITEWORD *)

(* BEGIN MODULE TOUPPER *)
(* END MODULE TOUPPER *)

(* BEGIN MODULE TOLOWER *)
(* END MODULE TOLOWER *)
 
(* BEGIN MODULE GETCHAR *)
(* END MODULE GETCHAR *)
 
(* BEGIN MODULE GETREAL *)
(* END MODULE GETREAL *)
 
(* BEGIN MODULE GETINTEGER *)
(* END MODULE GETINTEGER *)
 
(* BEGIN MODULE NUC *)
(* END MODULE NUC *)

(* BEGIN MODULE AA *)
(* END MODULE AA *)

(* BEGIN MODULE SAMEWORD *)
(* END MODULE SAMEWORD *)

(* BEGIN MODULE READSEQ *)
(* END MODULE READSEQ *)

(* BEGIN MODULE READPRO *)
(* END MODULE READPRO *)
 
(* BEGIN MODULE LINKED *)
(* END MODULE LINKED *)
 
(* BEGIN MODULE SORT *)
(* END MODULE SORT *)

(* BEGIN MODULE MV *)
(* END MODULE MV *)

(* BEGIN MODULE STARTUP *)
(* END MODULE STARTUP *)
(* END MODULE DUMMY *)
 
(********************* ACTUAL PROCEDURES ***************)
   (**********************************************************)
   (*  LINE I/O Procedures.                                  *)
   (**********************************************************)
(* BEGIN MODULE INPLINE *) 
    (*  Read a line from the console.    *)
    procedure INPLINE(var W:LINE);
      var I : integer;
         CH,BACKSPACE : char;
      begin
        BACKSPACE:=  chr(8);
(*!!!   if eoln(input) then readln; *)
(*!!!   get(input); *)
        with W do begin
          LEN:= 0;
          while not eoln(input) do begin
            read(CH);
(*!!!       if not eoln then*)
              if (CH = BACKSPACE) and (LEN > 0) then begin
                write(BACKSPACE); LEN:= LEN-1 end
              else if LEN < MAXLINE then begin
                LEN := LEN + 1;
                STR[LEN]:= CH
                end
            end;
          readln(input);
          for I := LEN+1 to MAXLINE do STR[I] := ' '
          end
      end; (* INPLINE *)
(* END MODULE INPLINE *)

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
(* END MODULE READLINE *)

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
(* END MODULE WRITELINE         VERSION= 'CSAPMODS     Version  4/13/88'; *)

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
(* END MODULE FILEARGS *)


(* BEGIN MODULE GETFILE *)
(* The c function access is found in sys/file.h. access must be declared in
   the outermost scope of the Pascal program. See access (2) manual pages *)
(*!!!*) function access(PATH:CHARARRAY; MODE:integer):integer; external c;

(************************************************************)
(* Open  files, checking for I/O errors.                    *)
(* Assumes user has already been prompted for filename.     *)
(* Syntax: GETFILE(F,'I',<filename>);                       *)
(*              opens file F for input                      *)
(*         GETFILE(F,'O',<filename>);                       *)
(*              opens file F for output                     *)
(* !!! NOTE: this non-Standard procedure must be rewritten  *)
(* for each different version of Pascal.                    *)
(************************************************************)
procedure GETFILE(var F:text; FTYPE:char; var FILENAME:LINE);

  var CH:char;
      OKAY:boolean;

  (* - - - - - - - - - - - - - - - - - - - - - - - - - - - *)
  (* =true if file exists, =false if it doesn't *)
  (* A comparable procedure is provided with ATT System V Pascal *)
  (* A special thank you to Mark Foster, CIS Rearch Computing, Univ. of Penn.
     and Glenn Gribble of Synaptics Inc. for showing me how to use access()
     within Pascal programs. B.F. *)
  function EXISTS(FILENAME:LINE):boolean;
    const F_OK=0; (* checks if file exists *)
    begin (* EXISTS *)
      with FILENAME do begin
        STR[LEN+1]:=chr(0); (* c strings terminate with null *)
        if access(STR,F_OK)=-1 then EXISTS:=false
        else  EXISTS:=true
        end (* with FILENAME *)
    end; (* EXISTS *)

  begin (* GETFILE - - - - - - - - - - - - - - - - - - - - - *)
      (* Read in filename.*)
      INPLINE(FILENAME);
      case FTYPE of
        (* Input file must exist before being opened.*)
        'I':begin (* input file *)
              while not(EXISTS(FILENAME)) do begin
                writeln('***FILE NOT FOUND');
                writeln('Enter another filename:');
                INPLINE(FILENAME);
                end; (* while *) 
(*!!!*)       reset(F,FILENAME.STR)
            end; (* I *)

        (* If output file already exists, ask user if he really wants
           to over write it. *)
        'O':begin
              repeat
                OKAY:=true;
                if EXISTS(FILENAME) then begin
                  repeat  
                    writeln('*** WARNING! File already exists. Overwrite?[Y|N]');
                    readln(CH)
                  until CH in ['Y','y','N','n'];
                  if not (CH in ['Y','y']) then begin
                    writeln('Enter another filename:');
                    INPLINE(FILENAME);
                    OKAY:=false
                    end
                  end (* EXISTS(FILENAME) *)
              until OKAY;
            rewrite(F,FILENAME.STR)  (* output file *)
            end (* O *)
        end (*case*)
  end;  (* GETFILE *)
(* END MODULE GETFILE         VERSION= 'SUNMODS     Version  3/21/90'; *)

   (**********************************************************)
   (*  WORD I/O Procedures.                                  *)
   (**********************************************************)
(* BEGIN MODULE INPWORD *)
  (* Read a WORD from the terminal. *)
   procedure INPWORD(var W:WORD);
      var I : integer;
         CH : char;
         BACKSPACE: char;
      begin
        with W do begin
          BACKSPACE:= chr(8);
          LEN:= 0;
          while input^ = ' ' do read(input,CH);
          while (input^ <> ' ') and (not eoln(input)) do begin
            read(CH);
(*!!!       if not eoln then*)
              if (CH= BACKSPACE) and (LEN > 0) then begin
                write(BACKSPACE);LEN:= LEN-1 end
              else if LEN < MAXWORD then begin
                LEN := LEN + 1;
                STR[LEN]:= CH
                end
            end;
        (*readln;*)
          for I := LEN+1 to MAXWORD do STR[I] := ' '
          end
      end; (* INPWORD *)
(* END MODULE INPWORD *)
 
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
(* END MODULE READWORD         VERSION= 'CSAPMODS     Version  4/13/88'; *)
 
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
(* END MODULE WRITEWORD         VERSION= 'CSAPMODS     Version  4/13/88'; *)
 
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
(* END MODULE TOUPPER *)

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
(* END MODULE TOLOWER *)

(* BEGIN MODULE GETCHAR *)
    (* Read a character from the console and check *)
    (*  for correct response.                      *)
    procedure GETCHAR(var CH:char;PNAME:LETTERS;ALLOWED:CHSET);
(* BEGIN MODULE TOUPPER *)
(* END MODULE TOUPPER *)
      begin
        writeln;
        repeat
          writeln('Type new value for ',PNAME,'  (CURRENT VALUE: ',CH,')');
          readln(CH); CH:=TOUPPER(CH);
          if CH in ALLOWED then
          else writeln('Inappropriate response: ',CH)
        until CH in ALLOWED
      end; (* GETCHAR *)
(* END MODULE GETCHAR *)
 
(* BEGIN MODULE GETREAL *)
 (* The procedure GETREAL has been adapted from the procedure 'rdr', from
    Jensen,K., and Wirth,N. (1974) Pascal User Manual and Report, 2nd
    Ed., Springer-Verlag, pp122-123.  The scope of real numbers considered
    legal by GETREAL includes all possible real numbers as defined on p111
    with two additional allowances: (i) numbers may begin with a decimal
    point (ii) 'E' or 'e' may be used to indicate exponentiation. *)
 
    procedure GETREAL(var VAL:real;LOW,HIGH:real);
 
    label 776,777,778;
 
(*!!!*)    const  HIGHEXP = 38;
(*!!!*)           LOWEXP = -38;
    type   EXPONENT = 0..HIGHEXP;
 
     var   CH: char;
           LIMIT: real;
           Y : real;
           ORDZERO,I,E : integer;
           LEGAL,INRANGE :boolean;
           NUMWORD  : WORD;
           NEG1,NEG2 : boolean;
           A,S : real;
           LCOUNT,RCOUNT : integer;
           DIGITS: set of '0'..'9';
 
 
      function TEN(E:EXPONENT):real; (* = 10**E *)
        var I:integer;
            T:real;
        begin
          I:=0; T:=1.0;
          repeat
            if odd(E) then
              case I of
                0: T:= T * 1.0E1;
                1: T:= T * 1.0E2;
                2: T:= T * 1.0E4;
                3: T:= T * 1.0E8;
                4: T:= T * 1.0E16;
                5: T:= T * 1.0E32;
        (*!!!   6: T:= T * 1.0E64; Max. exponent is 38
                7: T:= T * 1.0E128; 
                8: T:= T * 1.0E256; *)
                end;
            E:= E div 2; I:= I+1
          until E = 0;
          TEN:= T
        end; (* TEN *)
 
    begin
      LIMIT:= MAXREAL/10;
      DIGITS:=['0'..'9']; ORDZERO:= ord('0');
      repeat
       repeat
        LEGAL:= true; INRANGE:= true;
        INPWORD(NUMWORD);
        with NUMWORD do begin
          (* Evaluate sign, if any *)
          I:=1;
          CH:= STR[I];
          if CH= '-' then begin NEG1:= true; I:= I+1; CH := STR[I] end
          else begin NEG1:= false; if CH='+' then begin
                                         I:= I+1; CH:= STR[I] end end;
          if not (CH in ['0'..'9','.']) then goto 777;
 
          (* Evaluate whole number part (optional) *)
          A:= 0; E:= 0;
          LCOUNT:= 1;
          while CH in DIGITS do begin
            if A < LIMIT then A:= 10*A + (ord(CH)-ORDZERO)
                          else E:= E+1;
            I:= I+1; CH:= STR[I];LCOUNT:=LCOUNT+1
            end;
 
          (* Evaluate fractional part. *)
          RCOUNT:=0;
          if CH = '.' then begin
            I:= I+1; CH:= STR[I];
            while CH in DIGITS do begin
              if A < LIMIT then
                 begin A:= 10*A + (ord(CH)-ORDZERO); E:= E-1 end;
              I:=I+1; CH:= STR[I]; RCOUNT:= RCOUNT+1
              end
            end;
 
          (* Evaluate exponent *)
          if ((CH='E') or (CH='e')) then begin
            I:= I+1; CH:= STR[I]; S:= 0;
            if CH= '-' then begin NEG2:=true;I:= I+1;CH:=STR[I] end
                       else begin NEG2:=false;
                       if CH='+' then begin I:=I+1;CH:=STR[I] end
                       end;
            if CH in DIGITS then begin
              S:= ord(CH)-ORDZERO;
              I:=I+1; CH:= STR[I];
              while CH in DIGITS do begin
                if S < LIMIT then S:= 10*S + (ord(CH)-ORDZERO);
                I:=I+1; CH:= STR[I]
                end
              end
            else goto 777;
            if NEG2 then E:= E-round(S) else E:= E+round(S)
            end;
 
          (* Check for errors  *)
          if I <= LEN then (*illegal char*) goto 777;
          if ((E-RCOUNT<LOWEXP) or (E+LCOUNT>HIGHEXP)) then goto 776;
 
          (* Calculate final value *)
          Y:= A;
          if NEG1 then Y:= -Y;
          if E < 0 then VAL:= Y/(TEN(-E))
           else if E > 0 then VAL:= Y*(TEN(E))
           else VAL:=Y;
          LEGAL:= true;
          goto 778;
 
          (* Error handling statements *)
      776:LEGAL:= false;writeln;writeln(' Exponent Error.');
          writeln(' Enter new value:'); goto 778;
      777:LEGAL:=false;writeln;
          writeln(' Illegal character encountered.');
          write( 'Enter new value:' );
      778:; (* no error, do nothing *)
       end (* with *)
     until LEGAL;
 
 
     (* Make sure number is in range *)
     if ((VAL<LOW) or (VAL>HIGH)) then begin
       INRANGE:= false;
       writeln;
       writeln(' Number is out of range');
       writeln(' Enter new value:')
       end
     else INRANGE:= true
    until INRANGE
   end; (* GETREAL *)
(* END MODULE GETREAL *)
 
(* BEGIN MODULE GETINTEGER *)
  (* Prompts user for an integer, checks whether all characters are digits,*)
  (* and whether number is within desired range; harasses user until valid *)
  (* integer is received. *)
  procedure GETINTEGER(var NUM:integer;LBOUND,HBOUND:integer);
 
   var I,VAL,ORDZERO  : integer;
     LEGAL,INRANGE,NEGATIVE      : boolean;
     CH : char;
     NUMWORD : WORD;
 
  begin
    INRANGE:= false;
    LEGAL := false;
    NEGATIVE := false;
    ORDZERO:= ord('0');
    repeat
      repeat
        INPWORD(NUMWORD);
        with NUMWORD do begin
 
          (* Evaluate sign, if any *)
          I:=1;CH:=STR[I];
          if CH = '-' then begin
             NEGATIVE := true;
             I:=I+1; CH:= STR[I]
             end
          else begin NEGATIVE:= false; if CH='+' then begin
                         I:=I+1; CH:= STR[I] end end;
 
          (* Evaluate unsigned integer *)
          NUM:= 0;
          while CH in ['0'..'9'] do
           begin
             VAL:= ord(CH)-ORDZERO;
             NUM:= (NUM * 10) + VAL;
             I:= I+1; CH:= STR[I]
             end;
           if I > LEN then LEGAL:= true
           else  begin
              LEGAL:= false;
              writeln;
              writeln('Illegal character encountered.');
              writeln('Enter new value:  ');
             end
         end (* with *)
       until LEGAL;
 
       (* If the number entered was negative, multiply *)
       (* NUM by -1. Check range of number.            *)
       if NEGATIVE then NUM:= - NUM;
       if (NUM >= LBOUND) and (NUM <= HBOUND)
         then INRANGE := true
       else begin
         INRANGE := false;
         writeln;
         writeln('Number is out of range.');
         writeln('Please enter new value:')
         end;
     until INRANGE;
   end;   (* GETINTEGER *)
(* END MODULE GETINTEGER *)
 
(* BEGIN MODULE NUC *)
  (*****************************************************************)
  (*  Convert a character to the appropriate nucleotide symbol.    *)
  (*****************************************************************)
  function NUC(CH:char):NUCLEOTIDE;
    begin
      case CH of
       'A','a': NUC:= A;'C','c': NUC:= C;'G','g': NUC:= G; 
       'T','t','U','u': NUC:= T;'R','r': NUC:= R;'M','m': NUC:= M;
       'B','b': NUC:= B;'N','n': NUC:= N;'Y','y': NUC:= Y;
       'K','k': NUC:= K;'D','d': NUC:= D;'S','s': NUC:= S;
       'W','w': NUC:= W;'H','h': NUC:= H;'V','v': NUC:= V
       end
    end;
(* END MODULE NUC         VERSION= 'CSAPMODS     Version  4/13/88'; *)
 
(* BEGIN MODULE AA *)
  (*****************************************************************)
  (*  Convert a character to the appropriate amino acid symbol.    *)
  (*****************************************************************)
  function AA(CH:char):AMINOACID;
    begin
      case CH of
        'A':AA:=ALA; 'C':AA:=CYS; 'D':AA:=ASP;    'E':AA:=GLU; 'F':AA:=PHE;
        'G':AA:=GLY; 'H':AA:=HIS; 'I':AA:=ILE;    'K':AA:=LYS; 'L':AA:=LEU;
        'M':AA:=MET; 'N':AA:=ASN; 'P':AA:=PRO;    'Q':AA:=GLN; 'R':AA:=ARG;
        'S':AA:=SER; 'T':AA:=THR; 'V':AA:=VALINE; 'W':AA:=TRP; 'Y':AA:=TYR;
        'B':AA:=ASX; 'Z':AA:=GLX; 'X':AA:=UNKX;   '*':AA:=TERM
       end
    end; (* AA *)
(* END MODULE AA *)

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
(* END MODULE SAMEWORD  9/12/91 *)

(* BEGIN MODULE READSEQ *)
    (*  ******************************************* *)
    (*  Read a DNA or RNA sequence from SFILE       *)
    (*   and store it in S.                         *)
    (*  ******************************************* *)
    procedure READSEQ(var SFILE:text; var S:SEQUENCE; var SEQLEN:integer;
                      var NAME:WORD; var CIRCULAR:boolean);
      var CH,FILETYPE,LASTNUM: char;
          TOOLONG,BADFILE:boolean;
          ID,LOCUS,ORI,CIRC:WORD;
          I:integer;
      begin
        (* Prompt user for sequence file type *)    
        writeln('The following file formats can be read:');
        writeln('  F:free format   B:BIONET    G:GENBANK');
        repeat
          writeln('Type letter of format (F|B|G)');
          readln(CH)
        until CH in ['F','f','B','b','G','g'];
        case CH of
          'F','f':FILETYPE:='F';
          'B','b':FILETYPE:='B';
          'G','g':FILETYPE:='G'
          end; (* case *)
        writeln('Reading input file...');
        BADFILE:=false;
        NAME.LEN:=0;
        
        (* For BIONET or GENBANK, read in sequence name and topology *)
        (* Advance to beginning of sequence *)
        if FILETYPE= 'B' then begin
          (* First non-comment line is name. Name may be blank.*)
          while SFILE^=';' do readln(SFILE);
          while (SFILE^=' ') and (not eoln(SFILE)) do get(SFILE);
          if not eoln(SFILE) then READWORD(SFILE,NAME)
          else NAME.LEN:=0;
          if not eof(SFILE) then readln(SFILE)
          end (* BIONET *)
          
        else if FILETYPE='G' then begin
          (* Initialize identifiers *)
          with LOCUS do begin
            STR[1]:='L';STR[2]:='O';STR[3]:='C';STR[4]:='U';STR[5]:='S';
            LEN:=5
            end;
          with ORI do begin
            STR[1]:='O';STR[2]:='R';STR[3]:='I';STR[4]:='G';STR[5]:='I';
            STR[6]:='N';LEN:=6
            end;
          with CIRC do begin
            STR[1]:='c';STR[2]:='i';STR[3]:='r';STR[4]:='c';STR[5]:='u';
            STR[6]:='l';STR[7]:='a';STR[8]:='r';LEN:=8
            end;
           (* Advance to LOCUS line. Read in NAME and topology. *)
           while not((SFILE^='L') or (eof(SFILE))) do readln(SFILE);
           if not eof(SFILE) then begin
             READWORD(SFILE,ID);
             if SAMEWORD(ID,LOCUS) then begin
                if not eof(SFILE) then READWORD(SFILE,NAME);
                if not eof(SFILE) then READWORD(SFILE,ID); (* skip seq. length *)
                if not eof(SFILE) then READWORD(SFILE,ID); (* skip. 'bp' *)
                (* After 'bp', there is an optional field telling
                   the type of molecule (ss-RNA, ds-DNA etc.)
                   Since this field is optional, we must test the
                   next two tokens to see if they are 'circular' *)
                CIRCULAR:=false;
                if not eof(SFILE) then READWORD(SFILE,ID);
                if SAMEWORD(ID,CIRC) then CIRCULAR:=true
                else begin
                     if not eof(SFILE) then READWORD(SFILE,ID);
                     if SAMEWORD(ID,CIRC) then CIRCULAR:=true
                     end
                end (* SAMEWORD(ID,LOCUS) *)
              else BADFILE:=true
             end;
          
           (* Advance to ORIGIN line. Sequence begins on next line *)
           if not eof(SFILE) then begin
             repeat
               readln(SFILE);
               if SFILE^='O' then READWORD(SFILE,ID)
             until SAMEWORD(ID,ORI) or eof(SFILE);
             if SAMEWORD(ID,ORI) then readln(SFILE)
             end;
           if eof(SFILE) then BADFILE:= true
           end; (* GENBANK *)
        
        (* Read in sequence *)
        SEQLEN := 0;TOOLONG:=false;
        if not BADFILE then begin
          while not eof(SFILE) do begin
            while not eoln(SFILE) do begin
              read(SFILE,CH);
              if CH in ['A','a','C','c','G','g','T','t','U','u','N','n',
                        'R','r','Y','y','M','m','W','w','S','s','K','k',
                        'D','d','H','h','V','v','B','b'] then
                if SEQLEN < MAXSEQ-2 then begin
                (*write(CH);*)
                  SEQLEN := SEQLEN + 1;
                  S[SEQLEN]:= NUC(CH);
                  end
                else TOOLONG:= true
              else if CH=';' then (*begin    comment in input file *)
                                    readln(SFILE)(*;writeln end*)
              else if CH in ['1','2'] then LASTNUM:=CH                      
              end;
            readln(SFILE);(*writeln*)
            end;
          if TOOLONG then writeln(
            '>>> WARNING! Sequence length exceeds MAXSEQ-2. Seq. truncated.');

          if FILETYPE='F' then begin
            repeat
              writeln('Is sequence circular or linear? (Type C or L)');
              readln(CH)
            until CH in ['C','c','L','l'];
            case CH of
              'C','c':CIRCULAR:= true;
              'L','l':CIRCULAR:= false
              end
            end  
          else if FILETYPE='B' then 
            if LASTNUM='1' then CIRCULAR:=false
            else CIRCULAR:=true
          end
         else writeln('>>> ERROR: Not a GENBANK file. No sequence read.');
(*!!!   CLOSE(SFILE);*)
    end; (* READSEQ *)
(* END MODULE READSEQ *)

(* BEGIN MODULE READPRO *)
    (*  ******************************************* *)
    (*  Read a protein    sequence from SFILE       *)
    (*   and store it in PR.                        *)
    (*  ******************************************* *)
    procedure READPRO(var SFILE:text; var PR:PROTEIN; var LEN:integer;
                      var NAME:WORD);
      label 86;

      var CH,ANSWER,FILETYPE: char;
          J: integer;
          TOOLONG:boolean;
 
      begin
        (* Prompt user for sequence file type *)
        writeln('The following file formats can be read:');
        writeln('  F:free format    N:NBRF');
        repeat
          writeln('Type letter of format (F|N)');
          readln(ANSWER)
        until ANSWER in ['F','f','N','n'];
        case ANSWER of
          'F','f':FILETYPE:='F';
          'N','n':FILETYPE:='N'
          end;
        writeln('Reading input file...');

        (* NBRF: read sequence name and title line *)
        (* 'Official' NBRF files have a name line followed by a title line.
            In the name line, four data characters preceed the name itself.
            These are deleted by the program.  Short NBRF files have a title
            following the name, on the same line. The name is not preceeded 
            by data characters. *)
        NAME.LEN:=0;
        if FILETYPE='N' then begin
          if SFILE^='>' then read(SFILE,CH);
          READWORD(SFILE,NAME);
          if NAME.STR[3]=';' then with NAME do begin
             for J:= 4 to LEN do STR[J-3]:=STR[J];
             LEN:=LEN-3;
             readln(SFILE)
             end; (* with NAME *)
          readln(SFILE)
          end; (* NBRF *)

        (* Read in the sequence *)
        J := 0;TOOLONG:=false;
        while not eof(SFILE) do begin
          while not eoln(SFILE) do begin
            read(SFILE,CH);
            if CH in ['G','A','V','L','I','M','F','P','S','T','C','N','Q',
                      'Y','W','D','E','H','K','R','B','Z','*','X'] then
              if J < MAXSEQ-MAXRANGE then begin
                J := J + 1;
                case CH of
                  'G':PR[J]:= GLY; 'A':PR[J]:= ALA; 'V':PR[J]:= VALINE;
                  'L':PR[J]:= LEU; 'I':PR[J]:= ILE; 'M':PR[J]:= MET;
                  'F':PR[J]:= PHE; 'P':PR[J]:= PRO; 'S':PR[J]:= SER;
                  'T':PR[J]:= THR; 'C':PR[J]:= CYS; 'N':PR[J]:= ASN;
                  'Q':PR[J]:= GLN; 'Y':PR[J]:= TYR; 'W':PR[J]:= TRP;
                  'D':PR[J]:= ASP; 'E':PR[J]:= GLU; 'H':PR[J]:= HIS;
                  'K':PR[J]:= LYS; 'R':PR[J]:= ARG; 'X':PR[J]:= UNKX; 
                  'B':PR[J]:= ASX; 'Z':PR[J]:= GLX; 
                  '*':begin  
                        PR[J]:= TERM;
                        if FILETYPE='N' then goto 86 (*ignore rest of file*)
                      end
                  end
                end
              else TOOLONG:= true
            else if CH=';' then readln(SFILE) (* comment in input file *)
            end;
          readln(SFILE)
          end;

(* !!! *) 86:(*CLOSE(SFILE); *) (* branch destination for end of NBRF seq. *)
        if TOOLONG then writeln(
          '>>> WARNING! Sequence exceeds MAXSEQ-MAXRANGE. Sequence truncated.');
        LEN := J
      end; (* READPRO *)
(* END MODULE READPRO *) 

(* BEGIN MODULE LINKED *)
  (*********************************************************)
  (*  Linked-list operations for restriction fragment list.*)
  (*********************************************************)
 
  (*Get a new fragment from freelist.*)
  procedure GETFRAG(var NEWFRAG:FRAG);
    begin
      if FREEFRAG = nil then new(NEWFRAG)
      else begin
        NEWFRAG:= FREEFRAG;
        FREEFRAG:= FREEFRAG^.NEXT
      end
    end;
 
  (*Add a fragment after DEST*)
  procedure ADDFRAG(var AFRAG,DEST:FRAG);
    var TEMP: FRAG;
    begin
      TEMP:= DEST^.NEXT;
      DEST^.NEXT:= AFRAG; AFRAG^.PREV:= DEST;
      AFRAG^.NEXT:= TEMP;TEMP^.PREV:= AFRAG
    end;
 
  (*Return a list to the top of freelist*)
  procedure RIDOF(var HEAD,TAIL:FRAG);
    var TEMPHEAD,TEMPTAIL:FRAG;
    begin
      if HEAD^.NEXT <> TAIL then begin
        TEMPHEAD:= HEAD^.NEXT;TEMPTAIL:= TAIL^.PREV;
        HEAD^.NEXT:= TAIL;TAIL^.PREV:= HEAD;
        TEMPHEAD^.PREV:= nil;TEMPTAIL^.NEXT:= FREEFRAG;
        FREEFRAG:= TEMPHEAD
        end
    end;
(* END MODULE LINKED *)
 
(* BEGIN MODULE SORT *)
  (*  Invariant:  The array elements > TOP are sorted.*)
  (*    TOP >= unsorted elements.                     *)
  procedure BUBBLESORT(TOP,BOTTOM:integer);
  var SMALLEST,NEXT :integer;
  procedure SWAP(var FIRST,SECOND:FRAG);
    var TEMP:FRAG;
    begin
      TEMP:= FIRST;FIRST:=SECOND;SECOND:=TEMP
    end;
  begin
    while TOP >= BOTTOM do begin
      (*bubble smallest unsorted number to the top of sorted list*)
      SMALLEST:= BOTTOM; NEXT:= BOTTOM+1;
      while NEXT <= TOP do begin
        if ORDER[SMALLEST]^.SIZE < ORDER[NEXT]^.SIZE then
          SWAP(ORDER[SMALLEST],ORDER[NEXT]);
        SMALLEST:= NEXT;
        NEXT:= SMALLEST+1
        end;
      TOP:= TOP - 1
      end;
  end; (* BUBBLESORT *)
(* END MODULE SORT *)

(* BEGIN MODULE MV *)
   (*****************************************************************)
   (*  Calculate the mean and variance of first N elements          *)
   (*  in array X.  The type VECTOR must be declared in a higher    *)
   (*  scope as   array[1..<positive constant>] of <real|integer>   *)
   (*****************************************************************)
   procedure MV(var X:VECTOR; N:integer; var MEAN,VARIANCE:real);
     var SX,SSX:real;
         I:integer;
     begin
       SX:=0; SSX:=0;
       for I:= 1 to N do begin
         SX:= SX + X[I];
         SSX:= SSX + sqr(X[I])
         end;
       MEAN:= SX/N;
       VARIANCE:= (SSX - (sqr(SX)/N))/N-1
     end; (* MV *)
(* END MODULE MV *)
  
  (* ----------------------------------------------------------  *)
  (* ----------------- MAIN  PROCEDURE  -----------------------  *)
    begin
(* BEGIN MODULE STARTUP *)
(* Peform operations which must be done at beginning of main
   procedure. *)
(*!!!   TERMIN(input);    Open input for interactive use *)
(*!!!   TERMOUT(output);   "   output "      "        "  *)
     writeln(VERSION:50);
     writeln;
(* END MODULE STARTUP *)
 
      writeln('This program tests the module procedures in CSAPMOD.');
      writeln;
      writeln('TEST OF INTERACTIVE INPUT');
      CH:= ' ';
      writeln('When prompted, type a letter between A and Z');
      GETCHAR(CH,'CHARACTER:',['A'..'Z']);
      writeln('The character is: ',CH);
      writeln('Type a single word and press RETURN');
      INPWORD(WO); readln; 
      write('The word you typed was: ');
      WRITEWORD(output,WO,WO.LEN);writeln;
      writeln('Type a line of text and press RETURN');
      INPLINE(LI);
      writeln('The line you typed was:');
      WRITELINE(output,LI,LI.LEN);
      writeln;
      writeln('TEST OF NUMERICAL INPUT');
      writeln('Type an integer between ',MININT,' and ',MAXINT,
              'and press RETURN');
      GETINTEGER(INT,MININT,MAXINT);readln;
      writeln('The integer you typed was: ',INT);
      writeln('Type a real number between ',-MAXREAL,' and ',MAXREAL);
      writeln('and press RETURN');
      GETREAL(RE,-MAXREAL,MAXREAL);readln;
      writeln('The real number you typed was: ',RE);
      writeln;
      writeln('TEST OF FILE I/O');
      writeln('Enter sequence filename:');
      GETFILE(INFILE,'I',LI);
      READSEQ(INFILE,SEQ,SEQLEN,WO,CIRCULAR);
      writeln('The sequence in');
      WRITELINE(output,LI,LI.LEN);writeln;
      writeln('is ',SEQLEN,' bases long.')
    end. (* CSAPMOD *)
 
