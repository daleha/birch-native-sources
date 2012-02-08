  (* ********************************************************  *)
  (*                                                           *)
  (*   D3HOM     VERSION  8/13/2001, Standard Pascal           *)
  (*             Brian Fristensky                              *)
  (*             Dept. of Plant Science                        *)
  (*             University of Manitoba                        *)
  (*             Winnipeg, MB R3T 2N2   CANADA                 *)
  (*                                                           *)
  (* Copyright (c) 1984,1986,1987, 1990 by Brian Fristensky.   *)
  (* !!! in comment indicates feature which may need change.   *)
  (*  *******************************************************  *)
 
(*!!!*)program D3HOM(input, output (*, SFILEX,SFILEY,OUTFILE*));
 (*!!! Some Pascals require file parameters in heading *)
 
  const MAXSEQ  =  750000;
(* BEGIN MODULE REALLIMITS *)
(*!!!*)MAXREAL = 1.7E38;
(* END MODULE REALLIMITS         VERSION= 'SUNMODS     Version  6/26/01'; *)
        MAXWORD = 25;
        MAXRANGE= 30;
        MAXLINE = 150;
        VERSION = 'D3HOM         Version  8/13/2001';
 
  type NUCLEOTIDE  = (T,C,A,G,N,R,Y,M,W,S,K,D,H,V,B,Z);
       SEQUENCE    = array[-MAXRANGE..MAXSEQ] of NUCLEOTIDE;
       NP          = ^NODE;
       NODE        = record POS:integer;
                            NEXT:NP
                            end;
(* BEGIN MODULE TYPE.WORD *)
       (*   <word>::= <non-blank char>[<non-blank char>] *)
       WORD    = record
                 LEN:integer;
                 STR:array[1..MAXWORD] of char
                 end;
(* END MODULE TYPE.WORD         VERSION= 'SUNMODS     Version  6/26/01'; *)

(* BEGIN MODULE TYPE.LINE *)
       CHARARRAY = packed array[1..MAXLINE] of char;
       LINE = record
                STR:CHARARRAY;
                LEN:integer
                end;
(* END MODULE TYPE.LINE         VERSION= 'SUNMODS     Version  6/26/01'; *)
 
  var SFILEX,SFILEY,OUTFILE: text;
      XFN,YFN,OFN,HLINE:LINE;
      SEQX,SEQY    : SEQUENCE;
      LENX,LENY    : integer;
      NAMEX,NAMEY  : WORD;
      TRIPLET      : array[T..G,T..G,T..G] of NP;
      FREELIST     : NP;
      NUCHAR       : array[NUCLEOTIDE] of char;
      STARTX,FINISHX,
      STARTY,FINISHY,
      HOMRANGE,MINPER,COMPFACT,LINESIZE,
      MAXSCORE,THREEMATCH,MINSCORE: integer;
      SCALEFAC,MAXFACTOR,TEMP:real;
      SUBSCORE     : array[0..MAXRANGE] of integer;
      KIND:char;
      PC: array[-1..100] of char; (*Printing characters for graph*)
      INVERSEX,CIRCX,CIRCY : boolean;
      CHOICE,J:integer;
      
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
(* END MODULE INPLINE         VERSION= 'SUNMODS     Version  6/26/01'; *)

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
(* END MODULE WRITELINE         VERSION= 'SUNMODS     Version  6/26/01'; *)
 
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
(* END MODULE GETFILE         VERSION= 'SUNMODS     Version  6/26/01'; *)
 
  (**************************************************)
  (*  WORD   I/O  PROCEDURES                        *)
  (**************************************************)
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
(* END MODULE INPWORD         VERSION= 'SUNMODS     Version  6/26/01'; *)

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
(* END MODULE READWORD         VERSION= 'SUNMODS     Version  6/26/01'; *)
 
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
(* END MODULE WRITEWORD         VERSION= 'SUNMODS     Version  6/26/01'; *)
 
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
(* END MODULE GETREAL         VERSION= 'SUNMODS     Version  6/26/01'; *)

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
(* END MODULE NUC         VERSION= 'SUNMODS     Version  6/26/01'; *)

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
(* END MODULE SAMEWORD         VERSION= 'SUNMODS     Version  6/26/01'; *)

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
(* END MODULE READSEQ         VERSION= 'SUNMODS     Version  6/26/01'; *)
 
  (*****************************************************)
  (* Initialization procedures.                        *)
  (*****************************************************)
  procedure INITIALIZE;
    var N1,N2,N3:NUCLEOTIDE;
    begin
      (*Set default values for parameters *)
      HOMRANGE:= 10; SCALEFAC:=0.95;
      MINPER:=60; COMPFACT:=10; KIND:='D'; LINESIZE:= 70;
 
      (*   Initialize NUCHAR array which holds the            *)
      (*   character values of the NUCLEOTIDES.               *)
      NUCHAR[A]:='A'; NUCHAR[C]:='C';NUCHAR[G]:='G'; NUCHAR[T]:='T';
      NUCHAR[R]:='R'; NUCHAR[Y]:='Y';NUCHAR[M]:='M'; NUCHAR[W]:='W';
      NUCHAR[S]:='S'; NUCHAR[K]:='K';NUCHAR[D]:='D'; NUCHAR[H]:='H';
      NUCHAR[V]:='V'; NUCHAR[B]:='B';NUCHAR[N]:='N'; NUCHAR[Z]:='N';
 
      (* Initialize the PC array which holds the       *)
      (* characers which symbolize percent identity.   *)
      PC[ -1]:=' ';PC[  0]:='.';
      PC[ 25]:='m';PC[ 26]:='l';PC[ 27]:='l';PC[ 28]:='k';PC[ 29]:='k';
      PC[ 30]:='j';PC[ 31]:='j';PC[ 32]:='i';PC[ 33]:='i';PC[ 34]:='h';
      PC[ 35]:='h';PC[ 36]:='g';PC[ 37]:='g';PC[ 38]:='f';PC[ 39]:='f';
      PC[ 40]:='e';PC[ 41]:='e';PC[ 42]:='d';PC[ 43]:='d';PC[ 44]:='c';
      PC[ 45]:='c';PC[ 46]:='b';PC[ 47]:='b';PC[ 48]:='a';PC[ 49]:='a';
      PC[ 50]:='Z';PC[ 51]:='Z';PC[ 52]:='Y';PC[ 53]:='Y';PC[ 54]:='X';
      PC[ 55]:='X';PC[ 56]:='W';PC[ 57]:='W';PC[ 58]:='V';PC[ 59]:='V';
      PC[ 60]:='U';PC[ 61]:='U';PC[ 62]:='T';PC[ 63]:='T';PC[ 64]:='S';
      PC[ 65]:='S';PC[ 66]:='R';PC[ 67]:='R';PC[ 68]:='Q';PC[ 69]:='Q';
      PC[ 70]:='P';PC[ 71]:='P';PC[ 72]:='O';PC[ 73]:='O';PC[ 74]:='N';
      PC[ 75]:='N';PC[ 76]:='M';PC[ 77]:='M';PC[ 78]:='L';PC[ 79]:='L';
      PC[ 80]:='K';PC[ 81]:='K';PC[ 82]:='J';PC[ 83]:='J';PC[ 84]:='I';
      PC[ 85]:='I';PC[ 86]:='H';PC[ 87]:='H';PC[ 88]:='G';PC[ 89]:='G';
      PC[ 90]:='F';PC[ 91]:='F';PC[ 92]:='E';PC[ 93]:='E';PC[ 94]:='D';
      PC[ 95]:='D';PC[ 96]:='C';PC[ 97]:='C';PC[ 98]:='B';PC[ 99]:='B';
      PC[100]:='A';
 
      (* Initialize TRIPLET table to nil *)
      for N1:= T to G do
        for N2:= T to G do
          for N3:= T to G do TRIPLET[N1,N2,N3]:= nil
    end; (* INITIALIZE *)
 
    (* **************************************************** *)
    (* Open a sequence file and read in sequence.           *)
    (* **************************************************** *)
    procedure NEWSEQ(var SFILE:text; var FN:LINE; var S:SEQUENCE; 
                     var LEN,START,FINISH:integer;var NAME:WORD;
                     var CIRC:boolean; AXIS:char);

      var ANSWER:char;
          J:integer;
      (* Initialize negative and last parts of the sequence*)
      procedure SEQENDS(var S:SEQUENCE; var LEN:integer; var CIRCULAR:boolean;
                       UNKNOWN:NUCLEOTIDE);
        var I,J: integer;
        begin
          if CIRCULAR then begin
             J:= LEN -(MAXRANGE+1);
             for I:= -MAXRANGE to 0 do begin S[I]:=S[J]; J:= J+1 end;
             J:=1;
             for I:= LEN+1 to LEN+MAXRANGE do begin S[I]:= S[J];J:=J+1 end
             end
          else begin
            for I:= -MAXRANGE to 0 do S[I]:= UNKNOWN;
            for I:= LEN+1 to LEN+MAXRANGE do S[I]:= UNKNOWN
            end
         end; (* SEQENDS *)

      begin
        writeln('Enter ',AXIS,'-axis sequence filename:');
        GETFILE(SFILE,'I',FN);
        NAME.LEN:=0;
        READSEQ(SFILE,S,LEN,NAME,CIRC);
        if NAME.LEN=0 then begin
          writeln('Type name for ',AXIS,'-axis sequence to appear on output:');
          INPWORD(NAME);readln
          end;
        START:=1; FINISH:=LEN; 
        case AXIS of
          'X':begin
                SEQENDS(S,LEN,CIRC,Z);            
                for J:= 1 to LEN do if S[J]=N then S[J]:=Z;
                repeat
                  writeln('Is seq. on X-axis the inverse strand? [Y/N]');
                  readln(ANSWER)
                until ANSWER in ['Y','y','N','n'];
                if ANSWER in ['Y','y'] then begin
                   INVERSEX:= true; START:= LEN; FINISH:=1 end
                else INVERSEX:=false   
              end;
          'Y':SEQENDS(S,LEN,CIRC,N);
          end (* case *)
      end; (* NEWSEQ *)

    (* **************************************************** *)
    (* Prompt user for parameters used by program.          *)
    (* **************************************************** *)
    procedure PARAMETERS;
      type LETTERS = packed array[1..10] of char;
           CHSET   = set of 'A'..'Z';
 
      var RESPONSE :integer;
          TEMP     :real;
 
      (*  Read an integer parameter from the console and check *)
      (*    that it is in range.                               *)
      procedure GETNUMBER(var P:integer;PNAME:LETTERS;LOW,HIGH:integer);
        var TEMP:real;
        begin
          writeln;
          writeln('Type new value for ',PNAME,'  (CURRENT VALUE: ',P,')');
          GETREAL(TEMP,LOW,HIGH);readln;
          P:= round(TEMP)
        end; (* GETNUMBER *)
 
(* BEGIN MODULE GETCHAR *)
    (* Read a character from the console and check *)
    (*  for correct response.                      *)
    procedure GETCHAR(var CH:char;PNAME:LETTERS;ALLOWED:CHSET);
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
(* END MODULE TOUPPER         VERSION= 'SUNMODS     Version  6/26/01'; *)
      begin
        writeln;
        repeat
          writeln('Type new value for ',PNAME,'  (CURRENT VALUE: ',CH,')');
          readln(CH); CH:=TOUPPER(CH);
          if CH in ALLOWED then
          else writeln('Inappropriate response: ',CH)
        until CH in ALLOWED
      end; (* GETCHAR *)
(* END MODULE GETCHAR         VERSION= 'SUNMODS     Version  6/26/01'; *)
 
    (* Display  parameters on screen *)
    procedure DISPLAY;
      procedure HEADLINE(NAME:WORD; CIRC:boolean; INVERSE:boolean; LEN:integer;
                       AXIS:char);
        begin
          write(AXIS,'-axis: ');WRITEWORD(output,NAME,20);
          case CIRC of
            true:write('Topology: CIRCULAR');
            false:write('Topology:   LINEAR')
            end;
          case INVERSE of 
            true:write(' (INVERSE)');
            false:write('          ')
            end;
          writeln('Length: ':10,LEN:11,' nt');
        end; (* HEADLINE *)
     begin
     page(output);
     HEADLINE(NAMEX,CIRCX,INVERSEX,LENX,'X');
     HEADLINE(NAMEY,CIRCY,false,LENY,'Y');
     WRITELINE(output,HLINE,80);writeln;
     writeln('Parameter   Description/Response                 Value':66);
     WRITELINE(output,HLINE,80);writeln;
     writeln(' 1)STARTX   first nucleotide position in SEQX   ':60,STARTX:6);
     writeln(' 2)FINISHX  last  nucleotide position in SEQX   ':60,FINISHX:6);
     writeln(' 3)STARTY   first nucleotide position in SEQY   ':60,STARTY:6);
     writeln(' 4)FINISHY  last  nucleotide position in SEQY   ':60,FINISHY:6);
     writeln(' 5)HOMRANGE dist.from central triplet in a match':60,HOMRANGE:6);
     writeln(' 6)SCALEFAC scale factor for exponential curve  ':60,
                 SCALEFAC:6:2);
     writeln(' 7)MINPER   minimum percent similarity printed  ':60,MINPER:6);
     writeln(' 8)COMPFACT graph compression factor            ':60,COMPFACT:6);
     writeln(' 9)KIND     D:DNA            R:RNA              ':60,KIND:6);
     writeln('10)LINESIZE width of output line (eg.70,120)    ':60,LINESIZE:6);
     WRITELINE(output,HLINE,80);writeln;
     writeln('Type number of parameter you wish to change',
             ' (0 to continue)')
     end; (* DISPLAY *)
 
      begin
        (* Prompt user for new parameter values *)
        repeat
          DISPLAY;
          GETREAL(TEMP,0,10);readln;
          RESPONSE:= round(TEMP);
          if RESPONSE in [1..10] then
            case RESPONSE of
              1:GETNUMBER(STARTX,'STARTX    ',1,LENX);
              2:GETNUMBER(FINISHX,'FINISHX   ',1,LENX);
              3:GETNUMBER(STARTY,'STARTY    ',1,LENY);
              4:GETNUMBER(FINISHY,'FINISHY   ',1,LENY);
              5:GETNUMBER(HOMRANGE,'HOMRANGE  ',1,MAXRANGE);
              6:begin
                  writeln('Type new value for SCALEFAC:');
                  GETREAL(SCALEFAC,0,1);readln
                end;
              7:GETNUMBER(MINPER,'MINPER    ',40,100);
              8:GETNUMBER(COMPFACT,'COMPFACT  ',1,500);
              9:GETCHAR(KIND,'KIND      ',['D','R']);
             10:GETNUMBER(LINESIZE,'LINESIZE  ',40,MAXLINE)
             end
         until RESPONSE= 0;
         case KIND of
           'D':NUCHAR[T]:='T';
           'R':NUCHAR[T]:='U'
           end
     end; (* PARAMETERS *)

  (****************************************************************)
  (* Calculate a table of scores for matches at each distance from*)
  (* the central triplet in a local homology.                     *)
  (****************************************************************)
  procedure CALCSCORES;
    var I,V:integer;
        S:real;
    begin
      V:=100;
      SUBSCORE[0]:= V; MAXSCORE:= V; S:=SCALEFAC;
      for I:= 1 to HOMRANGE do begin
        SUBSCORE[I]:= round(V*S);
        MAXSCORE:= MAXSCORE + (2*SUBSCORE[I]);
        S:= S * SCALEFAC
        end;
      THREEMATCH:= SUBSCORE[0] + (2*SUBSCORE[1]);
      MINSCORE:= round((MINPER/100)*MAXSCORE);
      MAXFACTOR:= 100/MAXSCORE
    end; (*CALCSCORES*)
 
  (* ******************************************************************* *)
  (*  Find similarities (2*HOMRANGE)+1 long with MINPER or better match. *)
  (* ******************************************************************* *)
  procedure QUICKSEARCH(var OUTFILE:text);
 
    var POSX,POSY,RIGHT,RIGHTLIM,LEFT,I,J,K,INDEX,NEXTLINE:integer;
        ONELINE,TENLINE,THISLINE: array[1..MAXLINE] of -1..100;
        N1,N2,N3:NUCLEOTIDE;

    (* Print the header listing graph parameters. *)
    procedure HEADER;
      begin
        writeln(OUTFILE,VERSION);
        write(OUTFILE,'X-axis: ');
        WRITEWORD(OUTFILE,NAMEX,NAMEX.LEN);
        writeln(OUTFILE);
        write(OUTFILE,'Y-axis: ');
        WRITEWORD(OUTFILE,NAMEY,NAMEY.LEN);
        writeln(OUTFILE);
        writeln(OUTFILE,
          'SIMILARITY RANGE:',HOMRANGE:4,'MIN.PERCENT SIMILARITY:':29,MINPER:4);
        writeln(OUTFILE,'SCALE FACTOR:',SCALEFAC:8:2,
               'COMPRESSION:':18,COMPFACT:15)
      end; (*HEADER*)

    (* Print a horizontal axis *)
    procedure HORAXIS(LEFT,RIGHT:integer);
      var NUMBER,DELTA,I,NUMPRINTED:integer;
      begin
        writeln(OUTFILE);
        (* Write numbers *)
        write(OUTFILE,' ':7);
        DELTA:= 10*COMPFACT;
        NUMPRINTED:= (RIGHT-LEFT+1) div DELTA;
        if INVERSEX then begin
          DELTA:= -DELTA;
          NUMBER:= LENX - LEFT + 2
          end
        else NUMBER:= LEFT-1;
        for I:= 1 to NUMPRINTED do begin
          NUMBER:= NUMBER + DELTA;
          write(OUTFILE,NUMBER:10)
          end;
        writeln(OUTFILE);
 
        (* Write sequence if COMPFACT = 1*)
        if COMPFACT = 1 then begin
          write(OUTFILE,' ':7);
          for I:= LEFT to RIGHT do write(OUTFILE,NUCHAR[SEQX[I]])
          end;
        writeln(OUTFILE)
      end; (*HORAXIS*)
 
   (*  Make a table of locations of trinucleotides in SEQX.  Each triplet*)
   (*  has a stack of nodes, each of which holds a location in SEQX      *)
   (*  at which the trinucleotide occurs.                                *)
   procedure MAKETABLE(LEFT,RIGHT:integer);
     var I:integer;
         N1,N2,N3: NUCLEOTIDE;
 
     procedure ADDNODE(var N:NP);
       var TEMP: NP;
       begin
         TEMP:= N;
         if FREELIST = nil then new(N)
         else begin N:= FREELIST; FREELIST:= FREELIST^.NEXT end;
         N^.NEXT:= TEMP;
         N^.POS:= I
       end;
 
     begin
       for I:= RIGHT downto LEFT do begin
         N1:= SEQX[I-1]; N2:= SEQX[I]; N3:= SEQX[I+1];
         if ((N1<N) and (N2<N) and (N3<N)) then ADDNODE(TRIPLET[N1,N2,N3])
         end
    end; (*MAKETABLE*)
 
   (* At each occurrence of the triplet in SEQX, compare a region  *)
   (* HOMRANGE bases on either side.  If the match is good enough, *)
   (* print a character at the corresponding point in the matrix.  *)
   procedure SEARCHFOR(N:NP);
 
     var LX,LY,RX,RY,DISTANCE,SCORE,PERCENT,X:integer;
 
     begin
       while N <> nil do begin
         POSX:= N^.POS;
         SCORE:= THREEMATCH;
         LX:= POSX-2; RX:= POSX+2; LY:=POSY-2; RY:=POSY+2;
         DISTANCE:= 2;
         while DISTANCE <= HOMRANGE do begin
           if SEQX[LX] = SEQY[LY] then SCORE:= SCORE + SUBSCORE[DISTANCE];
           if SEQX[RX] = SEQY[RY] then SCORE:= SCORE + SUBSCORE[DISTANCE];
           LX:= LX-1; RX:= RX+1; LY:= LY-1; RY:= RY+1;
           DISTANCE:= DISTANCE+1
           end;
         if SCORE >= MINSCORE then begin
           PERCENT:= round(SCORE*MAXFACTOR);
           X:= ((POSX-LEFT) div COMPFACT) +1;
           if PERCENT > THISLINE[X] then THISLINE[X]:= PERCENT
           end;
         N:= N^.NEXT
       end
     end; (* SEARCHFOR *)
 
  (* Push the linked-list (if there is one) for each TRIPLET to *)
  (* the top of the FREELIST and reset TRIPLET to nil.          *)
  procedure RIDOF;
    var N1,N2,N3:NUCLEOTIDE;
        HEAD,TAIL:NP;
    begin
      for N1:= T to G do
        for N2:= T to G do
          for N3:= T to G do
            if TRIPLET[N1,N2,N3] <> nil then begin
              HEAD:= TRIPLET[N1,N2,N3];
              TAIL:= HEAD;
              while TAIL^.NEXT <> nil do TAIL:= TAIL^.NEXT;
              TAIL^.NEXT:= FREELIST;
              FREELIST:= HEAD;
              TRIPLET[N1,N2,N3]:= nil
              end
    end; (*RIDOF*)
 
    begin
      (* Initialize LEFT & RIGHT, which define the part of SEQX *)
      (* to be compared with SEQY.                              *)
      if INVERSEX then begin
        LEFT:= LENX-STARTX+1;
        RIGHTLIM:= LENX-FINISHX+1
        end
      else begin
        LEFT:= STARTX;
        RIGHTLIM:= FINISHX
        end;
      RIGHT:= LEFT + (LINESIZE*COMPFACT) -1;
      if RIGHT > RIGHTLIM then RIGHT:= RIGHTLIM;
 
      (* Initialize line templates *)
      J:= 1;
      for I:= 1 to LINESIZE do begin
        TENLINE[I]:= 0; (* period *)
        if J = 10 then begin ONELINE[I]:= 0; J:=0 end
        else ONELINE[I]:= -1; (* blank *)
        J:= J + 1
        end;
 
      (* SIMILARITY SEARCH *)
      writeln('Search begins...');
      HEADER;
      while LEFT <= RIGHTLIM do begin
        HORAXIS(LEFT,RIGHT);
        MAKETABLE(LEFT,RIGHT);
        INDEX:=0; K:= STARTY mod 10;
        if K = 0 then THISLINE:= TENLINE else THISLINE:= ONELINE;
        for POSY:= STARTY to FINISHY do begin
          N1:=SEQY[POSY-1]; N2:=SEQY[POSY]; N3:=SEQY[POSY+1];
          if ((N1<N) and (N2<N) and (N3<N)) then
            SEARCHFOR(TRIPLET[N1,N2,N3]);
          INDEX:= INDEX+1;
          if INDEX = COMPFACT then begin
            if K < 9 then begin write(OUTFILE,' ':6); NEXTLINE:=1 end
            else if K=9 then begin write(OUTFILE,' ':6); NEXTLINE:=10 end
            else begin write(OUTFILE,POSY:6);NEXTLINE:=1;K:=0 end;
            if COMPFACT = 1 then write(OUTFILE,NUCHAR[SEQY[POSY]])
                            else write(OUTFILE,' ');
            for I:= 1 to LINESIZE do write(OUTFILE,PC[THISLINE[I]]);
            writeln(OUTFILE);
            if NEXTLINE = 1 then THISLINE:= ONELINE
            else THISLINE:= TENLINE;
            INDEX:=0; K:= K + 1
            end
          end;
          RIDOF;
          LEFT:= RIGHT+1; RIGHT:= RIGHT + (LINESIZE*COMPFACT);
          if RIGHT > RIGHTLIM then RIGHT:= RIGHTLIM
          end;
       writeln(OUTFILE); writeln(OUTFILE)
    end; (* QUICKSEARCH *)
 
 
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
(* END MODULE STARTUP         VERSION= 'SUNMODS     Version  6/26/01'; *)
      writeln;
 
      (* Read in X-axis sequence *)
      NEWSEQ(SFILEX,XFN,SEQX,LENX,STARTX,FINISHX,NAMEX,CIRCX,'X');
      
      (* Read in Y-axis sequence *)
      NEWSEQ(SFILEY,YFN,SEQY,LENY,STARTY,FINISHY,NAMEY,CIRCY,'Y');
       
      (* Open output file. *)
      writeln('Type output filename:');
      GETFILE(OUTFILE,'O',OFN);
 
      INITIALIZE;
 
      (* Initialize horizontal output line *)
      with HLINE do begin
        for J:= 1 to MAXLINE do STR[J]:='_';
        LEN:=MAXLINE
        end;
             
      (* MAIN LOOP *)
      repeat
        writeln;
        WRITELINE(output,HLINE,80);writeln;
        writeln('D3HOM   ','MAIN MENU':30);
        WRITELINE(output,HLINE,80);writeln;
        write('X-axis file:       ');WRITELINE(output,XFN,XFN.LEN);writeln;
        write('Y-axis file:       ');WRITELINE(output,YFN,YFN.LEN);writeln;
        write('Output file:       ');WRITELINE(output,OFN,OFN.LEN);writeln;
        WRITELINE(output,HLINE,80);writeln;
        writeln(' ':20,'1) Read in a new X-axis sequence');
        writeln(' ':20,'2) Read in a new Y-axis sequence');
        writeln(' ':20,'3) Open a new output file');
        writeln(' ':20,'4) Change parameters');
        writeln(' ':20,'5) Compare sequences and write output screen');
        writeln(' ':20,'6) Compare sequences and write output to file');
        WRITELINE(output,HLINE,80);writeln;
        writeln('Type the number of your choice  (0 to quit program)');
        GETREAL(TEMP,0,6);readln;
        CHOICE:=round(TEMP);

        case CHOICE of
          0:;
          1:NEWSEQ(SFILEX,XFN,SEQX,LENX,STARTX,FINISHX,NAMEX,CIRCX,'X');
          2:NEWSEQ(SFILEY,YFN,SEQY,LENY,STARTY,FINISHY,NAMEY,CIRCY,'Y');
          3:begin
(*!!!*)       CLOSE(OUTFILE); 
              writeln('Type output filename:');
              GETFILE(OUTFILE,'O',OFN);
            end;
          4:PARAMETERS;
          5:begin
              CALCSCORES;
              QUICKSEARCH(output);
              write('Press RETURN to continue');
              readln
            end;
          6:begin
              CALCSCORES;
              QUICKSEARCH(OUTFILE)
            end
          end (* case *)
        until CHOICE = 0;
(*!!!*)  CLOSE(OUTFILE) 
    end. (* D3HOM  *)
