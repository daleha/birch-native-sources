  (* ********************************************************  *)
  (*                                                           *)
  (*   TESTCODE  Version   8/13/2001 Standard Pascal           *)
  (*             Brian Fristensky                              *)
  (*             Dept. of Plant Science                        *)
  (*             University of Manitoba                        *)
  (*             Winnipeg, MB  R3T 2N2  CANADA                 *)
  (*                                                           *)
  (*  Copyright (c) 1984-1990  by Brian Fristensky             *)
  (*  !!! in comment indicates feature which may need change.  *)
  (*  *******************************************************  *)
  (* REVISION HISTORY
  Aug. 13, 2001 Rebuilt using new Readseq procedure which reads
  new GenBank format.
  
*)
 
(*!!!*)program TESTCODE(input, output (*,INFILE,OUTFILE*));
(*!!! Some Pascals require file parameters in program heading *)
 
  const MAXSEQ   =  500000;
        MINCODONS =  30;
        MAXLINE   =  70;
        MAXWORD   =  25;
        LINESIZE  =  70;
        VERSION   = 'TESTCODE          Version   8/13/2001';
 
  type (* NUCLEOTIDE types *)
       NUCLEOTIDE  = (T,C,A,G,N,R,Y,M,W,S,K,D,H,V,B,Z);
       NA          = array[T..N] of NUCLEOTIDE;
       RA          = array[T..G] of real;
       SEQUENCE    = array[1..MAXSEQ] of NUCLEOTIDE;
 
       (* table for probability values *)
       PARAM   =  record
                    LIMIT: real;
                    PROB : array[T..G] of real
                    end;
       TABLE   =  array[1..10] of PARAM;
 
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
 
  var INFILE,OUTFILE  : text;
      IFN,OFN,HLINE   : LINE;
      SEQ             : SEQUENCE;
      SEQLEN          : integer; (* sequence length *)
      SEQNAME         : WORD;
      POSTABLE,CONTABLE: TABLE;
      INDVAL,CODPROB:array[1..10] of real;
      INPSTRAND,COMPSTRAND :NA;
      TITLE:LINE;
      (* GLOBAL PARAMETERS *)
      START,FINISH,WINDOW,SKIP,I,CHOICE:integer;
      WHICH,FORMAT: char;
      CIRCULAR       : boolean;

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
(* END MODULE GETINTEGER         VERSION= 'SUNMODS     Version  6/26/01'; *)
 
  (******************************************************)
  (* Convert ambiguous nucleotides to N's.              *)
  (******************************************************)
  procedure AMBIGTON(var S:SEQUENCE;SEQLEN:integer);
    var I:integer;
    begin
      for I:= 1 to SEQLEN do
        if S[I] > N then S[I]:=N
    end; (* AMBIGTON *)

    (******************************************************)
    (* Initialize the probability tables.                 *)
    (******************************************************)
    procedure INITTABLES;
      (* Initialize the Position table *)
      procedure INITPOST;
      begin
      with POSTABLE[ 1] do begin
        LIMIT:=0.0;PROB[A]:=0.22;PROB[C]:=0.23;PROB[G]:=0.08;PROB[T]:=0.09 end;
      with POSTABLE[ 2] do begin
        LIMIT:=1.1;PROB[A]:=0.20;PROB[C]:=0.30;PROB[G]:=0.08;PROB[T]:=0.09 end;
      with POSTABLE[ 3] do begin
        LIMIT:=1.2;PROB[A]:=0.34;PROB[C]:=0.33;PROB[G]:=0.16;PROB[T]:=0.20 end;
      with POSTABLE[ 4] do begin
        LIMIT:=1.3;PROB[A]:=0.45;PROB[C]:=0.51;PROB[G]:=0.27;PROB[T]:=0.54 end;
      with POSTABLE[ 5] do begin
        LIMIT:=1.4;PROB[A]:=0.68;PROB[C]:=0.48;PROB[G]:=0.48;PROB[T]:=0.44 end;
      with POSTABLE[ 6] do begin
        LIMIT:=1.5;PROB[A]:=0.58;PROB[C]:=0.66;PROB[G]:=0.53;PROB[T]:=0.69 end;
      with POSTABLE[ 7] do begin
        LIMIT:=1.6;PROB[A]:=0.93;PROB[C]:=0.81;PROB[G]:=0.64;PROB[T]:=0.68 end;
      with POSTABLE[ 8] do begin
        LIMIT:=1.7;PROB[A]:=0.84;PROB[C]:=0.70;PROB[G]:=0.74;PROB[T]:=0.91 end;
      with POSTABLE[ 9] do begin
        LIMIT:=1.8;PROB[A]:=0.68;PROB[C]:=0.70;PROB[G]:=0.88;PROB[T]:=0.97 end;
      with POSTABLE[10] do begin
        LIMIT:=1.9;PROB[A]:=0.94;PROB[C]:=0.80;PROB[G]:=0.90;PROB[T]:=0.97 end
      end; (* INITPOST *)
 
      (* Initialize the Content table *)
      procedure INITCONT;
      begin
      with CONTABLE[ 1] do begin
        LIMIT:=0.00;PROB[A]:=0.21;PROB[C]:=0.31;PROB[G]:=0.29;PROB[T]:=0.58 end;
      with CONTABLE[ 2] do begin
        LIMIT:=0.17;PROB[A]:=0.81;PROB[C]:=0.39;PROB[G]:=0.33;PROB[T]:=0.51 end;
      with CONTABLE[ 3] do begin
        LIMIT:=0.19;PROB[A]:=0.65;PROB[C]:=0.44;PROB[G]:=0.41;PROB[T]:=0.69 end;
      with CONTABLE[ 4] do begin
        LIMIT:=0.21;PROB[A]:=0.67;PROB[C]:=0.43;PROB[G]:=0.41;PROB[T]:=0.56 end;
      with CONTABLE[ 5] do begin
        LIMIT:=0.23;PROB[A]:=0.49;PROB[C]:=0.59;PROB[G]:=0.73;PROB[T]:=0.75 end;
      with CONTABLE[ 6] do begin
        LIMIT:=0.25;PROB[A]:=0.62;PROB[C]:=0.59;PROB[G]:=0.64;PROB[T]:=0.55 end;
      with CONTABLE[ 7] do begin
        LIMIT:=0.27;PROB[A]:=0.55;PROB[C]:=0.64;PROB[G]:=0.64;PROB[T]:=0.40 end;
      with CONTABLE[ 8] do begin
        LIMIT:=0.29;PROB[A]:=0.44;PROB[C]:=0.51;PROB[G]:=0.47;PROB[T]:=0.39 end;
      with CONTABLE[ 9] do begin
        LIMIT:=0.31;PROB[A]:=0.49;PROB[C]:=0.64;PROB[G]:=0.54;PROB[T]:=0.24 end;
      with CONTABLE[10] do begin
        LIMIT:=0.33;PROB[A]:=0.28;PROB[C]:=0.82;PROB[G]:=0.40;PROB[T]:=0.28 end;
      end; (* INITCONT *)
 
      (* Initialize the indicator arrays *)
      procedure INITIND;
      begin
      INDVAL[ 1]:=0.00; CODPROB[ 1]:=0.00;
      INDVAL[ 2]:=0.43; CODPROB[ 2]:=0.04;
      INDVAL[ 3]:=0.53; CODPROB[ 3]:=0.07;
      INDVAL[ 4]:=0.64; CODPROB[ 4]:=0.29;
      INDVAL[ 5]:=0.74; CODPROB[ 5]:=0.40;
      INDVAL[ 6]:=0.84; CODPROB[ 6]:=0.77;
      INDVAL[ 7]:=0.95; CODPROB[ 7]:=0.92;
      INDVAL[ 8]:=1.05; CODPROB[ 8]:=0.98;
      INDVAL[ 9]:=1.16; CODPROB[ 9]:=1.00;
      INDVAL[10]:=1.26; CODPROB[10]:=1.00
      end; (* INITIND *)
 
    begin
      INITPOST; INITCONT; INITIND
    end; (* INITTABLES *)
 
    (* **************************************************** *)
    (* Set default values for global parameters             *)
    (* **************************************************** *)
    procedure INITPARAM;
      begin
        START:= 1;FINISH:= SEQLEN;WHICH:='I';FORMAT:='G';
        WINDOW:= 67; SKIP:=10
      end; (* INITPARAM *)
 
    (* **************************************************** *)
    (* Prompt user for parameters used by program.          *)
    (* **************************************************** *)
    procedure PARAMETERS;
      type LETTERS = packed array[1..10] of char;
           CHSET   = set of 'A'..'Z';
      var RESPONSE :integer;
 
      (*  Read an integer parameter from the console and check *)
      (*    that it is in range.                               *)
      procedure GETNUMBER(var P:integer;PNAME:LETTERS;LOW,HIGH:integer);
        begin
          writeln;
          writeln('Type new value for ',PNAME,'  (CURRENT VALUE: ',P,')');
          GETINTEGER(P,LOW,HIGH);readln
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
      begin
      write('Name: ');WRITEWORD(output,SEQNAME,20);
      case CIRCULAR of 
        true:write('Topology: CIRCULAR');
        false:write('Topology:   LINEAR')
        end;
      writeln('Length: ':13,SEQLEN:11,' nt');
      WRITELINE(output,HLINE,80);writeln;
      writeln(' ':12,'Parameter   Description/Response                 Value');
      WRITELINE(output,HLINE,80);writeln;
      writeln(' ':12,' 1)START    first nucleotide evaluated',START:16);
      writeln(' ':12,' 2)FINISH   last  nucleotide evaluated',FINISH:16);
      writeln(' ':12,' 3)WHICH    I: input strand  O: opposite strand',WHICH:7);
      writeln(' ':12,' 4)FORMAT   T:tabular output G:graphic output ',FORMAT:8);
      writeln(' ':12,' 5)WINDOW   #codons in search window',WINDOW:18);
      writeln(' ':12,' 6)SKIP     #codons to skip for each window',SKIP:11);
      WRITELINE(output,HLINE,80);
      writeln(output)
      end; (* DISPLAY *)
 
      begin
        (* Prompt user for new parameter values *)
        repeat
          page(output);
          DISPLAY;
          if WHICH='O' then begin
            writeln('Be sure START and FINISH values are appropriate for');
            writeln('WHICH=O')end;
          writeln('Type number of parameter you wish to change',
                  ' (0 to continue)');
          GETINTEGER(RESPONSE,0,6);readln;
          if RESPONSE in [1..6] then
            case RESPONSE of
               1:GETNUMBER(START,'START     ',1,SEQLEN);
               2:GETNUMBER(FINISH,'FINISH    ',1,SEQLEN);
               3:GETCHAR(WHICH,'WHICH     ',['I','O']);
               4:GETCHAR(FORMAT,'FORMAT    ',['T','G']);
               5:GETNUMBER(WINDOW,'WINDOW    ',10,SEQLEN div 3);
               6:GETNUMBER(SKIP,'SKIP      ',1,SEQLEN div 3)
             end
         until RESPONSE= 0
     end; (* PARAMETERS *)
 
    (* **************************************************** *)
    (*   Initialize INPSTRAND and COMPSTRAND, the arrays    *)
    (*   which hold the valaues of the NUCLEOTIDES.         *)
    (* **************************************************** *)
    procedure INITNUCLEOTIDES;
      begin
        INPSTRAND[A]:=A;COMPSTRAND[A]:=T;
        INPSTRAND[C]:=C;COMPSTRAND[C]:=G;
        INPSTRAND[G]:=G;COMPSTRAND[G]:=C;
        INPSTRAND[T]:=T;COMPSTRAND[T]:=A;
        INPSTRAND[N]:=N;COMPSTRAND[N]:=N
      end; (* INITNUCLEOTIDES *)
 
  (* *********************************************************** *)
  (*  Test an open reading frame specified by the user to        *)
  (*  to determine whether or not it is a protein coding region  *)
  (*  using the algorithm found in:                              *)
  (*   Fickett,James, "Recognition of protein coding regions in  *)
  (*   DNA sequences", Nucleic Acids Research 10 No.17,5303-5318.*)
  (* *********************************************************** *)
  procedure TEST(var OUTFILE:text; STRAND:NA;START:integer);
    var  I,POSITION,DISTANCE,CODONSDONE,NUMCYCLES,CYCLESDONE,
         CENTER,SENSE,S2,S3:integer;
         NUC:NUCLEOTIDE;
         FREQ :array[T..N,1..3] of integer;
         POS,CONT :RA;
         INDICATOR:real;
         OUTLINE,TEMPLATE:array[0..70] of char;

    procedure SETPARAM;
      begin
      (* Set parameters dependant on the strand *)
      case WHICH of
        'I': begin
               SENSE:= 1;
               if START < FINISH then DISTANCE:= (FINISH - START + 1) div 3
               else DISTANCE:= (SEQLEN - START + FINISH + 1) div 3
               end;
        'O': begin
               SENSE:= -1;
               if START > FINISH then DISTANCE:= (START - FINISH + 1) div 3
               else DISTANCE:=(START + (SEQLEN - FINISH + 1)) div 3
               end
         end
      end; (* SETPARAM *)
 
    (* Set up variables and print header for graph *)
    procedure SETUPGRAPH;
      var I:integer;
          INDEX:real;
      begin
        writeln(OUTFILE,' ':10,'WINDOW=',WINDOW:5,'SKIP=':15,SKIP:5);
        writeln(OUTFILE);
        writeln(OUTFILE,'NON-CODING':35,'NO OPINION':22,'CODING':14);
        writeln(OUTFILE);
        INDEX:=0.0;
        for I:= 0 to 7 do begin
          write(OUTFILE,INDEX:10:1);INDEX:=INDEX+0.2 end;
        writeln(OUTFILE);
        write(OUTFILE,' ':10);
        for I:= 1 to 7 do write(OUTFILE,'---------+');
        writeln(OUTFILE);

        (* Initialize line templates *)
        for I:= 1 to 70 do TEMPLATE[I]:=' ';
        TEMPLATE[0]:='|'; TEMPLATE[37]:='|'; TEMPLATE[47]:='|';
        OUTLINE:= TEMPLATE
      end; (* SETUPGRAPH *)

    (* Compute the next nucleotide position *)
    function NEXT(POS:integer):integer;
      begin
        POS:= POS + SENSE;
        if CIRCULAR then
          if POS > SEQLEN then NEXT:= 1
          else if POS < 1 then NEXT:= SEQLEN
          else NEXT:= POS
        else if (POS > SEQLEN) or (POS < 1) then NEXT:= 0
             else NEXT:= POS
      end; (* NEXT *)
 
    (*Print the intermediate results *)
    procedure INTERMEDIATE;
      var I:integer;
          NUC:NUCLEOTIDE;
 
      procedure PRINTARRAY(AR:RA);
        var NUC:NUCLEOTIDE;
        begin
          for NUC:= T to G do write(OUTFILE,AR[NUC]:10:2);
          writeln(OUTFILE)
        end; (* PRINTARRAY *)
 
      begin
        writeln(OUTFILE);
        writeln(OUTFILE,'Open reading frame from ',START,' to ',FINISH);
        writeln(OUTFILE,'T':20,'C':10,'A':10,'G':10);
        writeln(OUTFILE,'Pos.Freq.');
        for I:= 1 to 3 do begin
          write(OUTFILE,I:10);
          for NUC:= T to G do write(OUTFILE,FREQ[NUC,I]:10);
          writeln(OUTFILE)
          end;
        for I:= 1 to 60 do write(OUTFILE,'-');
        writeln(OUTFILE);
        write(OUTFILE,'Cont.Param');PRINTARRAY(CONT);
        write(OUTFILE,'Posn.Param');PRINTARRAY(POS)
      end; (* INTERMEDIATE *)
 
    (* Calculate the maximum of three integers *)
    function MAX(I1,I2,I3:integer):integer;
      var TEMP:integer;
      begin
        if I1 > I2 then TEMP := I1
        else TEMP:= I2;
        if I3 > TEMP then MAX:= I3
        else MAX:= TEMP
      end; (* MAX *)
 
    (* Calculate the minimum of three integers *)
    function MIN(I1,I2,I3:integer):integer;
      var TEMP:integer;
      begin
        if I1 < I2 then TEMP := I1
        else TEMP:= I2;
        if I3 < TEMP then MIN:= I3
        else MIN:= TEMP
      end; (* MIN *)
 
    (* Look up the probability of a given parameter in the table specified *)
    procedure FINDPROB(var P:real; T:TABLE);
      var I:integer;
      begin
        I:=10;
        while P < T[I].LIMIT do I:= I - 1;
        P:= T[I].PROB[NUC]
      end; (* FINDPROB *)
 
    (* Print the results *)
    procedure RESULTS;
      var NUMCHAR,J:integer;
      begin
        case FORMAT of
          'T':begin
                writeln(OUTFILE,'TESTCODE indicator: ',INDICATOR);
                writeln(OUTFILE,'Probability of coding: ',CODPROB[I]);
                writeln(OUTFILE,'Prediction: ');
                case I of
                  1,2,3,4: writeln(OUTFILE,'NONCODING');
                  5,6: writeln(OUTFILE,'NO OPINION');
                  7,8,9,10:writeln(OUTFILE,'CODING')
                  end
              end;
          'G':begin
                NUMCHAR:= round((INDICATOR/1.4) * LINESIZE);
                if NUMCHAR > LINESIZE then NUMCHAR:= LINESIZE;
                write(OUTFILE,CENTER:9);
                for J:= 1 to NUMCHAR do OUTLINE[J]:='=';
                for J:= 0 to LINESIZE do write(OUTFILE,OUTLINE[J]);
                writeln(OUTFILE)
                end
          end (* case *)
      end; (* RESULTS *)

    begin (* TEST *)
      SETPARAM;

      (* Only test regions >= MINCODONS *)
      if DISTANCE < MINCODONS then 
        writeln('Reading frame must be >= ',MINCODONS,' codons (',
          (MINCODONS * 3),' bp)')
      else begin
        (* Print the header *)
        writeln(OUTFILE); writeln(OUTFILE);
        writeln(OUTFILE,VERSION:45);
        writeln(OUTFILE);
        write(OUTFILE,' ':10);
        WRITELINE(OUTFILE,TITLE,TITLE.LEN);writeln(OUTFILE);

        if FORMAT='G' then begin
          SETUPGRAPH;
          NUMCYCLES:= (DISTANCE-WINDOW) div SKIP+1
          end
        else begin
          WINDOW:= DISTANCE;
          NUMCYCLES:=1
          end; 
        CYCLESDONE:=0;

        while CYCLESDONE < NUMCYCLES do begin
          (* Initialize the frequency table *)
          for NUC:= T to N do
            for I:= 1 to 3 do FREQ[NUC,I]:=0;
 
          (* For the reading frame specified, calculate the sums of *)
          (* Ai, Gi, Ci, Ti for positions i= 1 to 3 in the codons   *)
          POSITION:= START;
          CENTER:= WINDOW div 2;
          CODONSDONE:=0;
          while CODONSDONE < WINDOW do begin
            S2:= NEXT(POSITION); S3:= NEXT(S2);
            FREQ[STRAND[SEQ[POSITION]],1]:= FREQ[STRAND[SEQ[POSITION]],1]+1;
            FREQ[STRAND[SEQ[S2]],2]:= FREQ[STRAND[SEQ[S2]],2]+1;
            FREQ[STRAND[SEQ[S3]],3]:= FREQ[STRAND[SEQ[S3]],3]+1;
            (* CENTER determines coordinate for center of winwow *)
            CODONSDONE:= CODONSDONE+1;
            if CODONSDONE = CENTER then CENTER:= POSITION;
            POSITION:= NEXT(S3)
            end; (* CODONSDONE < WINDOW *)
 
          (* Calculate the position and content parameters *)
          (* Content =  % A,C,G,T                          *)
          (*                   MAX(A1,A2,A3)               *)
          (* A-Posn. =        -----------------     etc.   *)
          (*                   MIN(A1,A2,A3) + 1           *)
          for NUC:= T to G do begin
            CONT[NUC]:= (FREQ[NUC,1]+FREQ[NUC,2]+FREQ[NUC,3]) / (DISTANCE*3);
            POS[NUC]:=  MAX(FREQ[NUC,1],FREQ[NUC,2],FREQ[NUC,3]) /
            (MIN(FREQ[NUC,1],FREQ[NUC,2],FREQ[NUC,3])+1)
            end;
 
          if FORMAT='T' then INTERMEDIATE; (*Print intermediate results *)
 
          (* Look up the probabilities in Table 1.  Substitute the eight  *)
          (* parameters with their corresponding probabilities            *)
          for NUC:= T to G do begin
            FINDPROB(POS[NUC],POSTABLE);
            FINDPROB(CONT[NUC],CONTABLE)
            end;
 
          (* Multiply by weighting factors in Table 2 *)
          POS[A]:= POS[A] * 0.26; CONT[A]:=CONT[A] * 0.11;
          POS[C]:= POS[C] * 0.18; CONT[C]:=CONT[C] * 0.12;
          POS[G]:= POS[G] * 0.31; CONT[G]:=CONT[G] * 0.15;
          POS[T]:= POS[T] * 0.33; CONT[T]:=CONT[T] * 0.14;
 
          (* Calculate the value of the TESTCODE indicator as *)
          (* the sum of the weighted products.                *)
          INDICATOR:=0;
          for NUC:= T to G do INDICATOR:= INDICATOR + POS[NUC] + CONT[NUC];
          I:= 10;
          while INDICATOR < INDVAL[I] do I:= I - 1;
          
          (* Print the results *)    
          RESULTS;

          CYCLESDONE:= CYCLESDONE + 1;
          OUTLINE:=TEMPLATE;
          for I:= 1 to SKIP*3 do START:= NEXT(START)
          end (* while CYCLESDONE < NUMCYCLES *)
        end
      end; (* TEST *)

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
      writeln('Type sequence filename:');
      GETFILE(INFILE,'I',IFN);
      READSEQ(INFILE,SEQ,SEQLEN,SEQNAME,CIRCULAR);
      writeln('Type output filename:');
      GETFILE(OUTFILE,'O',OFN);
      writeln(OUTFILE);
      writeln;
      AMBIGTON(SEQ,SEQLEN);
      INITPARAM;
      INITTABLES;
      INITNUCLEOTIDES;

      (* Initialize horizontal output line *)
      with HLINE do begin
        for I:= 1 to MAXLINE do STR[I]:='_';
        LEN:=MAXLINE
        end;
             
      (* MAIN LOOP *)
      repeat
        writeln;
        WRITELINE(output,HLINE,80);writeln;
        writeln('TESTCODE','MAIN MENU':30);
        WRITELINE(output,HLINE,80);writeln;
        write('Input file:        ');WRITELINE(output,IFN,IFN.LEN);writeln;
        write('Output file:       ');WRITELINE(output,OFN,OFN.LEN);writeln;
        write('Title:             ');WRITELINE(output,TITLE,TITLE.LEN);writeln;
        WRITELINE(output,HLINE,80);writeln;
        writeln(' ':20,'1) Read in a new sequence');
        writeln(' ':20,'2) Open a new output file');
        writeln(' ':20,'3) Type in a title line for output');
        writeln(' ':20,'4) Change parameters');
        writeln(' ':20,'5) Search sequence (output to screen)');
        writeln(' ':20,'6) Search sequence (output to file)');
        WRITELINE(output,HLINE,80);writeln;
        writeln('Type the number of your choice  (0 to quit program)');
        GETINTEGER(CHOICE,0,6);readln;

        case CHOICE of
          0:;
          1:begin
              writeln('Enter sequence filename:');
              GETFILE(INFILE,'I',IFN);
              READSEQ(INFILE,SEQ,SEQLEN,SEQNAME,CIRCULAR);
              AMBIGTON(SEQ,SEQLEN);
              INITPARAM
            end;
          2:begin
(*!!!*)       CLOSE(OUTFILE); 
              writeln('Type output filename:');
              GETFILE(OUTFILE,'O',OFN);
            end;
          3:begin
              writeln('Type a title to appear in output (<RETURN> for blank)');
              INPLINE(TITLE)
            end;
          4:PARAMETERS;
          5:begin
              case WHICH of
                'I':TEST(output,INPSTRAND,START);
                'O':TEST(output,COMPSTRAND,START)
                end;
              writeln(output);
              write('Press RETURN to continue');
              readln
            end;
          6:begin
              case WHICH of
                'I':TEST(OUTFILE,INPSTRAND,START);
                'O':TEST(OUTFILE,COMPSTRAND,START)
                end;
              writeln(OUTFILE);writeln(OUTFILE);
              write('Press RETURN to continue');
              readln
            end
          end (* case *)
        until CHOICE = 0;
(*!!!*) CLOSE(OUTFILE) 
    end. (* TESTCODE *)
 
