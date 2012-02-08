(*$DEBUG- *)
  (* ********************************************************  *)
  (*                                                           *)
  (*   P2HOM     Version  5/13/91, Standard Pascal             *)
  (*             Brian Fristensky                              *)
  (*             Dept. of Plant Science                        *)
  (*             University of Manitoba                        *)
  (*             Winnipeg, MB R3T 2N2  CANADA                  *)
  (*                                                           *)
  (* Copyright (c) 1984,1986,1988,1990 by Brian Fristensky.    *)
  (* !!! in comment indicates feature which may need change.   *)
  (*  *******************************************************  *)
 
  (*!!!*)program P2HOM(input, output (*,SFILE1,SFILE2,OUTFILE*));
 (*!!!*)(* Some Pascals require file parameters in program heading.*)
 
  const MAXSEQ  = 9000;
(* BEGIN MODULE REALLIMITS *)
(*!!!*)MAXREAL = 1.7E38;
(* END MODULE REALLIMITS         VERSION= 'SUNMODS     Version  9/12/91'; *)
        MAXWORD = 25;
        MAXRANGE= 30;
        MAXLINE =  150;
        VERSION = 'P2HOM         Version  5/13/91 ';
 
  type AMINOACID = (GLY,ALA,VALINE,LEU,ILE,MET,PHE,PRO,SER,THR,CYS,ASN,
                    GLN,TYR,TRP,ASP,GLU,HIS,LYS,ARG,ASX,GLX,TERM,UNKX,UNKY);
                    
       PROTEIN   = array[-MAXRANGE..MAXSEQ] of AMINOACID;
       NP        = ^NODE;
       NODE      = record POS:integer;
                            NEXT:NP
                            end;
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
 
  var SFILEX,SFILEY,OUTFILE: text;
      XFN,YFN,OFN,HLINE:LINE;
      SEQX,SEQY    : PROTEIN;
      LENX,LENY    : integer;
      NAMEX,NAMEY  : WORD;
      DIPEPTIDE    : array[GLY..ARG,GLY..ARG] of NP;
      FREELIST     : NP;
      AACHAR       : array[AMINOACID] of char;
      STARTX,FINISHX,
      STARTY,FINISHY,
      HOMRANGE,MINPER,COMPFACT,LINESIZE,
      MAXSCORE,TWOMATCH,MINSCORE: integer;
      SCALEFAC,MAXFACTOR,TEMP:real;
      SUBSCORE     : array[1..MAXRANGE] of integer;
      PC: array[-1..100] of char; (* Printing characters for graph *)
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
(* END MODULE INPLINE         VERSION= 'SUNMODS     Version  9/12/91'; *)
 
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
(* END MODULE GETFILE         VERSION= 'SUNMODS     Version  9/12/91'; *)
 
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
(* END MODULE INPWORD         VERSION= 'SUNMODS     Version  9/12/91'; *)

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
(* END MODULE GETREAL         VERSION= 'SUNMODS     Version  9/12/91'; *)
 
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
(* END MODULE READPRO         VERSION= 'SUNMODS     Version  9/12/91'; *)

  (*****************************************************)
  (* Initialization procedures.                        *)
  (*****************************************************)
  procedure INITIALIZE;
    var A1,A2:AMINOACID;
    begin
      (* Initialize the PC array which holds the       *)
      (* characers which symbolize percent identity.   *)
      PC[ -1]:=' ';PC[  0]:='.';
      PC[  5]:='w';PC[  6]:='v';PC[  7]:='v';PC[  8]:='u';PC[  9]:='u';
      PC[ 10]:='t';PC[ 11]:='t';PC[ 12]:='s';PC[ 13]:='s';PC[ 14]:='r';
      PC[ 15]:='r';PC[ 16]:='q';PC[ 17]:='q';PC[ 18]:='p';PC[ 19]:='p';
      PC[ 20]:='o';PC[ 21]:='o';PC[ 22]:='n';PC[ 23]:='n';PC[ 24]:='m';
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
      (*   Initialize AACHAR array which holds the            *)
      (*   character values of the AMINOACIDs.                *)
      AACHAR[GLY]:='G'; AACHAR[ALA]:='A';AACHAR[VALINE]:='V'; AACHAR[LEU]:='L';
      AACHAR[ILE]:='I'; AACHAR[MET]:='M';AACHAR[PHE]:='F'; AACHAR[PRO]:='P';
      AACHAR[SER]:='S'; AACHAR[THR]:='T';AACHAR[CYS]:='C'; AACHAR[ASN]:='N';
      AACHAR[GLN]:='Q'; AACHAR[TYR]:='Y';AACHAR[TRP]:='W'; AACHAR[ASP]:='D';
      AACHAR[GLU]:='E'; AACHAR[HIS]:='H';AACHAR[LYS]:='K'; AACHAR[ARG]:='R';
      AACHAR[ASX]:='B'; AACHAR[GLX]:='Z';AACHAR[TERM]:='*';
      AACHAR[UNKX]:='X';AACHAR[UNKY]:='X';
 
      (*Set default values for parameters *)
      HOMRANGE:= 10; SCALEFAC:=0.90;
      MINPER:=50; COMPFACT:=10; LINESIZE:= 70;
 
      (* Initialize DIPEPTIDE table to nil *)
      for A1:= GLY to ARG do
        for A2:= GLY to ARG do DIPEPTIDE[A1,A2]:= nil
      end; (* INITIALIZE *)
 
    (********************************************************)
    (* Open a sequence file and read in sequence.           *)
    (********************************************************)
    procedure NEWSEQ(var SFILE:text; var FN:LINE; var PR:PROTEIN;
                     var LEN,START,FINISH:integer; var NAME:WORD; AXIS:char);
      var ANSWER:char;
          J:integer;
      procedure SEQENDS(var PR:PROTEIN; LEN:integer; UNKNOWN:AMINOACID);
        var J:integer;
        begin
          (* Initialize negative and last parts of the sequence*)
          for J:= -MAXRANGE to 0 do PR[J]:= UNKNOWN;
          for J:= LEN+1 to LEN+MAXRANGE do PR[J]:= UNKNOWN
        end; (* SEQENDS *)
      begin
        writeln('Enter filename for sequence on ',AXIS,'-axis:');
        GETFILE(SFILE,'I',FN);
        NAME.LEN:=0;
        READPRO(SFILE,PR,LEN,NAME);
        case AXIS of
         'X':SEQENDS(PR,LEN,UNKX);
         'Y':begin
               SEQENDS(PR,LEN,UNKY);
               for J:= 1 to LEN do if PR[J]=UNKX then PR[J]:=UNKY
             end  
          end;
        if NAME.LEN=0 then begin
          writeln('Type name for ',AXIS,'-axis sequence to appear on output:');
          INPWORD(NAME);readln
          end;
        START:=1; FINISH:=LEN
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
 
  (* Display  parameters on screen *)
  procedure DISPLAY;
     begin
       page(output);
       write('X-axis sequence: ');WRITEWORD(output,NAMEX,20);
       writeln(' Length: ':21,LENX:10,' aa');
       write('Y-axis sequence: ');WRITEWORD(output,NAMEY,20);
       writeln(' Length: ':21,LENY:10,' aa');
       WRITELINE(output,HLINE,80);writeln;
       writeln(' ':12,'Parameter   Description/Response                 Value');
       WRITELINE(output,HLINE,80);writeln;
       writeln(' ':12,' 1)STARTX   first amino acid position in SEQX   ',
               STARTX:6);
       writeln(' ':12,' 2)FINISHX  last  amino acid position in SEQX   ',
               FINISHX:6);
       writeln(' ':12,' 3)STARTY   first amino acid position in SEQY   ',
               STARTY:6);
       writeln(' ':12,' 4)FINISHY  last  amino acid position in SEQY   ',
               FINISHY:6);
       writeln(' ':12,' 5)HOMRANGE dist.from central a.a. in a match.  ',
               HOMRANGE:6);
       writeln(' ':12,' 6)SCALEFAC scale factor for exponential curve  ',
               SCALEFAC:6:2);
       writeln(' ':12,' 7)MINPER   minimum percent similarity printed  ',
               MINPER:6);
       writeln(' ':12,' 8)COMPFACT graph compression factor            ',
               COMPFACT:6);
       writeln(' ':12,' 9)LINESIZE width of output line (ex. 70,120)   ',
               LINESIZE:6);
       WRITELINE(output,HLINE,80);writeln
     end; (* DISPLAY *)
 
      begin
        (* Prompt user for new parameter values *)
        repeat
          DISPLAY;
          writeln('Type number of parameter you wish to change',
                  ' (0 to continue)');
          GETREAL(TEMP,0,9);readln;
          RESPONSE:= round(TEMP);
          if RESPONSE in [1..9] then
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
              7:GETNUMBER(MINPER,'MINPER    ',5,100);
              8:GETNUMBER(COMPFACT,'COMPFACT  ',1,100);
              9:GETNUMBER(LINESIZE,'LINESIZE  ',40,MAXLINE)
             end
         until RESPONSE= 0
     end; (* PARAMETERS *)
 
  (****************************************************************)
  (* Calculate a table of scores for matches at each distance from*)
  (* the central dipeptide in a local homology.                   *)
  (****************************************************************)
  procedure CALCSCORES;
    var I,V:integer;
        S:real;
    begin
      V:=100;
      SUBSCORE[1]:= V; MAXSCORE:= 2*V; S:=SCALEFAC;
      for I:= 2 to HOMRANGE do begin
        SUBSCORE[I]:= round(V*S);
        MAXSCORE:= MAXSCORE + (2*SUBSCORE[I]);
        S:= S * SCALEFAC
        end;
      TWOMATCH:= 2*SUBSCORE[1];
      MINSCORE:= round((MINPER/100)*MAXSCORE);
      MAXFACTOR:= 100/MAXSCORE
    end; (*CALCSCORES*)
 
  (* ***************************************************************** *)
  (*  Find similarities (2*HOMRANGE) long with MINPER or better match. *)
  (* ***************************************************************** *)
  procedure QUICKSEARCH(var OUTFILE:text);
 
    var POSX,POSY,RIGHT,RIGHTLIM,LEFT,I,J,K,INDEX,NEXTLINE:integer;
        ONELINE,TENLINE,THISLINE: array[1..MAXLINE] of -1..100;
        A1,A2:AMINOACID;
         
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
        NUMBER:= LEFT-1;
        for I:= 1 to NUMPRINTED do begin
          NUMBER:= NUMBER + DELTA;
          write(OUTFILE,NUMBER:10)
          end;
        writeln(OUTFILE);
 
        (* Write sequence if COMPFACT = 1*)
        if COMPFACT = 1 then begin
          write(OUTFILE,' ':7);
          for I:= LEFT to RIGHT do write(OUTFILE,AACHAR[SEQX[I]])
          end;
        writeln(OUTFILE)
      end; (*HORAXIS*)
 
   (*  Make a table of locations of dipeptides in SEQX.  Each dipeptide  *)
   (*  has a stack of nodes, each of which holds a location in SEQX      *)
   (*  at which the dipeptide occurs.                                    *)
   procedure MAKETABLE(LEFT,RIGHT:integer);
     var I:integer;
         A1,A2: AMINOACID;
 
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
         A1:= SEQX[I]; A2:= SEQX[I+1];
         if  ((A1<=ARG) and (A2<=ARG)) then ADDNODE(DIPEPTIDE[A1,A2])
         end
    end; (*MAKETABLE*)
 
   (* At each occurrence of the dipeptide in SEQX, compare a region  *)
   (* HOMRANGE bases on either side.   If the match is good enough,  *)
   (* print a character at the corresponding point in the matrix.    *)
   procedure SEARCHFOR(N:NP);
 
     var LX,LY,RX,RY,DISTANCE,SCORE,PERCENT,X:integer;
         CH:char;
 
     begin
       while N <> nil do begin
         POSX:= N^.POS;
         SCORE:= TWOMATCH;
         LX:= POSX-1; RX:= POSX+2; LY:=POSY-1; RY:=POSY+2;
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
 
  (* Push the linked-list (if there is one) for each DIPEPTIDE to *)
  (* the top of the FREELIST and reset DIPEPTIDE to nil.          *)
  procedure RIDOF;
    var A1,A2:AMINOACID;
        HEAD,TAIL:NP;
    begin
      for A1:= GLY to ARG do
        for A2:= GLY to ARG do
          if DIPEPTIDE[A1,A2] <> nil then begin
             HEAD:= DIPEPTIDE[A1,A2];
             TAIL:= HEAD;
             while TAIL^.NEXT <> nil do TAIL:= TAIL^.NEXT;
             TAIL^.NEXT:= FREELIST;
             FREELIST:= HEAD;
             DIPEPTIDE[A1,A2]:= nil
             end
    end; (*RIDOF*)
 
    begin
      (* Initialize LEFT & RIGHT, which define the part of SEQX *)
      (* to be compared with SEQY.                              *)
      LEFT:= STARTX;
      RIGHTLIM:= FINISHX;
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
          A1:=SEQY[POSY]; A2:=SEQY[POSY+1];
          if (A1<=ARG) and (A2<=ARG) then
            SEARCHFOR(DIPEPTIDE[A1,A2]);
          INDEX:= INDEX+1;
          if INDEX = COMPFACT then begin
            if K < 9 then begin write(OUTFILE,' ':6); NEXTLINE:=1 end
            else if K = 9 then begin write(OUTFILE,' ':6);NEXTLINE:=10 end
            else begin write(OUTFILE,POSY:6);NEXTLINE:=1;K:=0 end;
            if COMPFACT = 1 then write(OUTFILE,AACHAR[SEQY[POSY]])
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
(* END MODULE STARTUP         VERSION= 'SUNMODS     Version  9/12/91'; *)
      writeln;
 
      (* Read in X-axis sequence *)
      NEWSEQ(SFILEX,XFN,SEQX,LENX,STARTX,FINISHX,NAMEX,'X');
      
      (* Read in Y-axis sequence *)
      NEWSEQ(SFILEY,YFN,SEQY,LENY,STARTY,FINISHY,NAMEY,'Y');
      
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
        writeln('P2HOM   ','MAIN MENU':30);
        WRITELINE(output,HLINE,80);writeln;
        write('X-axis file:       ');WRITELINE(output,XFN,XFN.LEN);writeln;
        write('Y-axis file:       ');WRITELINE(output,YFN,YFN.LEN);writeln;
        write('Output file:       ');WRITELINE(output,OFN,OFN.LEN);writeln;
        WRITELINE(output,HLINE,80);writeln;
        writeln(' ':20,'1) Read in a new X-axis sequence');
        writeln(' ':20,'2) Read in a new Y-axis sequence');
        writeln(' ':20,'3) Open a new output file');
        writeln(' ':20,'4) Change parameters');
        writeln(' ':20,'5) Compare sequences and write output to screen');
        writeln(' ':20,'6) Compare sequences and write output to file');
        WRITELINE(output,HLINE,80);writeln;
        writeln('Type the number of your choice  (0 to quit program)');
        GETREAL(TEMP,0,6);readln;
        CHOICE:=round(TEMP);

        case CHOICE of
          0:;
          1:NEWSEQ(SFILEX,XFN,SEQX,LENX,STARTX,FINISHX,NAMEX,'X');
          2:NEWSEQ(SFILEY,YFN,SEQY,LENY,STARTY,FINISHY,NAMEY,'Y');
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
              QUICKSEARCH(OUTFILE);
            end
          end (* case *)
        until CHOICE = 0;
(*!!!*)  CLOSE(OUTFILE) 
    end. (* P2HOM  *)

