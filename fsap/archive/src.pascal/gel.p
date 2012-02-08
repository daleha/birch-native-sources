(*$DEBUG-*)
  (* ********************************************************  *)
  (*                                                           *)
  (*   GEL       Version  5/13/91  Standard Pascal             *)
  (*             Brian Fristensky                              *)
  (*             Dept. of Plant Science                        *)
  (*             University of Manitoba                        *)
  (*             Winnipeg, MB R3T 2N2  CANADA                  *)
  (*                                                           *)
  (*  Copyright (c) l984-1990           by  Brian Fristensky   *)
  (*                                                           *)
  (*  !!! indicates feature which may need change              *)
  (*  *******************************************************  *)
  (*  Given the lengths and mobilities of a set of molecular   *)
  (*  weight standards on a gel, this program will estimate    *)
  (*  the sizes of unknown fragments whose mobilities are      *)
  (*  known, using the least squares approach found in:        *)
  (*                                                           *)
  (*  Schaffer & Sederoff, ANAL.BIOCHEM 115,p113-122 (1981)    *)
  (*************************************************************)
 
(*!!!*)program GEL(input, output (*,OUTFILE*));
(*!!! Some Pascals require file parameters in program heading*)
 
  const MAXFRAGS= 50;
        MAXLINE = 150;
        MAXWORD = 20;
(* BEGIN MODULE REALLIMITS *)
(*!!!*)MAXREAL = 1.7E38;
(* END MODULE REALLIMITS         VERSION= 'SUNMODS     Version  9/12/91'; *)
        VERSION = 'GEL            Version   5/13/91';
 
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
 
  var OUTFILE         : text;
      OFN,HLINE:LINE;
      ANSWER: char;
 
      (* GLOBAL PARAMETERS *)
      L0, M0, CCAP, SD, SC: real;
      L,M,PREDLEN,D,C : array[1..MAXFRAGS] of real;
      I,N:integer;
 
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
 
  (*********************************************************)
  (*  Read in standard fragments from console.             *)
  (*********************************************************)
  procedure READSTANDARDS;
 
  var LENGTH,MOBILITY: real;
  begin
    page(output);
    writeln('Type in the size of each fragment and distances migrated,');
    writeln('one fragment at a time. Press <RETURN> after each fragment.');
    writeln('Type 0 0  <RETURN> when finished');
    writeln(' Length    Mobility');
    writeln('--------+------------');
    GETREAL(LENGTH,0,MAXREAL);
    GETREAL(MOBILITY,0,MAXREAL);
    readln;
    N:=0;
    while (not((LENGTH=0) or (MOBILITY=0))) and (N < MAXFRAGS) do begin
       N:= N+1;
       L[N]:= LENGTH;
       M[N]:= MOBILITY;
       GETREAL(LENGTH,0,MAXREAL);
       GETREAL(MOBILITY,0,MAXREAL);
       readln
       end
     end; (* READSTANDARDS *)
 
   (*****************************************************)
   (* Display fragment lengths and mobilities of screen *)
   (*****************************************************)
   procedure DISPLAYFRAGS;
     var I:integer;
     begin
       page(output);
       WRITELINE(output,HLINE,80);writeln;
       writeln('Standard fragments:      Length  Mobility');
       for I:= 1 to N do writeln(I:20,')',L[I]:10:3,M[I]:10:3);
       WRITELINE(output,HLINE,80);writeln
     end; (* DISPLAYFRAGS *)
 
   (****************************************************)
   (* Let user edit the standards.                     *)
   (****************************************************)
   procedure EDIT;
     var P:integer;
         TEMP:real;
     begin
       if N > 0 then
         repeat
           DISPLAYFRAGS;
           writeln('Type number of a fragment you wish to change,');
           writeln('or 0 if all are correct.');
           GETREAL(TEMP,0,N);readln;
           P:= round(TEMP);
 
           if P > 0 then
             repeat
               writeln('Enter the correct length and mobility:');
               GETREAL(L[P],0,MAXREAL);
               GETREAL(M[P],0,MAXREAL);
               readln
             until (L[P] > 0) and (M[P] > 0)
         until P=0
       else writeln('There are  0 fragments.  Type in new standards.')
     end;  (* EDIT *)
 
 
  (*************************************)
  (* Let user add standard fragments.  *)
  (*************************************)
  procedure ADDFRAGS;
    var I,POSITION: integer;
        TEMP:real;
    begin
      if N > 0 then
        repeat
          if N < MAXFRAGS then begin
            DISPLAYFRAGS;
            writeln('Add a fragment before which number? (0 to quit) ');
            GETREAL(TEMP,0,N+1);readln;
            POSITION:= round(TEMP);
            if POSITION > 0 then begin
            (* Push up the fragment stack to make room for a new fragment *)
              for I:= N+1 downto POSITION+1 do begin
                L[I]:= L[I-1]; M[I]:= M[I-1]
                end;
              (* Add a new fragment *)
              repeat
                writeln('Type new fragment length and mobility');
                GETREAL(L[POSITION],0,MAXREAL);
                GETREAL(M[POSITION],0,MAXREAL);
                readln
              until (L[POSITION] > 0) and (M[POSITION] > 0);
              N:= N + 1
              end
            end
          else begin writeln('!!! Too many fragments !!!'); POSITION:=0 end
        until POSITION = 0
      else writeln('There are  0 fragments.  Type in new standards.')
    end; (*ADDFRAGS*)
 
  (****************************************)
  (* Let user delete standard fragments.  *)
  (****************************************)
  procedure DELETEFRAGS;
    var I,POSITION: integer;
        TEMP:real;
    begin
      if N > 0 then
        repeat
          DISPLAYFRAGS;
          writeln('Delete which fragment? (0 to quit)');
          GETREAL(TEMP,0,N);readln;
          POSITION:= round(TEMP);
          if POSITION > 0 then begin
            for I:= POSITION to N-1 do begin
              L[I]:= L[I+1]; M[I]:= M[I+1]
              end;
            N:= N-1
            end
        until POSITION = 0
      else writeln('!!! There are  0 fragments.  Type in new standards!!! ')
    end; (* DELETEFRAGS*)

  (**********************************************************************)
  (* Calculate parameters L0,M0,CCAP,SD,& SC  as described in the ref.  *)
  (**********************************************************************)
  procedure CALCULATE(var L0,M0,CCAP,SD,SC:real);

    var  MBAR,LBAR,MLBAR,    (* means of M,L, & M*L *)
         CSSM,CSSL,CSCPML,       (* corrected sums of squares of M,L, & M*L*)
         CSPMLL,CSPMLM,DELTA,
         SCi,SCiS,SDi,SDiS:real;  (* Sums & sums of sqrs of Ci & Di *)
         PROD,DWT,DDIST,DPROD : array[1..MAXFRAGS] of real;
         I: integer;
 
    begin
      (* Calculate means *)
      MBAR:=0;LBAR:=0;MLBAR:=0;
      for I:= 1 to N do begin
        MBAR:= MBAR + M[I];
        LBAR:= LBAR + L[I];
        PROD[I]:= (M[I]*L[I]);
        MLBAR:= MLBAR + PROD[I]
        end;
      MBAR:= MBAR/N; LBAR:=LBAR/N; MLBAR:= MLBAR/N;
 
      (* Calculate deviations *)
      for I:= 1 to N do begin
        DWT[I]:= L[I] - LBAR;
        DDIST[I]:= M[I] - MBAR;
        DPROD[I]:= PROD[I] - MLBAR
        end;
 
      (* Calculate intermediate values *)
      CSSM:=0;CSSL:=0;CSCPML:=0;CSPMLL:=0;CSPMLM:=0;
      for I:= 1 to N do begin
        CSSL:= CSSL + sqr(DWT[I]);
        CSSM:= CSSM + sqr(DDIST[I]);
        CSCPML:= CSCPML + (DWT[I] * DDIST[I]);
        CSPMLL:= CSPMLL + (DPROD[I] * DWT[I]);
        CSPMLM:= CSPMLM + (DPROD[I] * DDIST[I])
        end;
 
      (* Calculate L0, M0, and CCAP *)
      DELTA:= CSSM * CSSL - sqr(CSCPML);
      M0:= ((CSSM * CSPMLL) - (CSCPML * CSPMLM)) / DELTA;
      L0:= ((-CSCPML * CSPMLL) + (CSSL * CSPMLM)) / DELTA;
 
      (* Standard deviation of Ci *)
      SCi:=0; SCiS:=0; CCAP:=0;
      for I:= 1 to N do begin
        C[I]:= (M[I]-M0) * (L[I]-L0);
        SCi:= SCi + (C[I] - C[1]);
        SCiS:= SCiS + sqr(C[I] - C[1])
        end;
      CCAP:= SCi/N + C[1];
      SC:= sqrt((SCiS- sqr(SCi)/N) /(N-1));
 
      (* Calculate predicted fragment sizes and deviations from input sizes *)
      SDi:=0;SDiS:=0;
      for I:= 1 to N do begin
        PREDLEN[I]:= CCAP/(M[I]-M0) + L0;
        D[I]:= L[I] - PREDLEN[I];
        SDi:= SDi + D[I];
        SDiS:= SDiS + sqr(D[I])
        end;
 
      (* Standard deviation of fragment sizes *)
      SD:= sqrt((SDiS - sqr(SDi)/N) / (N-3))
 
    end; (* CALCULATE *)
 
  (******************************************************************)
  (* Print a report giving the standards and statistics.            *)
  (******************************************************************)
  procedure REPORT(var OUTFILE:text);
 
    var I : integer;
        TITLE:LINE;
    begin
      writeln('Type title to appear on output (<RETURN> for blank):');
      INPLINE(TITLE);
      WRITELINE(OUTFILE,TITLE,TITLE.LEN);
      if TITLE.LEN > 0 then writeln(OUTFILE);
 
      writeln(OUTFILE,'STD LEN':10,'DIST':10,'PRED LEN':10,'DEVIATION':10,
        '%DEV':10,'C[I]':10);
      for I:= 1 to N do
        writeln(OUTFILE,L[I]:10:2,M[I]:10:3,PREDLEN[I]:10:2,D[I]:10:3,
           D[I]/L[I]*100:10:3,C[I]:10:3);
      writeln(OUTFILE,'M0= ',M0:10:4,'  L0= ',L0:10:4,'  CCAP= ',CCAP:10:4);
      writeln(OUTFILE,'SC= ',SC:10:4,'  SD= ',SD:10:4);
      writeln(OUTFILE)
    end; (* REPORT *)
 
  (****************************************************************)
  (* Prompt user for migration distances of unknown fragments and *)
  (*  calculate the sizes of unknowns.                            *)
  (****************************************************************)
  procedure UNKNOWNS(var OUTFILE:text);
 
  var DISTANCE,LENGTH: real;
      NAME: LINE;
  begin
    writeln(OUTFILE,'UNKNOWN FRAGMENTS:');
    writeln(OUTFILE,'FRAGMENT':10,'DISTANCE':10,' PREDICTED LENGTH');
    writeln('Type an identifier (<=10 letters) for an unknown fragment');
    writeln('  (<RETURN> to quit)');
    INPLINE(NAME);
    while NAME.LEN > 0 do begin
      writeln('Distance migrated?');
      GETREAL(DISTANCE,0,MAXREAL);readln;
      LENGTH:= (CCAP/(DISTANCE-M0)) + L0;
      WRITELINE(OUTFILE,NAME,10);
      writeln(OUTFILE,DISTANCE:10:2,LENGTH:17:3);
      writeln('Type an identifier (<=10 letters) for an unknown fragment');
      writeln('  (<RETURN> to quit)');
      INPLINE(NAME)
      end;
    writeln(OUTFILE)
  end; (* UNKNOWNS *)
 
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

      (* Initialize horizontal output line *)
      with HLINE do begin
        for I:= 1 to MAXLINE do STR[I]:='_';
        LEN:=MAXLINE
        end;
     
      (* Open output file *)   
      writeln('Type output filename:');
      GETFILE(OUTFILE,'O',OFN);
      N:=0;
      
      (* MAIN LOOP *)
      repeat
        writeln;
        WRITELINE(output,HLINE,80);writeln;
        writeln('GEL','MAIN MENU':40);
        WRITELINE(output,HLINE,80);writeln;
        write('Output file: ':20);WRITELINE(output,OFN,OFN.LEN);writeln;
        WRITELINE(output,HLINE,80);writeln;
        writeln;
        writeln(' ':20,'1) Type in a set of standard fragments');
        writeln(' ':20,'2) Edit values of standard fragments');
        writeln(' ':20,'3) Add fragments');
        writeln(' ':20,'4) Delete fragments');
        writeln(' ':20,'5) Calculate sizes of unknowns (output to screen)');
        writeln(' ':20,'6) Calculate sizes of unknowns (output to file)');
        writeln(' ':20,'7) Open a new output file');
        WRITELINE(output,HLINE,80);writeln;
        writeln('Type number of your choice  (0 to quit)');
        readln(ANSWER);
        if ANSWER in ['1'..'7'] then
          case ANSWER of
            '1':READSTANDARDS;
            '2':EDIT;
            '3':ADDFRAGS;
            '4':DELETEFRAGS;
            '5':if N >= 4 then begin
                  CALCULATE(L0,M0,CCAP,SD,SC);
                  REPORT(output);
                  UNKNOWNS(output)
                  end
                else writeln('>>> There must be at least 4 standards. <<<');
            '6':if N >= 4 then begin
                  CALCULATE(L0,M0,CCAP,SD,SC);
                  REPORT(OUTFILE);
                  UNKNOWNS(OUTFILE)
                  end
                else writeln('>>> There must be at least 4 standards. <<<');
            '7':begin
(*!!!*)           CLOSE(OUTFILE); 
                  writeln('Type new output filename:');
                  GETFILE(OUTFILE,'O',OFN)
                end
            end; (* case *)
      until ANSWER='0';
(*!!!*)CLOSE(OUTFILE) 
    end. (* GEL *)
 
