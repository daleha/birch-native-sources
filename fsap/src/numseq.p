  (***********************************************************)
  (*                                                         *)
  (*  NUMSEQ    VERSION   6/26/2001  Standard Pascal         *)
  (*            Brian Fristensky                             *)
  (*            Dept. of Plant Science                       *)
  (*            University of Manitoba                       *)
  (*            Winnipeg, MB R3T 2N2  CANADA                 *)
  (*                                                         *)
  (*  Copyright (c) l982,1986,1988,1990  by Brian Fristensky *)
  (*  !!! in comment indicates feature which may need change *)
  (***********************************************************)

(* Revision history
06/26/2001 - Changed parameter menu to allow sequence lengths
             up to 11 digits for compatibility with GenBank.
11/05/2000 - In previous versions, if START and FINISH imply
             circularity in a linear molecule, NUMSEQ would
             just print through the ends. For example, if
             START > FINISH on the input strand, NUMSEQ
             would incorrectly print to the end, and then
             resume printing at 1. Now, NUMSEQ would correctly
             "fall off the end" of a linear molecule.
            
*)

  program NUMSEQ(input, output (*,INFILE,OUTFILE,GCFILE*));
(*!!!  Some Pascals require file parameters in program heading *)

  const MAXSEQ = 750000;
        MAXLINE = 150;
        MAXWORD = 25;
        VERSION = 'NUMSEQ            Version   6/26/2001';
 
  type NUCLEOTIDE  = (T,C,A,G,R,Y,M,W,S,K,D,H,V,B,N);
       SS          = set of NUCLEOTIDE;
       AMINOACID   = (ALA,CYS,ASP,GLU,PHE,GLY,HIS,ILE,LYS,LEU,MET,ASN,PRO,
                      GLN,ARG,SER,THR,VALINE,TRP,TYR,ASX,GLX,TERM,UNKX);
                   (* Note: VAL is a reserved word in some Pascals *)

       NA          = array[NUCLEOTIDE] of NUCLEOTIDE;
       SEQUENCE    = array[1..MAXSEQ] of NUCLEOTIDE;
 
       LETTERS     = packed array[1..3] of char;

       (* This table holds the amino acid assignments for the 64 codons.*)
       GENETICCODE = array[T..G,T..G,T..G] of AMINOACID;

       RULE = array[1..3] of NUCLEOTIDE;
      

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
 
  var (* File variables *)
      INFILE,                  (* input sequence file *)
      OUTFILE,                 (* output sequence file *)
      GCFILE          : text;  (* genetic code file *)
      IFN,OFN,GCFN    : LINE;  (* input, output, and genetic code filenames*)

      (* Sequence variables *)
      SEQ             : SEQUENCE;
      SEQLEN,CHOICE   : integer;
      SEQNAME         : WORD;  (* Req'd by READSEQ but not used *)
      INUC,COMP       : NA;
      NUCSET          : array[NUCLEOTIDE] of SS;
      NUCHAR          : array[NUCLEOTIDE] of char;
   
      (* Genetic code variables *)
      GC              : GENETICCODE;
      RULELIST        : array[1..64] of RULE;
      AALIST          : array[1..64] of AMINOACID;
      NUMRULES        : integer; 
      AAS             : array[AMINOACID] of LETTERS;

      (* GLOBAL PARAMETERS *)
      START,FINISH,STARTNO,
      GROUP,GPL,SENSE,FRAMES,STRANDS,I: integer;
      NUCCASE,
      COORD,WHICH,FORM,KIND,
      NUMBERS,NUCS,PEPTIDES: char;
      CIRCULAR       : boolean;
      HLINE          : LINE;
      ANSWER         : char; 

 
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
(* END MODULE AA         VERSION= 'SUNMODS     Version  6/26/01'; *)

 (*******************************************************)
 (* Initialize tables of nucleotides, amino acids etc.  *)
 (* This is done only once when the program begins.     *)
 (*******************************************************)
 procedure INITTABLES;
   begin
     (*   Initialize INUC and COMP, the arrays which *)
     (*   hold the  values of the NUCLEOTIDEs.      *)
        INUC[T]:=T;COMP[T]:=A;
        INUC[C]:=C;COMP[C]:=G;
        INUC[A]:=A;COMP[A]:=T;
        INUC[G]:=G;COMP[G]:=C;
        INUC[R]:=R;COMP[R]:=Y;
        INUC[D]:=D;COMP[D]:=H;
        INUC[V]:=V;COMP[V]:=B;
        INUC[M]:=M;COMP[M]:=K;
        INUC[K]:=K;COMP[K]:=M;
        INUC[B]:=B;COMP[B]:=V;
        INUC[H]:=H;COMP[H]:=D;
        INUC[Y]:=Y;COMP[Y]:=R;
        INUC[W]:=W;COMP[W]:=W;
        INUC[S]:=S;COMP[S]:=S;
        INUC[N]:=N;COMP[N]:=N;

     (* Initialize the nucleotide sets *)
        NUCSET[N]:=[T..N];
        NUCSET[R]:=[A,G,R];NUCSET[Y]:=[C,T,Y];NUCSET[M]:=[A,C,M];
        NUCSET[W]:=[A,T,W];NUCSET[S]:=[C,G,S];NUCSET[K]:=[G,T,K];
        NUCSET[D]:=[T,A,G,W,K,R];NUCSET[H]:=[T,C,A,Y,W,M];
        NUCSET[V]:=[C,A,G,M,S,R];NUCSET[B]:=[T,C,G,Y,K,S];
        NUCSET[T]:=[T];NUCSET[C]:=[C];NUCSET[A]:=[A];NUCSET[G]:=[G];
     
     (* Initialize GC using universal genetic code *)
        GC[T,T,T]:= PHE; GC[T,C,T]:=SER; GC[T,A,T]:=TYR; GC[T,G,T]:=CYS;
        GC[T,T,C]:= PHE; GC[T,C,C]:=SER; GC[T,A,C]:=TYR; GC[T,G,C]:=CYS;
        GC[T,T,A]:= LEU; GC[T,C,A]:=SER; GC[T,A,A]:=TERM;GC[T,G,A]:=TERM;
        GC[T,T,G]:= LEU; GC[T,C,G]:=SER; GC[T,A,G]:=TERM;GC[T,G,G]:=TRP;
        GC[C,T,T]:= LEU; GC[C,C,T]:=PRO; GC[C,A,T]:=HIS; GC[C,G,T]:=ARG;
        GC[C,T,C]:= LEU; GC[C,C,C]:=PRO; GC[C,A,C]:=HIS; GC[C,G,C]:=ARG;
        GC[C,T,A]:= LEU; GC[C,C,A]:=PRO; GC[C,A,A]:=GLN; GC[C,G,A]:=ARG;
        GC[C,T,G]:= LEU; GC[C,C,G]:=PRO; GC[C,A,G]:=GLN; GC[C,G,G]:=ARG;
        GC[A,T,T]:= ILE; GC[A,C,T]:=THR; GC[A,A,T]:=ASN; GC[A,G,T]:=SER;
        GC[A,T,C]:= ILE; GC[A,C,C]:=THR; GC[A,A,C]:=ASN; GC[A,G,C]:=SER;
        GC[A,T,A]:= ILE; GC[A,C,A]:=THR; GC[A,A,A]:=LYS; GC[A,G,A]:=ARG;
        GC[A,T,G]:= MET; GC[A,C,G]:=THR; GC[A,A,G]:=LYS; GC[A,G,G]:=ARG;
        GC[G,T,T]:= VALINE; GC[G,C,T]:=ALA; GC[G,A,T]:=ASP; GC[G,G,T]:=GLY;
        GC[G,T,C]:= VALINE; GC[G,C,C]:=ALA; GC[G,A,C]:=ASP; GC[G,G,C]:=GLY;
        GC[G,T,A]:= VALINE; GC[G,C,A]:=ALA; GC[G,A,A]:=GLU; GC[G,G,A]:=GLY;
        GC[G,T,G]:= VALINE; GC[G,C,G]:=ALA; GC[G,A,G]:=GLU; GC[G,G,G]:=GLY
   end; (* INITTABLES *)

  (* **************************************************** *)
  (* Initialize parameters for a new sequence.            *)
  (* **************************************************** *)
  procedure INITPARAM;
    begin
    (* Set default values for global parameters *)
      START:= 1;FINISH:= SEQLEN;STARTNO:=START;SENSE:=1;GROUP:=10;GPL:=7;
      NUCCASE:='U';FRAMES:=1;FORM:='L';COORD:='S';WHICH:='I';STRANDS:=1;
      KIND:='D';NUMBERS:='Y';NUCS:='Y';PEPTIDES:='N'
    end; (* INITPARAM *)
 
  (* ************************************************************ *)
  (* Initialize nucleotide and amino acid output strings          *)
  (* according to parameters.                                     *)
  (* ************************************************************ *)
    procedure INITSTRING;
         begin

           (* Initialize nucleotide output characters *)
           case NUCCASE of
             'U':begin
                 if KIND='D' then NUCHAR[T]:='T'
                 else NUCHAR[T]:='U';
                 NUCHAR[C]:='C';NUCHAR[A]:='A';NUCHAR[G]:='G';
                 NUCHAR[R]:='R';NUCHAR[D]:='D';NUCHAR[V]:='V';NUCHAR[M]:='M';
                 NUCHAR[K]:='K';NUCHAR[B]:='B';NUCHAR[H]:='H';NUCHAR[Y]:='Y';
                 NUCHAR[W]:='W';NUCHAR[S]:='S';NUCHAR[N]:='N'
                 end;
             'L':begin
                 if KIND='D' then NUCHAR[T]:='t'
                 else NUCHAR[T]:='u';
                 NUCHAR[C]:='c';NUCHAR[A]:='a';NUCHAR[G]:='g';
                 NUCHAR[R]:='r';NUCHAR[D]:='d';NUCHAR[V]:='v';NUCHAR[M]:='m';
                 NUCHAR[K]:='k';NUCHAR[B]:='b';NUCHAR[H]:='h';NUCHAR[Y]:='y';
                 NUCHAR[W]:='w';NUCHAR[S]:='s';NUCHAR[N]:='n'
                 end
             end; (* case NUCCASE *)

           (* Initialize amino acid output strings. *)
           case FORM of
             'S': begin
                    AAS[PHE]:='F  ';AAS[LEU]:='L  ';AAS[ILE]:='I  ';
                    AAS[MET]:='M  ';AAS[VALINE]:='V  ';AAS[SER]:='S  ';
                    AAS[PRO]:='P  ';AAS[THR]:='T  ';AAS[ALA]:='A  ';
                    AAS[TYR]:='Y  ';AAS[HIS]:='H  ';AAS[GLN]:='Q  ';
                    AAS[ASN]:='N  ';AAS[LYS]:='K  ';AAS[ASP]:='D  ';
                    AAS[GLU]:='E  ';AAS[CYS]:='C  ';AAS[TRP]:='W  ';
                    AAS[ARG]:='R  ';AAS[GLY]:='G  ';AAS[ASX]:='B  ';
                    AAS[GLX]:='Z  ';AAS[TERM]:='*  ';AAS[UNKX]:='X  '
                  end;
             'L': begin                                      
                    AAS[PHE]:='Phe';AAS[LEU]:='Leu';AAS[ILE]:='Ile';
                    AAS[MET]:='MET';AAS[VALINE]:='Val';AAS[SER]:='Ser';
                    AAS[PRO]:='Pro';AAS[THR]:='Thr';AAS[ALA]:='Ala';
                    AAS[TYR]:='Tyr';AAS[HIS]:='His';AAS[GLN]:='Gln';
                    AAS[ASN]:='Asn';AAS[LYS]:='Lys';AAS[ASP]:='Asp';
                    AAS[GLU]:='Glu';AAS[CYS]:='Cys';AAS[TRP]:='Trp';
                    AAS[ARG]:='Arg';AAS[GLY]:='Gly';AAS[ASX]:='Asx';
                    AAS[GLX]:='Glx';AAS[TERM]:='   ';AAS[UNKX]:='---'
                  end
              end (* case FORM *)
          end; (*INITSTRING*)

  (*************************************************************)
  (*  Read in an alternative genetic code from a datafile.     *)
  (*  The datafile takes the form of a two-dimensional table   *)
  (*  of amino acids as found in Watson, J.D. (1976), Molecular*)
  (*  Biology of the Gene pg.356, Table 13-7. The file         *)
  (*  should contain a single CAPITAL letter amino acid code   *)
  (*  for each position corresponding to a given codon.  This  *)
  (*  program will read in the first 64 characters that could  *)
  (*  be amino acid or stop ( * ) symbols. Blanks or lowercase *)
  (*  letters may be included for clarity.                     *)
  (*************************************************************)
  procedure READCODE(var F:text; var GCFN:LINE; var GC:GENETICCODE);
    var N1,N2,N3:NUCLEOTIDE;
        AACHAR:char;
    begin
      writeln('Type name of genetic code file:');
      GETFILE(F,'I',GCFN);
      for N1:= T to G do
        for N3:= T to G do     (* This is correct. Think about it. *)
          for N2:= T to G do begin
            AACHAR:= ' ';
            while (not eof(F)) and (not(AACHAR in ['A','C','D','E','F','G','H',
                 'I','K','L','M','N','P','Q','R','S','T','V','W','Y','*'])) do 
              read(F,AACHAR);
            GC[N1,N2,N3]:= AA(AACHAR)
            end;
(*!!!*) CLOSE(F) 
    end; (* READCODE *)
 
  (*************************************************************)
  (* Using the 64 unambiguous codons stored in GC, make the    *)
  (* the optimal list of ambiguity rules, assigning codons to  *)
  (* amino acids. All possible ambiguous codons are taken into *)
  (* account.                                                  *)
  (*************************************************************)
  procedure MAKERULES(GC:GENETICCODE);
   
    var N1,N2,N3 : NUCLEOTIDE; (* current nucleotides *)
        AMA      : AMINOACID; (* current amino acid *)
        RA,RB,
        FIRST,LASTOLDRULE    : integer; (* All are indices of rules *)
        S1,S2,S3:SS;

    (* Given a set, return the nucleotide symbol corresponding to the minimal
       set which contains all members of that set *)
    function NUCSYM(NS:SS):NUCLEOTIDE;
      var N1:NUCLEOTIDE;
        begin
          N1:= T;
          while not(NS <= NUCSET[N1]) do N1:= succ(N1);
          NUCSYM:= N1        
        end; (* NUCSYM *)

    (* Add a rule to the end of RULELIST *)
    procedure ADDRULE(N1,N2,N3:NUCLEOTIDE;AMA:AMINOACID);
      begin
        NUMRULES:= NUMRULES + 1;
        RULELIST[NUMRULES][1]:=N1;
        RULELIST[NUMRULES][2]:=N2;
        RULELIST[NUMRULES][3]:=N3;
        AALIST[NUMRULES]:= AMA
      end; (* ADDRULE *)

    (* Delete a rule from the RULELIST. Move the other rules down the list.*)
    procedure DELRULE(I:integer);
      var J:integer;
        begin
          J:= I + 1;
          while J <= NUMRULES do begin
            RULELIST[I]:= RULELIST[J]; AALIST[I]:= AALIST[J];
            I:= I + 1; J:= J + 1
            end;
          NUMRULES:= NUMRULES - 1
        end; (* DELRULE *)

    begin (* MAKERULES *)
      (* For each amino acid, make the minimal set of rules such that all
         possible codons for that amino acid are included in some rule *)
      NUMRULES:= 0;
      for AMA:= ALA to TERM do begin
        FIRST:= NUMRULES + 1;
       
        (* STEP 1: Add to RULELIST all codons for a given amino acid *)
        for N2:= T to G do        (* This is correct, too. *)
          for N1:= T to G do 
            for N3:= T to G do
              if (GC[N1,N2,N3]=AMA) then ADDRULE(N1,N2,N3,AMA);

        (* STEP2: Generalization of the rules.  If any two positions in rule A
           are both equal to the corresponding positions in rule B, then update
           the remaining position of rule A to include both members at that 
           position. Delete rule B. *)
        RA:= FIRST;
        while RA < NUMRULES do begin
          RB:= RA + 1;
          while RB <= NUMRULES do 
            if (RULELIST[RA][1] = RULELIST[RB][1]) and
               (RULELIST[RA][2] = RULELIST[RB][2]) then begin
               RULELIST[RA][3]:= NUCSYM(NUCSET[RULELIST[RA][3]] +
                                 NUCSET[RULELIST[RB][3]]);
               DELRULE(RB)
               end
            else if (RULELIST[RA][1] = RULELIST[RB][1]) and
                    (RULELIST[RA][3] = RULELIST[RB][3]) then begin
                    RULELIST[RA][2]:= NUCSYM(NUCSET[RULELIST[RA][2]] + 
                    NUCSET[RULELIST[RB][2]]);
                    DELRULE(RB)
                    end
            else if (RULELIST[RA][2] = RULELIST[RB][2]) and
                    (RULELIST[RA][3] = RULELIST[RB][3]) then begin
                    RULELIST[RA][1]:= NUCSYM(NUCSET[RULELIST[RA][1]] + 
                    NUCSET[RULELIST[RB][1]]);
                    DELRULE(RB)
                      end
            else RB:= RB + 1;
          RA:= RA + 1
          end; 

        (* STEP 3: Generalization of the rules.  If any two positions in rule A 
           both intersect with the corresponding positions in rule B, then
           create a new rule comprised of the intersection at those two
           positions and the sum of the remaining positions of rule A and
           rule B. *)
        RA:= FIRST;
        LASTOLDRULE:=NUMRULES; (* The last rule made in STEP 2 *)
        while RA < LASTOLDRULE do begin
          RB:= RA + 1;
          while RB <= LASTOLDRULE do begin
            S1:= NUCSET[RULELIST[RA][1]] * NUCSET[RULELIST[RB][1]];
            S2:= NUCSET[RULELIST[RA][2]] * NUCSET[RULELIST[RB][2]];
            S3:= NUCSET[RULELIST[RA][3]] * NUCSET[RULELIST[RB][3]];
            if (S1<>[]) and (S2<>[]) then ADDRULE(NUCSYM(S1),NUCSYM(S2),
                NUCSYM(NUCSET[RULELIST[RA][3]]+NUCSET[RULELIST[RB][3]]),
                AMA);
            if (S1<>[]) and (S3<>[]) then ADDRULE(NUCSYM(S1),
               NUCSYM(NUCSET[RULELIST[RA][2]]+NUCSET[RULELIST[RB][2]]),
                NUCSYM(S3),AMA);
            if (S2<>[]) and (S3<>[]) then ADDRULE(
                NUCSYM(NUCSET[RULELIST[RA][1]]+NUCSET[RULELIST[RB][1]]),
                NUCSYM(S2),NUCSYM(S3),AMA);
            RB:= RB + 1
            end;
          RA:= RA + 1
          end;

        (* STEP 4: Delete any rule which is entirely contained within 
           another rule *)
        RA:= FIRST;
        while RA < NUMRULES do begin
          RB:= RA + 1;
          while RB <= NUMRULES do 
            if (RULELIST[RA][1] in NUCSET[RULELIST[RB][1]]) and
               (RULELIST[RA][2] in NUCSET[RULELIST[RB][2]]) and
               (RULELIST[RA][3] in NUCSET[RULELIST[RB][3]]) 
               then  DELRULE(RA)
            else RB:= RB + 1;
          RA:= RA + 1
          end 

        end;  (* for AMA *)
     
      (* STEP 5: Add a rule for unknown codons *)
      NUMRULES:= NUMRULES + 1;
      RULELIST[NUMRULES][1]:=N;
      RULELIST[NUMRULES][2]:=N;
      RULELIST[NUMRULES][3]:=N;
      AALIST[NUMRULES]:= UNKX;

(*!!!*)(* Print rule list *)
(*    FORM:='L';INITSTRING;
      writeln('The following rules implement the genetic code:'); 
      for RA:= 1 to NUMRULES do begin
        for RB:= 1 to 3 do write(NUCHAR[RULELIST[RA][RB]]);
        write(' ':7);
        writeln(AAS[AALIST[RA]]:10)
        end;
      writeln('Press RETURN to continue:');
      readln *)
  end; (* MAKERULES *)

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
          writeln('Type new value for ',PNAME,'  (CURRENT VALUE: ',P,')');
          GETINTEGER(P,LOW,HIGH);
          readln
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
        true:write('Topology:   CIRCULAR');
        false:write('Topology:     LINEAR')
        end;      
      writeln('Length: ':14,SEQLEN:11,' nt');
      WRITELINE(output,HLINE,80);writeln;
      writeln(' ':12,'Parameter   Description/Response                 Value');
      WRITELINE(output,HLINE,80);writeln;
      writeln(' ':12,' 1)START    first nucleotide printed',START:16);
      writeln(' ':12,' 2)FINISH   last  nucleotide printed',FINISH:16);
      writeln(' ':12,' 3)NUCCASE   U:(A,G,C,T...), l:(a,g,c,t...)',NUCCASE:9);
      writeln(' ':12,' 4)STARTNO  number of starting nucleotide',STARTNO:11);
      writeln(' ':12,' 5)GROUP    number every GROUP nucleotides',GROUP:10);
      writeln(' ':12,' 6)GPL      number of GROUPs printed per line',GPL:7);
      writeln(' ':12,' 7)WHICH    I: input strand  O: opposite strand',WHICH:5);
      writeln(' ':12,' 8)STRANDS  1: one  strand,  2:both strands',STRANDS:9);
      writeln(' ':12,' 9)KIND     R:RNA            D:DNA',KIND:18);
      writeln(' ':12,'10)NUMBERS  Number  the sequence    (Y or N)',NUMBERS:8);
      writeln(' ':12,'11)NUCS     Print nucleotide seq.   (Y or N)',NUCS:8);
      writeln(' ':12,'12)PEPTIDES Print amino acid seq.   (Y or N)',PEPTIDES:8);
      writeln(' ':12,'13)FRAMES   1 for this frame, 3 for 3 frames',FRAMES:8);
      writeln(' ':12,'14)FORM     L:3 letter amino acid, S: 1 letter',FORM:6);
      WRITELINE(output,HLINE,80);writeln;
      end; (* DISPLAY *)
 
      begin
        (* Prompt user for new parameter values *)
        repeat
          page(output);
          DISPLAY;
          if WHICH='O' then begin
            writeln('Be sure START and FINISH values are appropriate for',
            ' WHICH=O')
            end;
          writeln('Type number of parameter you wish to change ',
            '(0 to continue)');
          GETINTEGER(RESPONSE,0,14);
          readln;
          if RESPONSE in [1..14] then
            case RESPONSE of
               1:begin
                   GETNUMBER(START,'START     ',1,SEQLEN);
                   STARTNO:=START
                   end;
               2:GETNUMBER(FINISH,'FINISH    ',1,SEQLEN);
               3:GETCHAR(NUCCASE,'NUCCASE   ',['U','L']);
               4:GETNUMBER(STARTNO,'STARTNO   ',-MAXSEQ,MAXSEQ);
               5:GETNUMBER(GROUP,'GROUP     ',3,160);
               6:GETNUMBER(GPL,'GPL       ',1,160 div GROUP);
               7:GETCHAR(WHICH,'WHICH     ',['I','O']);
               8:GETNUMBER(STRANDS,'STRANDS   ',1,2);
               9:GETCHAR(KIND,'KIND      ',['R','D']);
              10:GETCHAR(NUMBERS,'NUMBERS   ',['Y','N']);
              11:GETCHAR(NUCS,'NUCS      ',['Y','N']);
              12:GETCHAR(PEPTIDES,'PEPTIDES  ',['Y','N']);
              13:GETNUMBER(FRAMES,'FRAMES    ',1,3);
              14:GETCHAR(FORM,'FORM      ',['S','L'])
              end;
          if (PEPTIDES='Y') and (GROUP mod 3 > 0) then begin
            writeln('>>>> GROUP must be divisible by 3! Press <RETURN>');
            readln;
            RESPONSE:= -1
            end
        until RESPONSE= 0;
        if STARTNO=START then COORD:='S' else COORD:='U'
     end; (* PARAMETERS *)
 

  (*************************************************************)
  (*  Write the sequence to OUTFILE in the specified format.   *)
  (*************************************************************)
  procedure PRINTSEQ(var F:text);
    var  I,POS,POSITION,GROUPSLEFT,NUCPERLINE,GROUPWIDTH,DELTA,
         NUCLEFT,SKIP,NUCSKIP,FRAMESHIFT,THISLINE,NUMBER,
         NUMCODONS,S2,S3:integer;
         LASTLINE:boolean;
         NAME:LINE;
 
    (* Skip the indicated number of spaces *)
    procedure ADVANCE(NUM:integer);
      var I:integer;
      begin
        for I:= 1 to NUM do write(F,' ')
      end; (* ADVANCE *)
 
    (* Write a line of numbers *)
    procedure WRITENUM;
      var WIDTH,NEXTNUM:integer;
      begin
        WIDTH:= 0;
        while WIDTH < THISLINE do begin
          (* 0 position not permitted in output. *)
          NEXTNUM:= NUMBER + DELTA;
          case COORD of
            'S': if (WHICH= 'I') and (NEXTNUM > SEQLEN) then
                          NUMBER:= DELTA - (SEQLEN - NUMBER)
                 else if (WHICH= 'O') and (NEXTNUM < 1) then
                          NUMBER:= SEQLEN + NEXTNUM
                 else NUMBER:= NEXTNUM;
            'U': if NEXTNUM= 0 then NUMBER:= SENSE
                 else if NUMBER/NEXTNUM < 0 then NUMBER:= NEXTNUM + SENSE
                 else NUMBER:= NEXTNUM
            end;
          write(F,NUMBER:GROUPWIDTH);
          ADVANCE(SKIP);
          WIDTH:= WIDTH+1
          end;
        writeln(F)
      end;   (* WRITENUM *)
 
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
 
    (* Write a line of bases *)
    procedure WRITENUC(NUC:NA; var POS:integer);
      var WIDTH,I : integer;
      begin
        I:= 0; WIDTH:= 0;
        while WIDTH < THISLINE do begin
          POS:= NEXT(POS);
          if POS>0 then write(F,NUCHAR[NUC[SEQ[POS]]]);
          WIDTH:= WIDTH + 1;
          I:= I+1; if I = NUCSKIP then begin ADVANCE(SKIP);I:=0 end
          end;
        writeln(F)
      end; (* WRITENUC *)
 
    (* Write a line of amino acids *)
    procedure WRITEAA(NUC:NA;POS,NUMCODONS:integer);
      var WIDTH,P2,P3 : integer;
      (* Given a codon, return an amino acid *)
      function AAMATCH(N1,N2,N3:NUCLEOTIDE):AMINOACID;
        var I:integer;    (* rule currently pointed to *)
            NOTFOUND:boolean;
        begin
          I:=1;
          NOTFOUND:= true;
          while NOTFOUND do
            if N1 in NUCSET[RULELIST[I][1]] then
              if N2 in NUCSET[RULELIST[I][2]] then
                if N3 in NUCSET[RULELIST[I][3]] then NOTFOUND:= false
                else I:= I + 1
              else I:= I + 1
            else I:= I + 1;
          AAMATCH:= AALIST[I]
        end; (* AAMATCH *)

      begin
        WIDTH:= 0;
        POS:= NEXT(POS)+(FRAMESHIFT*SENSE);
        ADVANCE(FRAMESHIFT);
        while (WIDTH < NUMCODONS) do begin
          P2:= NEXT(POS);P3:= NEXT(P2);
          if (P2 > 0) and (P3 > 0) then
            write(F,AAS[AAMATCH(NUC[SEQ[POS]],NUC[SEQ[P2]],NUC[SEQ[P3]])]);
          ADVANCE(SKIP);
          WIDTH:= WIDTH + 1;
          POS:= NEXT(P3)
          end;
        writeln(F)
      end; (* WRITEAA *)
 
    begin
      writeln('Type title to appear on output (<RETURN> for blank):');
      INPLINE(NAME);
      WRITELINE(F,NAME,NAME.LEN);
      if NAME.LEN > 0 then writeln(F);
      (* Set parameters dependant on the strand *)
      case WHICH of
        'I': begin
               SENSE:= 1;
               DELTA:= GROUP;
               if START <= FINISH then NUCLEFT:= FINISH - START + 1
               else if CIRCULAR then NUCLEFT:= SEQLEN - START + FINISH + 1
                    else NUCLEFT:= SEQLEN - START + 1
               end;
        'O': begin
               SENSE:= -1;
               if  COORD= 'S' then DELTA:= -GROUP
               else DELTA:= GROUP;
               if START >= FINISH then NUCLEFT:= START - FINISH + 1
               else if CIRCULAR then NUCLEFT:=START + (SEQLEN - FINISH + 1)
                    else NUCLEFT:= START
               end
         end;
      NUCPERLINE:= GROUP * GPL;
      (* Set parameters dependant upon amino acid seq. *)
      if PEPTIDES= 'Y' then begin
        if FRAMES=3 then begin
          SKIP:=0; NUCSKIP:= NUCPERLINE+1 end
        else begin
          SKIP:=1; NUCSKIP:= 3 end;
        GROUPWIDTH:= GROUP + (SKIP *(GROUP div 3))-SKIP;
        S2:= SENSE*2;S3:= SENSE*3;
        end
      else begin
        SKIP:= 1;
        NUCSKIP:= GROUP;
        GROUPWIDTH:= GROUP
        end;
 
      POSITION:= START- SENSE;
      if COORD='S' then NUMBER:= START-SENSE
      else NUMBER:= STARTNO-1;
      GROUPSLEFT:= NUCLEFT div GROUP;
      LASTLINE:= false;
 
      repeat
        (* Write numbers *)
        if GROUPSLEFT < GPL then THISLINE:= GROUPSLEFT
                            else THISLINE:= GPL;
        if NUMBERS ='Y' then WRITENUM;
        GROUPSLEFT:= GROUPSLEFT - THISLINE;
 
        (* Write bases *)
        if NUCLEFT < NUCPERLINE then THISLINE:= NUCLEFT
                                else THISLINE:= NUCPERLINE;
        POS:= POSITION;
        if NUCS='Y' then begin
          case WHICH of
            'I':WRITENUC(INUC,POS);
            'O':WRITENUC(COMP,POS)
            end;
          if STRANDS=2 then begin
            POS:= POSITION;
            case WHICH of
              'I':WRITENUC(COMP,POS);
              'O':WRITENUC(INUC,POS)
              end
            end
          end
        else for I:= 1 to THISLINE do POS:= NEXT(POS);
        NUCLEFT:= NUCLEFT- THISLINE;
 
        (* Write amino acids *)
        if PEPTIDES='Y' then begin
          FRAMESHIFT:=0;
          while FRAMESHIFT < FRAMES do begin
            NUMCODONS:= THISLINE div 3;
            case WHICH of
              'I':WRITEAA(INUC,POSITION,NUMCODONS);
              'O':WRITEAA(COMP,POSITION,NUMCODONS)
              end;
            FRAMESHIFT:= FRAMESHIFT + 1
            end
          end;
        POSITION:= POS;
 
        (* Write a blank line *)
        if NUMBERS='Y' then writeln(F);
        if NUCLEFT= 0 then LASTLINE:= true
      until LASTLINE
      end; (* PRINTSEQ *)
 
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
      
      (* Set up initial values for nucleotide tables and genetic code *)
      INITTABLES;
      MAKERULES(GC);

      (* Initialize horizontal output line *)
      with HLINE do begin
        for I:= 1 to MAXLINE do STR[I]:='_';
        LEN:=MAXLINE
        end;
      
      (* Read in initial sequence and set up parameters. *)
      writeln('Enter sequence filename:');
      GETFILE(INFILE,'I',IFN);
      READSEQ(INFILE,SEQ,SEQLEN,SEQNAME,CIRCULAR);
      INITPARAM;

      (* Open initial output file *)
      writeln('Type output filename:');
      GETFILE(OUTFILE,'O',OFN);

      (* MAIN MENU *)
      GCFN.LEN:=0;
      repeat
        writeln;
        WRITELINE(output,HLINE,80);writeln;
        writeln('NUMSEQ','MAIN MENU':30);
        WRITELINE(output,HLINE,80);writeln;
        write('Input file:        ');WRITELINE(output,IFN,IFN.LEN);writeln;
        write('Output file:       ');WRITELINE(output,OFN,OFN.LEN);writeln;
        write('Genetic code file: ');WRITELINE(output,GCFN,GCFN.LEN);writeln;
        WRITELINE(output,HLINE,80);writeln;
        writeln(' ':20,'1) Read in a new sequence');
        writeln(' ':20,'2) Open a new output file');
        writeln(' ':20,'3) Read in an alternative genetic code');
        writeln(' ':20,'4) Change parameters');
        writeln(' ':20,'5) Write output to screen');
        writeln(' ':20,'6) Write output to file');
        WRITELINE(output,HLINE,80);writeln;
        writeln('Type the number of your choice  (0 to quit program)');
        GETINTEGER(CHOICE,0,6);readln;

        case CHOICE of
          0:;
          1:begin
            writeln('Enter sequence filename:');
            GETFILE(INFILE,'I',IFN);
            READSEQ(INFILE,SEQ,SEQLEN,SEQNAME,CIRCULAR);
            START:=1; FINISH:=SEQLEN; STARTNO:=START
            end;
          2:begin
(*!!!*)       CLOSE(OUTFILE); 
              writeln('Type output filename:');
              GETFILE(OUTFILE,'O',OFN);
            end;
          3:begin
              READCODE(GCFILE,GCFN,GC);
              MAKERULES(GC)
            end;
          4:PARAMETERS;
          5:begin
              INITSTRING;
              PRINTSEQ(output);
              write('Press RETURN to continue');
              readln
            end;
          6:begin
              INITSTRING;
              PRINTSEQ(OUTFILE)
            end
          end (* case *)
        until CHOICE = 0;
(*!!!*)  CLOSE(OUTFILE) 
    end. (* NUMSEQ *)
