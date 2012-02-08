  (* ********************************************************  *)
  (*                                                           *)
  (*   PROT2NUC  Version  8/10/94, Standard Pascal             *)
  (*             Brian Fristensky                              *)
  (*             Dept. of Plant Science                        *)
  (*             University of Manitoba                        *)
  (*             Winnipeg, MB R3T 2N2 CANADA                   *)
  (*                                                           *)
  (* SYNOPSIS                                                  *)
  (*    prot2nuc [-l<n> -g<n>]                                 *)
  (*                                                           *)
  (* DESCRIPTION                                               *)
  (*    Reverse translates protein to nucleic acid sequence.   *)
  (*    File is Pearson (.wrp) format file. see Lipman et al   *)
  (*                                                           *)
  (*    -l<n>  integer, length of output line in codons        *)
  (*           (default = 20)                                  *)
  (*                                                           *)
  (*    -g<n>  integer; number ever n codons (default = 5)     *)
  (*                                                           *)
  (* Copyright (c) 1994 by Brian Fristensky.                   *)
  (* !!! in comment indicates feature which may need change.   *)
  (*  *******************************************************  *)
 
  program  PROT2NUC(input,output);
 (*!!! Some Pascals require file parameters in program heading *)
 
  const VERSION = 'PROT2NUC       Version  8/10/94';
        MAXSEQ  = 9000;
        MAXWORD = 25;
        MAXLINE = 70; (* only used for ARGUMENT *)
(* BEGIN MODULE STARTARGNUM *)
	STARTARGNUM=1;    (* SUN Pascal: ARG(1) is 1st command line argument*)
      (*STARTARGNUM=2;       ATT Pascal: ARG(2) is 1st command line argument*)
(* END MODULE STARTARGNUM         VERSION= 'SUNMODS     Version  8/ 9/94'; *)

 
  type  AMINOACID = (GLY,ALA,VALINE,LEU,ILE,MET,PHE,PRO,SER,THR,CYS,ASN,
                    GLN,TYR,TRP,ASP,GLU,HIS,LYS,ARG,ASX,GLX,TERM,UNKX);
        PROTEIN   = array[1..MAXSEQ] of AMINOACID;
        TRIPLET   = array[1..3] of char;
        DEGENERATE = array[AMINOACID] of TRIPLET;

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

   var 
      (* Command line argument variables. *)
      ARGUMENT     :CHARARRAY; (*     command line argument *)
      ARGNUM,                  (* # of     "     "      "   *)
      POS          : integer;  (* posn. in "     "      "   *)
      UNREADOPTIONS:boolean;

      (* Amino acid sequence variables *)
      NAME         : WORD;
      PR           : PROTEIN;
      SEQLEN       : integer;

      (* Variables for printing *)
      CODONSPERLINE,
      GROUP        : integer;
      AACHAR       : array[AMINOACID] of char;
      CODON1,
      CODON2       : DEGENERATE;
      
    (* **************************************************** *)
    (*              INITIALIZATION  PROCEDURES              *)
    (* **************************************************** *)
    procedure INITIALIZE;

      const C2STR = '   ';

      var AA1:AMINOACID;

      begin
      (*   Initialize AACHAR array which holds the            *)
      (*   character values of the AMINOACIDs.                *)
      AACHAR[GLY]:='G'; AACHAR[ALA]:='A';AACHAR[VALINE]:='V'; AACHAR[LEU]:='L';
      AACHAR[ILE]:='I'; AACHAR[MET]:='M';AACHAR[PHE]:='F'; AACHAR[PRO]:='P';
      AACHAR[SER]:='S'; AACHAR[THR]:='T';AACHAR[CYS]:='C'; AACHAR[ASN]:='N';
      AACHAR[GLN]:='Q'; AACHAR[TYR]:='Y';AACHAR[TRP]:='W'; AACHAR[ASP]:='D';
      AACHAR[GLU]:='E'; AACHAR[HIS]:='H';AACHAR[LYS]:='K'; AACHAR[ARG]:='R';
      AACHAR[ASX]:='B'; AACHAR[GLX]:='Z';AACHAR[TERM]:='*';
      AACHAR[UNKX]:='X';

      (* Initialize CODON1 & CODON2. All aminoacids have an equivalent
         TRIPLET in the CODON1 array. For amino acids that require two
         degenerate codons (eg. Leu = TTR, CTN), a triplet is needed
         for CODON2. For amino acids that can be represented in a single
         degenerate codon (eg. Cys = TGY), the TRIPLET in CODON2 is
         the value of the constant C2STR. *)
      CODON1[PHE]:='TTy'; CODON1[LEU]:='CTn'; CODON1[SER]:='TCn';
      CODON1[TYR]:='TAy'; CODON1[TERM]:='TAr'; CODON1[CYS]:='TGy';
      CODON1[TRP]:='TGG'; CODON1[PRO]:='CCn'; CODON1[HIS]:='CAy';
      CODON1[GLN]:='CAr'; CODON1[ARG]:='CGn'; CODON1[ILE]:='ATh';
      CODON1[MET]:='ATG'; CODON1[THR]:='ACn'; CODON1[ASN]:='AAy';
      CODON1[LYS]:='AAr'; CODON1[VALINE]:='GTn'; CODON1[ALA]:='GCn'; 
      CODON1[ASP]:='GAy'; CODON1[GLU]:='GAr'; CODON1[GLY]:='GGn';
      CODON1[ASX]:='rAy'; CODON1[GLX]:='sAr'; CODON1[UNKX]:='nnn';

      for AA1:= GLY to UNKX do CODON2[AA1]:=C2STR;
      CODON2[LEU]:='TTr'; CODON2[TERM]:='TGA'; CODON2[ARG]:='AGr';
      CODON2[SER]:='AGy'
      end; (* INITIALIZE *)

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
(* END MODULE AA         VERSION= 'SUNMODS     Version  8/ 9/94'; *)

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

    (*  ******************************************* *)
    (*  Read a protein    sequence from SEQFILE     *)
    (*   and store it in PR.                        *)
    (*  ******************************************* *)
    procedure READPRO(var SEQFILE:text; var PR:PROTEIN; var SEQLEN:integer;
                      var NAME:WORD);
    var CH:char;
        TOOLONG:boolean;
    begin
      SEQLEN := 0; TOOLONG:=false;
      (* Read in the sequence name, if present. Name may be preceeded by
         comment lines. *)
      while SEQFILE^=';' do readln(SEQFILE);
      if SEQFILE^='>' then begin
         read(SEQFILE,CH); READWORD(SEQFILE,NAME);readln(SEQFILE) end;
      
      (* Read in the sequence *)
      while (not eof(SEQFILE)) and (not (SEQFILE^='>')) do begin
        if SEQFILE^=';' then (* comment line, ignore *)
        else while not eoln(SEQFILE) do begin
          read(SEQFILE,CH);
          if CH in ['G','A','V','L','I','M','F','P','S','T','C','N','Q',
                      'Y','W','D','E','H','K','R','B','Z','*','X'] then
            if SEQLEN < MAXSEQ-2 then begin
              SEQLEN := SEQLEN + 1;
              PR[SEQLEN]:=AA(CH);
              end
            else TOOLONG:= true
          end; (* eoln *)
          readln(SEQFILE)
          end; (* eof *)
      if TOOLONG then begin
         writeln(
        '>>> WARNING! Sequence length exceeds MAXSEQ-2. Seq. truncated.')
        end (* TOOLONG *)
      end; (* READPRO *)

    (*  ******************************************* *)
    (*  Print a header listing the amino acid       *)
    (*  and nucleic acid one-letter codes.          *)
    (*  ******************************************* *)
    procedure HEADER;
    begin
      writeln(VERSION);
      writeln;
      writeln('     IUPAC-IUP AMINO ACID SYMBOLS');
      writeln('     [J. Biol. Chem. 243, 3557-3559 (1968)]');
      writeln;
      writeln('          Phe         F          Leu         L          Ile         I');
      writeln('          Met         M          Val         V          Ser         S');
      writeln('          Pro         P          Thr         T          Ala         A');
      writeln('          Tyr         Y          His         H          Gln         Q');
      writeln('          Asn         N          Lys         K          Asp         D');
      writeln('          Glu         E          Cys         C          Trp         W');
      writeln('          Arg         R          Gly         G          STOP        *');
      writeln('          Asx         B          Glx         Z          UNKNOWN     X');
      writeln;
      writeln;
      writeln('     IUPAC-IUB SYMBOLS FOR NUCLEOTIDE NOMENCLATURE');
      writeln('     [Cornish-Bowden (1985) Nucl. Acids Res. 13: 3021-3030.]');
      writeln;
      writeln('     Symbol         Meaning              | Symbol         Meaning');
      writeln('     ------------------------------------+---------------------------------');
      writeln('     G              Guanine              | k              G or T');
      writeln('     A              Adenine              | s              G or C');
      writeln('     C              Cytosine             | w              A or T');
      writeln('     T              Thymine              | h              A or C or T');
      writeln('     U              Uracil               | b              G or T or C');
      writeln('     r              Purine (A or G)      | v              G or C or A');
      writeln('     y              Pyrimidine (C or T)  | d              G or T or A');
      writeln('     m              A or C               | n              G or A or T or C');
      writeln
    end; (* HEADER *)

  (*************************************************************)
  (*  Write the sequence to OUTFILE in the specified format.   *)
  (*************************************************************)
  procedure PRINTSEQ(var F:text; PR:PROTEIN; SEQLEN:integer; NAME:WORD);
    var  POS,POSITION,GROUPWIDTH,
         CODONSLEFT,THISLINE,NUMBER:integer;
         LASTLINE:boolean;
 
    (* Write a line of numbers *)
    procedure WRITENUM(THISLINE:integer; var NUMBER:integer);
      var WIDTH,NEXTNUM:integer;
      begin
        WIDTH:= 0;
        while WIDTH < THISLINE do begin
          NEXTNUM:= NUMBER + GROUP;
          if (NEXTNUM <= SEQLEN) then begin
             NUMBER:= NEXTNUM;
             write(F,NUMBER:GROUPWIDTH)
             end;
          WIDTH:= WIDTH+GROUP
          end;
        writeln(F)
      end;   (* WRITENUM *)
 
    (* Write a line of amino acids *)
    procedure WRITEAA(THISLINE:integer; var POS:integer);
      var WIDTH : integer;
      begin
        WIDTH:= 0;
        while (WIDTH < THISLINE) do begin
          POS:= POS+1;
          write(F,AACHAR[PR[POS]],'  ');
          WIDTH:= WIDTH + 1;
          end;
        writeln(F)
      end; (* WRITEAA *)
 
    (* Write a line of degenerate triplets *)
    procedure WRITENUC(THISLINE:integer; NUCSEQ:DEGENERATE; var POS:integer);
      var WIDTH: integer;
      begin
        WIDTH:= 0;
        while WIDTH < THISLINE do begin
          POS:= POS +1;
          write(F,NUCSEQ[PR[POS]]);
          WIDTH:= WIDTH + 1
          end;
        writeln(F)
      end; (* WRITENUC *)

     begin
      WRITEWORD(F,NAME,NAME.LEN);
      if NAME.LEN > 0 then writeln(F);
      CODONSLEFT:=SEQLEN;
      POSITION:= 0;
      NUMBER:= 0;
      GROUPWIDTH:= GROUP * 3; (* width in characters, rather than codons *)
      LASTLINE:= false;
 
      repeat
        if CODONSLEFT < CODONSPERLINE then THISLINE:= CODONSLEFT
                                      else THISLINE:= CODONSPERLINE;

        (* Write numbers *)
        WRITENUM(THISLINE,NUMBER);

        (* Write amino acids *)
        POS:= POSITION;
        WRITEAA(THISLINE,POS);

        (* Write bases *)
        POS:= POSITION;
        WRITENUC(THISLINE,CODON1,POS);
        POS:= POSITION;
        WRITENUC(THISLINE,CODON2,POS);
        POSITION:=POS;
 
        (* Write a blank line *)
        writeln(F);
        CODONSLEFT:= CODONSLEFT - THISLINE;
        if CODONSLEFT= 0 then LASTLINE:= true
      until LASTLINE
      end; (* PRINTSEQ *)
   
  (* ----------------------------------------------------------  *)
  (* ----------------- MAIN  PROCEDURE  -----------------------  *)
    begin
      INITIALIZE;
      CODONSPERLINE:=25;
      GROUP:=5;

      (* Read options from the command line. *)
      ARGNUM:= STARTARGNUM;
      UNREADOPTIONS:=true;
      while UNREADOPTIONS do
        if ARGNUM < argc then begin
          argv(ARGNUM,ARGUMENT);
          if ARGUMENT[1]='-' then begin
            if ARGUMENT[2] in ['l','g'] then begin
              POS:=3;
              case ARGUMENT[2] of
                'l':CODONSPERLINE:=NUMBER(ARGUMENT,POS);
                'g':GROUP:=NUMBER(ARGUMENT,POS)
                end (* case*)
              end (* l,g *)
            end (* '-' *)
          else UNREADOPTIONS:=false;
          ARGNUM:=ARGNUM+1
          end (* if ARGNUM <= argc *)
        else UNREADOPTIONS:=false;

      (* GROUP must divide evenly into CODONSPERLINE *)
      if not (CODONSPERLINE mod GROUP = 0) then begin
             GROUP:=5;
             CODONSPERLINE:=25
             end;

      (* Read the protein and print the reverse translation. *)
      READPRO(input,PR,SEQLEN,NAME);
      HEADER;
      PRINTSEQ(output,PR,SEQLEN,NAME)
    end. (* PROT2NUC *)

