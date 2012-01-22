   (***********************************************************)
  (*                                                         *)
  (*  RIBOSOME  VERSION   8/12/93  Standard Pascal           *)
  (*            Brian Fristensky                             *)
  (*            Dept. of Plant Science                       *)
  (*            University of Manitoba                       *)
  (*            Winnipeg, MB R3T 2N2   CANADA                *)
  (*                                                         *)
  (*  Copyright (c) l989, 1990 by Brian Fristensky           *)
  (*  !!! in comment indicates feature which may need change *)
  (***********************************************************)
  program RIBOSOME(input, output (*,GCFILE*));
(*!!!  Some Pascals require file parameters in program heading *)

  const MAXSEQ =  32700;
        MAXWORD = 25;
        STARTARGNUM=1;
 
  type NUCLEOTIDE  = (T,C,A,G,R,Y,M,W,S,K,D,H,V,B,N);
       SS          = set of NUCLEOTIDE;
       AMINOACID   = (ALA,CYS,ASP,GLU,PHE,GLY,HIS,ILE,LYS,LEU,MET,ASN,PRO,
                      GLN,ARG,SER,THR,VALINE,TRP,TYR,ASX,GLX,TERM,UNKX);
                   (* Note: VAL is a reserved word in some Pascals *)

       SEQUENCE    = array[1..MAXSEQ] of NUCLEOTIDE;
 
       (* This table holds the amino acid assignments for the 64 codons.*)
       GENETICCODE = array[T..G,T..G,T..G] of AMINOACID;

       RULE = array[1..3] of NUCLEOTIDE;
      

(* BEGIN MODULE TYPE.WORD *)
       (*   <word>::= <non-blank char>[<non-blank char>] *)
       WORD    = record
                 LEN:integer;
                 STR:array[1..MAXWORD] of char
                 end;
(* END MODULE TYPE.WORD         VERSION= 'SUNMODS     Version  9/12/91'; *)
 

  var (* File variables *)
      GCFILE          : text;  (* genetic code file *)
      ARGNUM          : integer;
      ARGUMENT        : packed array[1..132] of char;

      (* Sequence variables *)
      SEQ             : SEQUENCE;
      SEQLEN          : integer;
      SEQNAME         : WORD;
      NUCSET          : array[NUCLEOTIDE] of SS;
      NUCHAR          : array[NUCLEOTIDE] of char;
   
      (* Genetic code variables *)
      GC              : GENETICCODE;
      RULELIST        : array[1..64] of RULE;
      AALIST          : array[1..64] of AMINOACID;
      NUMRULES        : integer; 
      AAS             : array[AMINOACID] of char;

      (* GLOBAL PARAMETERS *)
      START,FINISH    : integer;

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
(* END MODULE FILEARGS         VERSION= 'SUNMODS     Version  9/12/91'; *)
   
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
(* END MODULE NUC         VERSION= 'SUNMODS     Version  9/12/91'; *)

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
(* END MODULE AA         VERSION= 'SUNMODS     Version  9/12/91'; *)

 (*******************************************************)
 (* Initialize tables of nucleotides, amino acids etc.  *)
 (* This is done only once when the program begins.     *)
 (*******************************************************)
 procedure INITTABLES;
   begin
     (*   Initialize NUCHAR the array which *)
     (*   hold the  values of the NUCLEOTIDEs.      *)
     NUCHAR[T]:='T'; NUCHAR[C]:='C'; NUCHAR[A]:='A'; NUCHAR[G]:='G';
     NUCHAR[R]:='R'; NUCHAR[D]:='D'; NUCHAR[V]:='V'; NUCHAR[M]:='M';
     NUCHAR[K]:='K'; NUCHAR[B]:='B'; NUCHAR[H]:='H'; NUCHAR[Y]:='Y';
     NUCHAR[W]:='W'; NUCHAR[S]:='S'; NUCHAR[N]:='N';

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

  (* ************************************************************ *)
  (* Initialize amino acid output strings according to parameters.*)
  (* ************************************************************ *)
    procedure INITSTRING;
         begin
           AAS[PHE]:='F';AAS[LEU]:='L';AAS[ILE]:='I';
           AAS[MET]:='M';AAS[VALINE]:='V';AAS[SER]:='S';
           AAS[PRO]:='P';AAS[THR]:='T';AAS[ALA]:='A';
           AAS[TYR]:='Y';AAS[HIS]:='H';AAS[GLN]:='Q';
           AAS[ASN]:='N';AAS[LYS]:='K';AAS[ASP]:='D';
           AAS[GLU]:='E';AAS[CYS]:='C';AAS[TRP]:='W';
           AAS[ARG]:='R';AAS[GLY]:='G';AAS[ASX]:='B';
           AAS[GLX]:='Z';AAS[TERM]:='*';AAS[UNKX]:='X'
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
  procedure READCODE(var F:text; var GC:GENETICCODE);
    var N1,N2,N3:NUCLEOTIDE;
        AACHAR:char;
    begin
      for N1:= T to G do
        for N3:= T to G do     (* This is correct. Think about it. *)
          for N2:= T to G do begin
            AACHAR:= ' ';
            while (not eof(F)) and (not(AACHAR in ['A','C','D','E','F','G','H',
                 'I','K','L','M','N','P','Q','R','S','T','V','W','Y','*'])) do 
              read(F,AACHAR);
            GC[N1,N2,N3]:= AA(AACHAR)
            end;
(*!!! CLOSE(F) *)
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
(*    INITSTRING;
      writeln('The following rules implement the genetic code:'); 
      for RA:= 1 to NUMRULES do begin
        for RB:= 1 to 3 do write(NUCHAR[RULELIST[RA][RB]]);
        write(' ':7);
        writeln(AAS[AALIST[RA]]:10)
        end
      writeln;
      readln *)
  end; (* MAKERULES *)

  (*************************************************************)
  (*  Read a nucleic acid sequence from SEQFILE.               *)
  (*************************************************************)
  procedure READSEQUENCE(var SEQFILE:text; var NAME:WORD; var SEQ:SEQUENCE;
                         var START,FINISH:integer);
    var CH:char;
        TOOLONG:boolean;
    begin
      SEQLEN := 0; TOOLONG:=false;
      (* Read in the sequence name, if present. Name may be preceeded by
         comment lines. *)
      while SEQFILE^=';' do readln(SEQFILE);
      if SEQFILE^='>' then begin READWORD(SEQFILE,NAME);readln(SEQFILE) end;
      
      (* Read in the sequence *)
      while (not eof(SEQFILE)) and (not (SEQFILE^='>')) do begin
        if SEQFILE^=';' then (* comment line, ignore *)
        else while not eoln(SEQFILE) do begin
          read(SEQFILE,CH);
          if CH in ['A','a','C','c','G','g','T','t','U','u','N','n',
                    'R','r','Y','y','M','m','W','w','S','s','K','k',
                    'D','d','H','h','V','v','B','b'] then
            if SEQLEN < MAXSEQ-2 then begin
              SEQLEN := SEQLEN + 1;
              SEQ[SEQLEN]:=NUC(CH);
              end
            else TOOLONG:= true
          end; (* eoln *)
          readln(SEQFILE);
          end; (* eof *)
      if TOOLONG then begin
         writeln(
        '>>> WARNING! Sequence length exceeds MAXSEQ-2. Seq. truncated.')
        end; (* TOOLONG *)
      START:=1; FINISH:=SEQLEN
    end; (* READSEQUENCE *)

  (*************************************************************)
  (*  Translate sequence and write to output.                  *)
  (*************************************************************)
  procedure TRANSLATE(var F:text; var NAME:WORD; var SEQ:SEQUENCE;
                      var START,FINISH:integer);
    const CODONSPERLINE=50;
    var  P1,P2,P3,CODONSLEFT,EXTRANUCS,
         THISLINE,TOTALDISTANCE,
         TOTALCODONS,WIDTH:integer;

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

    begin (* TRANSLATE *)

      if NAME.LEN >0 then begin WRITEWORD(F,NAME,NAME.LEN);writeln(F) end;
      TOTALDISTANCE:=FINISH-START+1;
      TOTALCODONS:=TOTALDISTANCE div 3;
      
      (* If sequence is truncated within last codon, add one or two N's to
         complete the codon. Sequence ambiguity may allow another aa to
         be assigned. *)
      EXTRANUCS:=TOTALDISTANCE mod 3;
      if EXTRANUCS > 0 then begin
        SEQ[SEQLEN+1]:=N;
        if EXTRANUCS=2 then SEQ[SEQLEN+2]:=N;
        TOTALCODONS:=TOTALCODONS+1
        end; (* EXTRANUCS > 0 *)

      (* Translate the sequence *)
      P1:=START;
      CODONSLEFT:=TOTALCODONS;
      while CODONSLEFT > 0 do begin
        if CODONSLEFT>=CODONSPERLINE then THISLINE:=CODONSPERLINE
        else THISLINE:=CODONSLEFT;
        for WIDTH:= 1 to THISLINE do begin
          P2:= P1+1;P3:=P2+1;
          write(F,AAS[AAMATCH(SEQ[P1],SEQ[P2],SEQ[P3])]);
          P1:= P3+1
          end; (* for WIDTH *)
        writeln(F);
        CODONSLEFT:=CODONSLEFT-THISLINE
        end (* CODONSLEFT>0 *)           
    end; (* TRANSLATE *)
 
     
  (* ----------------------------------------------------------  *)
  (* ----------------- MAIN  PROCEDURE  -----------------------  *)
    begin
      (* Set up initial values for nucleotide tables and genetic code *)
      INITTABLES;
      INITSTRING;
      
      (* Read in a genetic code file, if specified by command line option -g*)
      if argc >1 then begin
        ARGNUM:=2;
        argv(ARGNUM,ARGUMENT);
        if (ARGUMENT[1]='-') and (ARGUMENT[2]='g') then begin 
          ARGNUM:=3;
          FILEARGS(GCFILE,'I',ARGNUM);
          READCODE(GCFILE,GC)
          end;      
        end;

      (* Derive ambiguity rules based on genetic code *)
      MAKERULES(GC);
       
      (* Read each sequence and translate while writing to output.*)
      while not eof(input) do begin
        READSEQUENCE(input,SEQNAME,SEQ,START,FINISH);
        TRANSLATE(output,SEQNAME,SEQ,START,FINISH);        
        (* Advance to first line beginning with '>' *)
        while not(input^='>') and not eof(input) do readln(input)
        end  
    end. (* RIBOSOME *)
