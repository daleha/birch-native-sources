  (**********************************************************)
  (*                                                        *)
  (*     INTREST Version  8/14/2001, Standard Pascal        *)
  (*             Brian Fristensky                           *)
  (*             Dept. of Plant Science                     *)
  (*             University of Manitoba                     *)
  (*             Winnipeg, MB R3T 2N2  CANADA               *)
  (*                                                        *)
  (* Copyright (c) 1986,1987,1988,1990 by Brian Fristensky  *)
  (* !!! in comment indicates feature which may need change *)
  (******************************************************** *)
  program INTREST(input, output (*,INFILE,OUTFILE*));
(*!!! Some Pascals may require file parameters in program heading *)
 
  const MAXSEQ  =  750000; 
 
(* BEGIN MODULE INTLIMITS *)
(*!!!  MAXINT =  2147483647; *)
(*!!!  MININT = -2147483647; *)
(* END MODULE INTLIMITS         VERSION= 'SUNMODS     Version  6/26/01'; *)
 
        MAXPAT  = 23;
        MAXFRAGS= 6000; (*max. # of fragmenst in a single digest *)
        MAXWORD = 23;
        MAXLINE = 150;
        VERSION = 'INTREST   Version  8/14/2001';
 
  type NUCLEOTIDE= (CANT,A,C,R,D,V,M,K,B,H,Y,G,T,W,S,N,WONT);
       SS          = set of NUCLEOTIDE;
       SEQUENCE    = array[-MAXPAT..MAXSEQ] of NUCLEOTIDE;
 
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
 
       (* RESTRICTION ENZYME *)
       ENZYME      = record RSEQ: array[1..MAXPAT] of NUCLEOTIDE;
                            PLEN: integer;
                            NAME,SITE: WORD;
                            SYMMETRIC:boolean;
                            CUT,ACUT: integer
                            end;
 
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
(* END MODULE TYPE.FRAG         VERSION= 'SUNMODS     Version  6/26/01'; *)
 
  var INFILE,OUTFILE  : text;
      IFN,OFN         : LINE;
      NAME            : WORD;
      COM    : array[NUCLEOTIDE] of NUCLEOTIDE;(*Complements of NUCLEOTIDE *)
      NUCSET          : array[NUCLEOTIDE] of SS;
      AMBIGUOUS       : SS;
      ENDNUC          : NUCLEOTIDE;
      ENZ             : ENZYME;                  (*Enzyme site to search for *)
      FOUND           : FRAGSFOUND;              (*List of sites found*)
      FREEFRAG        : FRAG;                    (*Points to freelist of sites*)
      SEQ             : SEQUENCE;                (*Sequence to be searched*)
      SEQLEN          : integer;                 (*Sequence length*)
      CIRCULAR        : boolean;                 (*=true if seq. is circular *)
      HLINE           : LINE; 
      I,CHOICE        : integer;
 
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
 
  (************************************************)
  (*  Initialization procedures.                  *)
  (************************************************)
  procedure INITIALIZE;
    var B1,B2:NUCLEOTIDE;
    begin
      (* Array NUCSET holds sets of nucleotides  *)
      (* using the conventions of IUPAC-IUB      *)
      NUCSET[A]:= [A];NUCSET[C]:= [C];NUCSET[G]:= [G];NUCSET[T]:= [T];
      NUCSET[R]:= [A,G,R];NUCSET[Y]:= [C,T,Y];NUCSET[S]:= [C,G,S];
      NUCSET[M]:= [A,C,M];NUCSET[K]:= [G,T,K];NUCSET[W]:= [A,T,W];
      NUCSET[B]:= [C,G,T,B,Y,K,S];NUCSET[D]:= [A,G,T,W,R,D,K];
      NUCSET[H]:= [A,C,T,M,Y,W,H];NUCSET[V]:= [A,C,G,M,R,S,V];
      NUCSET[N]:= [A..N];
      AMBIGUOUS:= NUCSET[N] - [A,C,G,T];
      (* Array  COM holds the symbol values of nucleotide complements.*)
      (* A complements T, C complements G, R complements Y etc.       *)
      B1:= A; B2:= T;
      while B1<= T do begin
        COM[B1]:= B2;
        B1:= succ(B1); B2:= pred(B2)
        end;
      (* N,S, and W complement themselves *)
      for B1:= W to N do COM[B1]:= B1;
      (* Initialize the linked list FOUND *)
      FREEFRAG:= nil;
      new(FOUND.HEAD); new(FOUND.TAIL);
      FOUND.HEAD^.NEXT:= FOUND.TAIL;FOUND.TAIL^.PREV:= FOUND.HEAD
    end; (* INITIALIZE *)

    (*******************************************************************)
    (* If sequence is circular, copy the last MAXPAT bases to the start*)
    (*******************************************************************)
    procedure INITENDS;
      var Km,Sk:integer;
      begin
      if CIRCULAR then begin
        Km:= 1-MAXPAT;
        for Sk:= SEQLEN-(MAXPAT-1) to SEQLEN do begin
          SEQ[Km]:= SEQ[Sk]; Km:= Km + 1 end
        end;
      end; (* INITENDS *)

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

 
  (*********************************************************)
  (*  Read in a restriction enzyme site from the terminal. *)
  (*********************************************************)
  procedure READSITE(var ENZ:ENZYME);
 
  var LOC: integer;
      LEGAL:boolean;
 
  begin  (* READSITE *)
    with ENZ do begin
      repeat
        LEGAL:= true;
        writeln('Enter name of enzyme:');
        INPWORD(NAME);readln;
        writeln('Enter recognition site using the symbols below:');
        writeln('Nucleotides:  A,C,G,T');
        writeln('Ambiguous nucleotides:');
        writeln('R = A,G    M = A,C    B = C,G,T   N = A,C,G,T');
        writeln('Y = C,T    K = G,T    D = A,G,T   V = A,C,G');
        writeln('S = C,G    W = A,T    H = A,C,T');
        writeln;
        writeln('Site must be <=',MAXPAT-1,' characters:');
        INPWORD(SITE);readln;
        (* Check for illegal symbols *)
        with SITE do for LOC:= 1 to LEN do
          if STR[LOC] in ['A','a','C','c','R','r','D','d','V','v','M','m',
             'K','k','B','b','H','h','Y','y','G','g','T','t','N','n','S','s',
             'W','w'] then RSEQ[LOC]:=NUC(STR[LOC])
          else begin
            writeln('Illegal character  ',STR[LOC]);
            LEGAL:= false
            end (* else *)
      until LEGAL;
      PLEN:= SITE.LEN;
      writeln('Position of cut?');
      GETINTEGER(CUT,-MAXINT,MAXINT);readln
    end (* with ENZ *)
  end; (* READSITE *)
 
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
(* END MODULE LINKED         VERSION= 'SUNMODS     Version  6/26/01'; *)
 
  (****************************************************************)
  (*  Search for site ENZ in SEQ, using the algorithm found in    *)
  (*    "Fast Pattern Matching in Strings"                        *)
  (*     Knuth, Morris, and Pratt, SIAM Journal of Computing,     *)
  (*     Vol.6, No.2, June 1977.  pp323-350.                      *)
  (****************************************************************)
  procedure KMP(var ENZ:ENZYME; var OUTFILE:text);
 
    var NEXT,GOBACK :array[1..MAXPAT] of integer;
        PATTERN :array[1..MAXPAT] of SS;
        Pj,Cj,PLEN1: integer;
        RIGHTMOST: FRAG;
        CRP: NUCLEOTIDE;
 
    (* The NEXT array is used when  PATTERN[Px] does not match  STR[Sk].*)
    (* After a mismatch is found, instead of backing up to PATTERN[1],  *)
    (* the search resumes at PATTERN[NEXT[Pj]].                         *)
    procedure MAKENEXT;
      label  1,2;
      var Pj,FPj,LEFTAM:integer;
      begin
        with ENZ do begin
          PLEN1:= PLEN+1;
          (* Ambiguous nucleotide in first position is a special case *)
          if (AMBIGUOUS * PATTERN[1]) <> [] then begin
            NEXT[1]:= 0; GOBACK[1]:=0;
            for Pj:= 2 to PLEN1 do begin
              NEXT[Pj]:=1; GOBACK[Pj]:= Pj-2 end
            end
          (*Set LEFTAM equal to position of leftmost ambiguous symbol*)
          (* after the first position in pattern.                    *)
          else begin
            LEFTAM:=PLEN1;Pj:=2;
            while Pj < PLEN1 do begin
              if (AMBIGUOUS * PATTERN[Pj]) <> [] then begin
                LEFTAM:= Pj; goto 1 end;
              Pj:= Pj+1
              end;
            (*Calculate NEXT & GOBACK tables*)
          1:Pj:=1; FPj:=0; NEXT[1]:=0;GOBACK[1]:=0;
            while Pj < PLEN1 do begin  (*FPj = F[Pj] *)
              2: if FPj > 0 then if not(PATTERN[Pj] <= PATTERN[FPj]) then
                   begin FPj:= NEXT[FPj]; goto 2 end;
              FPj:= FPj+1; Pj:= Pj+1;
              if Pj <= LEFTAM then begin
                if PATTERN[Pj] <= PATTERN[FPj]
                  then NEXT[Pj]:= NEXT[FPj]
                  else NEXT[Pj]:= FPj;
                GOBACK[Pj]:= 0
                end
              else begin
                NEXT[Pj]:= NEXT[LEFTAM];
                GOBACK[Pj]:= Pj-LEFTAM
                end
              end;
            GOBACK[PLEN1]:= Pj-LEFTAM
          end
          end
        end;  (*MAKENEXT*)
 
    (* Search the sequence.*)
    procedure SEARCH(CUTSITE:integer);
      label 10,20,30,99;
      var Pj,Sk,RESUME :integer;
          ENDSET:SS;
 
    (*  Add a node to the linked list FOUND, telling the cutting site.*)
    procedure PATTERNMATCHED;
      var NEWFRAG:FRAG;
          POSN:integer;
      begin
        with ENZ do begin
          POSN:= (Sk-PLEN)+CUTSITE+GOBACK[PLEN1];
          if CIRCULAR then
            if POSN < 1 then POSN:= SEQLEN + POSN
            else if POSN > SEQLEN then POSN:= POSN-SEQLEN;
          if (POSN >= 1) and (POSN <= SEQLEN) then begin
            GETFRAG(NEWFRAG);
            NEWFRAG^.START:= POSN;
            if RIGHTMOST^.START > POSN then RIGHTMOST:= FOUND.HEAD;
            while RIGHTMOST^.NEXT^.START < POSN do
              RIGHTMOST:= RIGHTMOST^.NEXT;
            ADDFRAG(NEWFRAG,RIGHTMOST);
            RIGHTMOST:= NEWFRAG;
            FOUND.LNUM:= FOUND.LNUM+1
            end
          end
      end;  (* PATTERNMATCHED *)
 
      begin (* SEARCH *)
        with ENZ do begin
          (* Initialize search parameters *)
          ENDSET:= PATTERN[1];
          SEQ[SEQLEN+1]:= WONT;(*WONT will never be found in CANT set*)
          RESUME:= NEXT[PLEN+1];(* Resume search at this position after match*)
          NEXT[PLEN+1]:= -1;
          SEQ[SEQLEN+2]:=ENDNUC;(*BASE value of PATTERN[1]*)
          (*Begin at leftmost chars of PATTERN and STR*)
          Pj:=1;
          if CIRCULAR then Sk:= 1-(PLEN-1) else Sk:= 1;
 
      10: (* GETSTARTED, Pj = 1 *)
          (*Advance until first match *)
          while not (SEQ[Sk] in ENDSET) do Sk:=Sk+1;
          if Sk > SEQLEN then goto 99; (*INPUTEXHAUSTED*)
 
      20: (* CHARMATCHED *) Pj:= Pj+1; Sk:= Sk+1;
 
      30: (* LOOP, Pj > 0 *)
          if SEQ[Sk] in PATTERN[Pj] then goto 20; (*CHARMATCHED*)
 
          Sk:= Sk-GOBACK[Pj]; Pj:= NEXT[Pj]; (* char not matched *)
          if Pj = 1 then goto 10; (*GETSTARTED*)
          if Pj = 0 then begin
               Pj:=1; Sk:= Sk+1; goto 10 (*GETSTARTED*) end;
          if Pj > 0 then goto 30; (*LOOP*)
 
          (* Pj = -1, PATTERN matched *)
          PATTERNMATCHED;
          Pj:= RESUME; goto 30; (*LOOP*)
 
      99: (* INPUTEXHAUSTED *) ;
      end
    end; (*SEARCH*)
 
 
    begin (* KMP *)
      FOUND.LNUM:=0;FOUND.TAIL^.START:=MAXSEQ+1;(*Always > POSN of site found*)
      with ENZ do begin
  
        (* Search SEQ using PATTERN derived from RSEQ  *)
        ENDNUC:= RSEQ[1];
        for Pj:= 1 to PLEN do PATTERN[Pj]:= NUCSET[RSEQ[Pj]];
        PATTERN[PLEN+1]:= [CANT]; (* Can not be matched *)
        RIGHTMOST:= FOUND.HEAD;
        MAKENEXT;
        SEARCH(CUT);
 
        (* Let PATTERN = the inverse complement of PATTERN*)
        (* Check for symmetry.                            *)
        Cj:= PLEN; SYMMETRIC:= true;
        for Pj:= 1 to PLEN do begin
          CRP:= COM[RSEQ[Pj]];
          PATTERN[Cj]:= NUCSET[CRP];
          if RSEQ[Cj] <> CRP then SYMMETRIC:= false;
          Cj:= Cj-1
          end;
        if not SYMMETRIC then begin
          writeln('Asymmetric Site!');
          writeln('Type position of cut on opposite strand');
          GETINTEGER(ACUT,-MAXINT,MAXINT);readln;
          ENDNUC:= CRP;
          MAKENEXT;
          RIGHTMOST:= FOUND.HEAD;
          SEARCH(-(ACUT-PLEN))
          end (* not SYMMETRIC *)
      end (* with ENZ *)
    end; (*KMP*)
 
 
    (*********************************************)
    (* Compile a report for output.              *)
    (*********************************************)
    procedure REPORT(var OUTFILE:text);
    var ORDER: array[1..MAXFRAGS] of FRAG;
        TOOMANY_FRAGS: boolean; (* =true if SITES > MAXFRAGS *)
        NUMSITES:integer;
    (*Calculate sizes and ends of fragments*)
    procedure CALCULATE;
 
      var CURRENTFRAG:FRAG;
          SITE:integer;
      begin  (* CALCULATE *)
        with FOUND do begin
          NUMSITES:=LNUM;
        (* if LINEAR add a fragment to head of list*)
        if (not CIRCULAR) and (HEAD^.NEXT^.START > 1) then begin
          GETFRAG(CURRENTFRAG);
          ADDFRAG(CURRENTFRAG,HEAD);
          CURRENTFRAG^.START:= 1;
          LNUM:= LNUM+1
          end
        else CURRENTFRAG:= HEAD^.NEXT;
 
        if LNUM > 0 then begin
        (*Calculate ends and size of each fragment and assign it *)
        (* a place in the ORDER array, to be sorted later.    *)
          SITE:=1;
          while CURRENTFRAG^.NEXT <> TAIL do begin
            with CURRENTFRAG^ do begin
              FINISH:= NEXT^.START-1; SIZE:= FINISH-START+1 end;
            if (SITE <= MAXFRAGS) then ORDER[SITE]:= CURRENTFRAG
            else TOOMANY_FRAGS:=true;  
            SITE:= SITE + 1;                   
            CURRENTFRAG:= CURRENTFRAG^.NEXT
            end;
          (*Last fragment in list is a special case*)
          with CURRENTFRAG^ do
            if CIRCULAR then begin
              FINISH:= HEAD^.NEXT^.START-1;
              SIZE:= (SEQLEN-START)+(FINISH+1);
              if FINISH=0 then FINISH:= SEQLEN end
            else begin
              FINISH:= SEQLEN; SIZE:= FINISH-START+1 end;
            if (SITE <= MAXFRAGS) then ORDER[SITE]:= CURRENTFRAG
            else TOOMANY_FRAGS:=true             
          end
     end
   end; (* CALCULATE *)  
 
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
(* END MODULE SORT         VERSION= 'SUNMODS     Version  6/26/01'; *)
 
    (*  Print the list of sites found.        *)
    procedure PRINTLIST;
      var THISFRAG: FRAG;
          SITE:integer;
      begin (* PRINTLIST *)
        with FOUND do begin
        (* Write a heading for the site *)
          with ENZ do begin
            WRITEWORD(OUTFILE,NAME,10);
            WRITEWORD(OUTFILE,SITE,MAXPAT);
            write(OUTFILE,CUT:3);
            if SYMMETRIC then write(OUTFILE,' ':6)
            else write(OUTFILE,' (',ACUT:3,')') 
            end;
          writeln(OUTFILE,NUMSITES:5);
          if NUMSITES = LNUM then THISFRAG:= HEAD
          else THISFRAG:= HEAD^.NEXT;
          SITE:=1;
          while THISFRAG^.NEXT <> TAIL do begin
            THISFRAG:= THISFRAG^.NEXT;
            write(OUTFILE,THISFRAG^.START:MAXPAT+30);
            with ORDER[SITE]^ do writeln(OUTFILE,SIZE:8,START:8,FINISH:8);
            SITE:= SITE+1
            end;
          (* Print last fragment. If no sites are found, a linear *)
          (* sequence will still yeild one fragment.              *)
          if NUMSITES < LNUM then with ORDER[SITE]^ do
            writeln(OUTFILE,SIZE:MAXPAT+38,START:8,FINISH:8)
        end;
        writeln(OUTFILE)
      end; (* PRINTLIST *)
 
 
  begin (* REPORT *)
    TOOMANY_FRAGS:=false;
    CALCULATE;
    if TOOMANY_FRAGS then begin
       writeln('>>> Number of fragments exceeds ',MAXFRAGS);
       writeln('>>> Increase MAXFRAGS and recompile');
       writeln(OUTFILE);
       writeln(OUTFILE, '>>> Number of fragments exceeds ',MAXFRAGS);
       writeln(OUTFILE, '>>> Increase MAXFRAGS and recompile');
       end
    else begin
      BUBBLESORT(FOUND.LNUM,1); 
      PRINTLIST
    end;
    with FOUND do RIDOF(HEAD,TAIL)
  end;  (*REPORT*)
 
  (***************************************************************)
  (* Prompt user for sites and call SEARCH.                      *) 
  (***************************************************************)
  procedure PROMPT(var OUTFILE:text);
    var ANSWER:char;
    begin
      (* Write heading *)
      writeln(OUTFILE,VERSION);
      WRITEWORD(OUTFILE,NAME,NAME.LEN);
      write(OUTFILE,'  Configuration: ');
      case CIRCULAR of
        true:write(OUTFILE,' CIRCULAR');
        false:write(OUTFILE,' LINEAR')
        end;
      writeln(OUTFILE,'  Length: ',SEQLEN,' bp');
      writeln(OUTFILE,'# of':MAXPAT+24);
      writeln(OUTFILE,'Cut     Sites Sites   Frags   Begin     End':MAXPAT+54);
      repeat
        READSITE(ENZ);
        KMP(ENZ,OUTFILE);
        REPORT(OUTFILE);
        repeat
          writeln('Type  S to search for a site, Q to return to main menu:');
          readln(ANSWER)
        until ANSWER in ['S','s','Q','q']
      until ANSWER in ['Q','q'];
    end; (* PROMPT *) 


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
      INITIALIZE;
      (* Initialize horizontal output line *)
      with HLINE do begin
        for I:= 1 to MAXLINE do STR[I]:='_';
        LEN:=MAXLINE
        end;
      
      (* Read in initial sequence and set up arrays. *)
      writeln('Enter sequence filename:');
      GETFILE(INFILE,'I',IFN);
      NAME.LEN:=0;
      READSEQ(INFILE,SEQ,SEQLEN,NAME,CIRCULAR);
      if NAME.LEN=0 then begin
        writeln('Type name to appear on output:');
        INPWORD(NAME);readln
        end;
      INITENDS;

      OFN.LEN:=0; (* indicates that OUTFILE is not yet open *)

      (* MAIN MENU *)
      repeat
        writeln;
        WRITELINE(output,HLINE,80);writeln;
        writeln('INTREST','MAIN MENU':30);
        WRITELINE(output,HLINE,80);writeln;
        write('Input file:        ');WRITELINE(output,IFN,IFN.LEN);writeln;
        write('Output file:       ');WRITELINE(output,OFN,OFN.LEN);writeln;
        WRITELINE(output,HLINE,80);writeln;
        writeln(' ':20,'1) Read in a new sequence');
        writeln(' ':20,'2) Open a new output file');
        writeln(' ':20,'3) Search for sites (output to screen)');
        writeln(' ':20,'4) Search for sites (output to file)');
        WRITELINE(output,HLINE,80);writeln;
        writeln('Type the number of your choice  (0 to quit program)');
        GETINTEGER(CHOICE,0,4);readln;

        case CHOICE of
          0:;
          1:begin
              writeln('Enter sequence filename:');
              GETFILE(INFILE,'I',IFN);
              NAME.LEN:=0;
              READSEQ(INFILE,SEQ,SEQLEN,NAME,CIRCULAR);
              if NAME.LEN=0 then begin
                writeln('Type name to appear on output:');
                INPWORD(NAME);readln
                end;
              INITENDS
            end;
          2:begin
(*!!!         CLOSE(OUTFILE); *)
              writeln('Type output filename:');
              GETFILE(OUTFILE,'O',OFN);
            end;
          3:PROMPT(output);
          4:if OFN.LEN>0 then PROMPT(OUTFILE)
            else begin writeln('>>> NO OUTPUT FILE CURRENTLY OPEN');
                       readln end
          end (* case *)
        until CHOICE = 0;
(*!!!*)CLOSE(OUTFILE)
    end. (* INTREST *)
 
