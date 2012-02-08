  (**********************************************************)
  (*                                                        *)
  (*    BACHREST Version  3/29/2006,  Standard Pascal       *)
  (*             Brian Fristensky                           *)
  (*             Dept. of Plant Science                     *)
  (*             University of Manitoba                     *)
  (*             Winnipeg, MB R3T 2N2  CANADA               *)
  (*                                                        *)
  (* Copyright (c) 1986, 1987, 1990  by Brian Fristensky    *)
  (* !!! in comment indicates feature which may need change *)
  (******************************************************** *)
(* REVISION HISTORY 
26 Mar 2006 Added PROT5, BLUNT and PROT3 parameters.
24 Mar 2006 Added SYMM and FRAGPRINT parameters.
20 Mar 2006 A digest will only be printed if the number of fragments
            is such that FRAGLEAST <= # of fragments <= FRAGMOST.
17 Mar 2006 Increased MAXSEQ to 20,000,000, which should handle 
            all prokaryotic and many eukaryotic chromosomes.
11 Aug 2001 Bachrest rebuild using readseq, modified for GenBank
            locus names of 18 char. Improved error checking for enzyme
            lists exceeding MAXFRAGS. Increased MAXFRAGS for working
            with larger sequences.
5 Mar 1998  Width of enzyme name increased to 15 to accommodate long
            names in REBASE. Also, a mandatory blank is included between
            enzyme name and recognition seq. in output, to make sure
            that DIGEST can still read them, if the name was ever 15
            characters or longer. Columns had to be adjusted acordingly
            for other output fields.
*)
  program BACHREST(input, output (*,INFILE,RESTFILE,OUTFILE*));
(*!!! Other Pascals may require file parameters in heading *)
  const MAXSEQ  =  20000000;
 
(* BEGIN MODULE INTLIMITS *)
(*!!!  MAXINT =  2147483647; *)
(*!!!  MININT = -2147483647; *)
(* END MODULE INTLIMITS         VERSION= 'SUNMODS     Version  6/26/01'; *)
 
        MAXPAT  = 23;
        MAXFRAGS= 6000;
        MAXWORD = 23;
        MAXLINE = 150;
        VERSION = 'BACHREST   Version  3/29/2006';
 
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
       ENZYME      = record NAME,PROTONAME: WORD;
                            RECSTR: WORD; (* string value of recog. seq.*)
                            PLEN: integer;
                            RSEQ: array[1..MAXPAT] of NUCLEOTIDE;
			    SYMMETRIC:boolean;
                            CUT,CUTOPP: integer;
			    FRAGENDS:char;
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
 
  var INFILE,RESTFILE,OUTFILE  : text;
      IFN,RFN,OFN:LINE; (* File names *)
      REBASE:boolean; (*=true if RESTFILE is REBASE, false if FSAP format *)
      NAME:WORD;
      COM    : array[NUCLEOTIDE] of NUCLEOTIDE;(*Complements of NUCLEOTIDE *)
      NUCSET           : array[NUCLEOTIDE] of SS;
      LEGALNUCS       : set of char;
      AMBIGUOUS       : SS;
      ENDNUC          : NUCLEOTIDE;
      ENZ             : ENZYME;                  (*Enzyme site to search for *)
      FOUND           : FRAGSFOUND;              (*List of sites found*)
      FREEFRAG        : FRAG;                    (*Points to freelist of sites*)
      SEQ             : SEQUENCE;                (*Sequence to be searched*)
      SEQLEN          : integer;                 (*Sequence length *)
      CIRCULAR        : boolean;                 (*=true if seq. is circular *)
      TEMPLINE,                                  (* dummy input line *)
      HLINE           : LINE;                    (* horizontal line for menus *)
      I,J,
      FIRSTITEM,                                 (* line number of first REBASE entry *)       
      CHOICE          : integer;                 (* # of menu choice *)

      (* Global search parameters *)
      SOURCE,                 (* Search for commercial or all *)
      PROTOTYPE,              (* Search for prototype, or all isoscizomers *)
      PROT3,                  (* Search for 3' protruding end sites *)
      BLUNT,                  (* Search for blunt end sites *)
      PROT5,                  (* Search for 5' protruding end sites *)
      SYMM:char;              (* Search for symmetric, assym. or both *)
      MINSITE,                (* Minimum length of RE seq. to search for *)
      MAXSITE:integer;        (* Maximum length of RE seq. to search for *)
      FRAGLEAST,              (* Minimum number of fragments in a digest *)
      FRAGMOST:integer;       (* Maximum number of fragments in a digest *)
      FRAGPRINT: integer;       (* Maximum number of fragments to print *)
      FITSPARAMS:boolean;     (* Last enzyme read fits parameter settings *)
                              (* for type of RE to search for*)


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


(* BEGIN MODULE SKIPBLANKS *)
   (* Skip to next non-blank character on the current line. *)
   procedure SKIPBLANKS(var F:text);
     var DONE:boolean;
         CH:char; 
     begin
       DONE:=false;
       repeat
         if (eoln(F) or eof(F)) then DONE:= true
                                else if F^ = ' ' then read(F,CH)
                                                 else DONE:=true
       until DONE   
     end; (* SKIPBLANKS *)
(* END MODULE SKIPBLANKS *)


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
(* END MODULE READLINE         VERSION= 'SUNMODS     Version  6/26/01'; *)

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
 
  (************************************************)
  (*  Initialize arrays.                          *)
  (************************************************)
  procedure INITIALIZE;
    var B1,B2:NUCLEOTIDE;
    begin
      (* Legal nucleotide characters *)
      LEGALNUCS:= ['A','a','C','c','R','r','D','d','V','v','M','m',
             'K','k','B','b','H','h','Y','y','G','g','T','t','N','n','S','s',
             'W','w']; 
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

      (* Initialize search parameters *)
      SOURCE:='C';  PROTOTYPE:='A'; 
      PROT3:='Y'; BLUNT:='Y'; PROT5:='Y';
      SYMM:='B';
      MINSITE:=4; MAXSITE:=MAXPAT;
      FRAGLEAST:=0; FRAGMOST:=MAXFRAGS; FRAGPRINT:= 30;
 
      (* Initialize the linked list FOUND *)
      FREEFRAG:= nil;
      new(FOUND.HEAD); new(FOUND.TAIL);
      FOUND.HEAD^.NEXT:= FOUND.TAIL;FOUND.TAIL^.PREV:= FOUND.HEAD;

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

    (*******************************************************************)
    (* Print a header to OUTFILE.                                      *)
    (*******************************************************************)
    procedure HEADER(var OUTFILE:text);
      const HLINE = '-----------------------------------------------------------';
      begin
        writeln(OUTFILE,HLINE);
        writeln(OUTFILE,VERSION:50);
	writeln(OUTFILE,'');
        WRITEWORD(OUTFILE,NAME,18);
        write(OUTFILE,'  Topology: ');
        case CIRCULAR of
          true:write(OUTFILE,'CIRCULAR');
          false:write(OUTFILE,'LINEAR')
          end;
        writeln(OUTFILE,'  Length: ',SEQLEN:8,' bp');
        writeln(OUTFILE,HLINE);	
	writeln(OUTFILE,'Search parameters:');
	writeln(OUTFILE,'   Recognition sequences between ',MINSITE:4,' and ',MAXSITE:4, ' bp');
	
	write(OUTFILE,'   Ends: ');
	if PROT5 = 'Y' then write(OUTFILE,'5'' protruding');
	if ((PROT5 = 'Y') and ((BLUNT = 'Y') or (PROT3 = 'Y'))) then write(OUTFILE,', ');
	if BLUNT = 'Y' then write(OUTFILE,'Blunt');
	if (BLUNT = 'Y') and (PROT3 = 'Y') then write(OUTFILE,', ');
	if PROT3 = 'Y' then write(OUTFILE,'3'' protruding');
	writeln(OUTFILE);
	
	write(OUTFILE,'   Type: ');
        if (SYMM = 'S') then writeln(OUTFILE,'Symmetric')		
        else if (SYMM = 'B') then writeln(OUTFILE,'Symmetric, Asymmetric')
	else writeln(OUTFILE,'Asymmetric');
	
	write(OUTFILE,'   Minimum fragments: ' ,FRAGLEAST:5,'     ');
	writeln(OUTFILE,'Maximum fragments: ',FRAGMOST:5);
	writeln(OUTFILE,'   Maximum fragments to print: ',FRAGPRINT:5);
		
        writeln(OUTFILE,HLINE);	
	writeln(OUTFILE,'');
        writeln(OUTFILE,'# of':MAXPAT+22);
        writeln(OUTFILE,'Enzyme          Recognition Sequence     Sites     Sites   Frags   Begin     End');
        writeln(OUTFILE,HLINE,'---------------------');
	writeln(OUTFILE,'');	
      end; (* HEADER *)


  (**************************************************************)
  (* Determine whether RESTFILE is in REBASE or FSAP format.    *)
  (* if REBASE = true, then RESTFILE^ points to the first entry.*)
  (* if REBASE = false, the end of file will be reached.        *)
  (**************************************************************)
  procedure CHECKFORMAT(var RESTFILE:text; var REBASE:boolean); 
    var I:integer;
    begin
      I:=1;
      while (RESTFILE^ <> '<') and (not eof(RESTFILE)) do begin
         readln(RESTFILE);
         I:=I+1
         end;
      if eof(RESTFILE) then REBASE:=false (* FSAP format *)
      else REBASE:=true; (* REBASE *)
         
    end; (* CHECKFORMAT *)

  (******************************************************)
  (*  Read in a restriction enzyme site from RESTFILE.  *)
  (******************************************************)
  procedure READSITE(var RESTFILE:text; var ENZ:ENZYME);
 
  var LOC: integer;
      LEGAL:boolean;
 
  begin
    with ENZ do begin
      repeat
        LEGAL:= true;
        (* Read enzyme name, recognition site, and cutting site from  file *)
        READWORD(RESTFILE,NAME);
        READWORD(RESTFILE,RECSTR);
        read(RESTFILE,CUT);
        SKIPBLANKS(RESTFILE);
        if eoln(RESTFILE) then readln(RESTFILE)
        else if not eof(RESTFILE) then begin
                     read(RESTFILE,CUTOPP);
                     if eoln(RESTFILE) then readln(RESTFILE)
                     end;

        (* Check for illegal symbols *)
        for LOC:= 1 to RECSTR.LEN do
          if RECSTR.STR[LOC] in LEGALNUCS then RSEQ[LOC]:=NUC(RECSTR.STR[LOC])
          else LEGAL:= false;
          if not LEGAL then writeln(OUTFILE,'>>>> Illegal input ignored')     
      until LEGAL;
      PLEN:= RECSTR.LEN
 

    end
  end; (* READSITE *)

  (******************************************************)
  (*  Read in a restriction enzyme site from RESTFILE.  *)
  (******************************************************)
  procedure READREBASE(var RESTFILE:text; var ENZ:ENZYME;
                       var FITSPARAMS:boolean);
 
  var COMSOURCE:WORD;
      CARET:boolean; (* G^AATTC format, rather than GAAGA(8/7) *)
      POSN,I,ORDZERO: integer;
      CH:char;

  (* Extract an integer from a WORD *)
  procedure GETNUM(W:WORD;var POSN,NUM:integer);
    var SIGN:integer;
    begin
      with W do begin
        SIGN:=1; NUM:=0;
        if STR[POSN] = '-' then begin
           SIGN:=-1;
           POSN:=POSN+1
           end;
        while (STR[POSN] in ['0'..'9']) and (POSN<=LEN) do begin
           NUM:= (10*NUM) + (ord(STR[POSN])-ORDZERO);
           POSN:=POSN+1
           end;
        if SIGN=-1 then NUM:= -NUM
      end (* with W *)
    end; (* GETNUM *)
 
  begin
    with ENZ do begin
      ORDZERO:=ord('0');
      CARET:=true;

        (* Name begins at 4th char. or 1st line of a record *)
           for I:= 1 to 3 do read(RESTFILE,CH);
           READWORD(RESTFILE,NAME);
           readln(RESTFILE);

        (* Prototype begins at 4th char. on 2nd line of a record *)
           for I:= 1 to 3 do read(RESTFILE,CH);
           if not eoln(RESTFILE) then READWORD(RESTFILE,PROTONAME)
           else PROTONAME.LEN:=0;
           readln(RESTFILE);
           if (PROTOTYPE='P') and (PROTONAME.LEN >0) then FITSPARAMS:=false;

        (* Site begins at 4th char. of line 3 *)
           for I:= 1 to 3 do read(RESTFILE,CH);
           READWORD(RESTFILE,RECSTR);
           readln(RESTFILE);

        (* Extract cut sites from RECSTR, and write the site to RSEQ *)
           POSN:= 1;PLEN:=0;

           CUT:=0; (* if no cut site specified, assume zero *)
           while POSN <= RECSTR.LEN do begin
               if RECSTR.STR[POSN] in LEGALNUCS then begin (* recognition seq. *)
                  PLEN:= PLEN+1;
                  RSEQ[PLEN]:= NUC(RECSTR.STR[POSN]);
                  POSN:=POSN+1
                  end
               else if RECSTR.STR[POSN] = '^' then begin (* symmetric cut *)
                 CUT:=POSN-1;
                 CUTOPP:=CUT;
                 POSN:=POSN+1
                 end
               else if RECSTR.STR[POSN] = '(' then begin (* asymmetric cut *)
                 CARET:=false;
                 POSN:= POSN+1;
                 GETNUM(RECSTR,POSN,CUT);
                 POSN:=POSN+1;
                 GETNUM(RECSTR,POSN,CUTOPP);
                 end
               else POSN:=POSN+1
               end; (* while *)

         (* To keep the specification of cut sites consistent, cut sites
         specified using caret (eg. G^AATTC) are converted to coordinates
         consistent with parenthesis notation (eg. GAAGA(8/7)). In caret
         notation, CUT is with reference to the position 5' to the start
         of the recognition sequence. In parenthesis notation, CUT is
         with reference to the 3' end of the sequence. *)
         if CARET then begin
            CUT:= -(PLEN-CUT);
            CUTOPP:=CUT
            end;     


        readln(RESTFILE); (* skip methylation data *)
       
        (* Read list of commercial suppliers. *)
           for I:= 1 to 3 do read(RESTFILE,CH);
           if not eoln(RESTFILE) then READWORD(RESTFILE,COMSOURCE)
           else COMSOURCE.LEN:=0;
           readln(RESTFILE);
           if (SOURCE='C') and (COMSOURCE.LEN =0) then FITSPARAMS:=false;
        readln(RESTFILE); (* skip reference data *)
        while (RESTFILE^ <> '<' ) and (not eof(RESTFILE)) do readln(RESTFILE)

    end (* with ENZ *)
  end; (* READREBASE *)
 

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
      WRITELINE(output,HLINE,80);writeln;
      writeln(' ':12,'Parameter   Description/Response                     Value');
      WRITELINE(output,HLINE,80);writeln;
      writeln(' ':12,' 1)SOURCE    C:Commercial only  A:all',SOURCE:20);
      writeln(' ':12,' 2)PROTOTYPE P:prototypes only  A:all',PROTOTYPE:20);
      writeln(' ':12,' 3)PROT3     3'' protruding end cutters (Y/N)',PROT3:13);
      writeln(' ':12,' 4)BLUNT     Blunt end cutters         (Y/N)',BLUNT:13);
      writeln(' ':12,' 5)PROT5     5'' protruding end cutters (Y/N)',PROT5:13);
      writeln(' ':12,' 6)SYMM      S:symmetric A:asymmetric B:both',SYMM:13);
      writeln(' ':12,' 7)MINSITE   Minimum RE site length           ',MINSITE:11);
      writeln(' ':12,' 8)MAXSITE   Maximum RE site length           ',MAXSITE:11);
      writeln(' ':12,' 9)FRAGLEAST Min. # of fragments to print a digest',FRAGLEAST:7);
      writeln(' ':12,'10)FRAGMOST  Max. # of fragments to print a digest',FRAGMOST:7);
      writeln(' ':12,'11)FRAGPRINT Print a number if > FRAGPRINT frags',FRAGPRINT:9);
      WRITELINE(output,HLINE,80);writeln;
      end; (* DISPLAY *)

    procedure REBASEMSG;
      begin
        writeln('>>> This parameter has no effect unless REBASE is used');
	writeln('>>> as the restriction site file. Press ENTER to continue');
	readln
      end; (* REBASEMSG *)

      begin
        (* Prompt user for new parameter values *)
        repeat
          page(output);
          DISPLAY;
          writeln('Type number of parameter you wish to change ',
            '(0 to continue)');
          GETINTEGER(RESPONSE,0,11);
          readln;
          if RESPONSE in [1..11] then
            case RESPONSE of
               1:if REBASE then GETCHAR(SOURCE,'SOURCE   ',['C','A'])
	         else REBASEMSG;
               2:if REBASE then GETCHAR(PROTOTYPE,'PROTOTYPE ',['P','A'])
	         else REBASEMSG;
               3:GETCHAR(PROT3,'PROT3     ',['Y','N']);
               4:GETCHAR(BLUNT,'BLUNT     ',['Y','N']);
               5:GETCHAR(PROT5,'PROT5     ',['Y','N']);	       
               6:GETCHAR(SYMM,'SYMM      ',['S','A','B']); 
               7:GETNUMBER(MINSITE,'MINSITE   ',4,MAXSITE);
               8:GETNUMBER(MAXSITE,'MAXSITE   ',MINSITE,MAXPAT);
	       9:GETNUMBER(FRAGLEAST,'FRAGLEAST   ',0,FRAGMOST);
               10:GETNUMBER(FRAGMOST,'FRAGMOST   ',FRAGLEAST,MAXFRAGS);
	       11:GETNUMBER(FRAGPRINT,'FRAGPRINT  ',1,FRAGMOST)
              end
        until RESPONSE= 0
     end; (* PARAMETERS *)
 

  (**************************************************************)
  (* Determine whether the recognition sequence is symmetric .  *)
  (**************************************************************)
  procedure TESTSYMMETRY(var ENZ:ENZYME); 
    var Cj,Pj:integer;
        CRP:NUCLEOTIDE;
    begin
      with ENZ do begin
           SYMMETRIC:= true;   
           Cj:= PLEN; 
           for Pj:= 1 to PLEN do begin
             CRP:= COM[RSEQ[Pj]];
             if RSEQ[Cj] <> CRP then SYMMETRIC:= false;
             Cj:= Cj-1
             end      
        end         
    end; (* TESTSYMMETRY *)

  (**************************************************************)
  (* Determine the type of cutting site.                        *)
  (**************************************************************)
  procedure TESTCUT(var ENZ:ENZYME); 
    var CUTSAFTER:integer;
        
    begin
      with ENZ do begin
        if SYMMETRIC then begin
	   if REBASE then CUTSAFTER:= PLEN+CUT
	   else CUTSAFTER:= CUT;
	   if CUTSAFTER < (PLEN/2) then FRAGENDS:= '5'
	   else if CUTSAFTER > (PLEN/2) then FRAGENDS:= '3'
	   else FRAGENDS:= 'B'
	   end
	else begin
	   if CUT < CUTOPP then FRAGENDS:= '5'
	   else if CUT > CUTOPP then FRAGENDS:= '3'
	   else FRAGENDS:= 'B'
	   end   
        end         
    end; (* TESTCUT *)
 
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
  (*                                                              *)
  (*     The algorithm works in O(m+n) units of time, where n is  *)
  (*     the length of the sequence searched, and m is the length *)
  (*     of the pattern to search for.                            *)
  (****************************************************************)
  procedure KMP(var ENZ:ENZYME);
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
            NEXT[1]:=0; GOBACK[1]:=0;
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
        end; (* MAKENEXT *)
 
    (* Search the sequence.*)
    procedure SEARCH(CUTSITE:integer);
      label 10,20,30,99;
      var Pj,Sk,RESUME :integer;
          ENDSET:SS;
 
    (*  Add a node to the linked list FOUND, telling the cutting site.*)
    procedure PATTERNMATCHED;
      var NEWFRAG:FRAG;
          POSN:integer; (* # of nucleotide 5' to the cut on input strand *)
      begin
        with ENZ do begin
          POSN:= (Sk-PLEN)+CUTSITE+GOBACK[PLEN1];
          if CIRCULAR then (* Test for enz. which cuts beyond recog. site*)
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
            FOUND.LNUM:= FOUND.LNUM + 1
            end
          end
      end; (* PATTERNMATCHED *)
 
      begin  (* SEARCH *)
        WRITEWORD(output,ENZ.NAME,10);writeln;
        with ENZ do begin
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
 
      99: (* INPUTEXHAUSTED *) writeln
      end
    end; (* SEARCH *)
 
 
    begin  (* KMP *)
      FOUND.LNUM:=0;
      FOUND.TAIL^.START:=MAXSEQ+1;(*Always > POSN of site found*)
      with ENZ do begin
        (* Search SEQ using PATTERN derived from RSEQ  *)
        ENDNUC:= RSEQ[1];
        for Pj:= 1 to PLEN do PATTERN[Pj]:= NUCSET[RSEQ[Pj]];
        PATTERN[PLEN+1]:= [CANT]; (* Can not match anything *)
        RIGHTMOST:= FOUND.HEAD;
        MAKENEXT;
        if REBASE then SEARCH(CUT+PLEN)
        else SEARCH(CUT);
 
        (* Let PATTERN = the inverse complement of PATTERN*)
        (* Check for symmetry.                            *)
        Cj:= PLEN; 
        for Pj:= 1 to PLEN do begin
          CRP:= COM[RSEQ[Pj]];
          PATTERN[Cj]:= NUCSET[CRP];
          if RSEQ[Cj] <> CRP then SYMMETRIC:= false;
          Cj:= Cj-1
          end;
 
        if not SYMMETRIC then begin
          ENDNUC:= CRP;
          MAKENEXT;
          RIGHTMOST:= FOUND.HEAD;
          if REBASE then SEARCH(-(CUTOPP))
          else SEARCH(-(CUTOPP-PLEN))
          end
      end (* with ENZ *)
    end; (* KMP *)
 
 
    (*********************************************)
    (* Compile a report for output.              *)
    (*********************************************)
    (* Note (3/4/93). The original BACHREST had ORDER as a global array
       within procedure REPORT, which could be directly referenced by
       CALCULATE, BUBBLESORT and PRINTLIST. Although this worked well
       for the better part of a decade, code compiled by SUN Pascal with
       any level of optimization (-O, -O3, -O4) under SUN-OS 4.1.1 caused
       one or more elemets of the ORDER array to be reset to nil infrequently
       in  what appears to be a data-dependent fashion. It appears that ORDER
       was okay until PRINTLIST was called, at which point you would see
       an array element reset to nil. This didn't happen often, but with a
       large sequence run against a long list of enzymes, it was detected
       to happen in digests giving >30 fragments. No occurrences were seen
       with short fragment lists. It is also worth noting that this error 
       was not seen with INTREST when the same sequence and enzymes were
       used.

       Recompiling without optimization eliminated the problem. Therefore,
       this is probably a bug in Sun Pascal's optimization of pointer usage.
       As a work around, I have redefined ORDER as type SORTEDLIST to allow
       this array to be passed as a parameter. Probably what this does is to
       defeat some of the optimization that would normally be done. 
       Regardless of the cause, it works. *)
 
    procedure REPORT;
    type  SORTEDLIST = array[1..MAXFRAGS] of FRAG;
          
    var ORDER:SORTEDLIST;
        TOOMANY_FRAGS: boolean; (* =true if SITES > MAXFRAGS *)
        NUMSITES:integer;
 
    (*Calculate sizes and ends of fragments*)
    procedure CALCULATE(var ORDER:SORTEDLIST);
      var CURRENTFRAG:FRAG;
          SITE,I:integer;
      begin
        (* Initialize order array, to prevent optimizer errors *)
        for I:= 1 to MAXFRAGS do ORDER[I]:= nil;

        with FOUND do begin
          NUMSITES:= LNUM;
          (* if LINEAR add a fragment to head of list*)
          if (not CIRCULAR) and (HEAD^.NEXT^.START > 1) then begin
            GETFRAG(CURRENTFRAG);
            ADDFRAG(CURRENTFRAG,HEAD);
            CURRENTFRAG^.START:= 1;
            LNUM:= LNUM + 1
            end
          else CURRENTFRAG:= HEAD^.NEXT;
 
          if LNUM > 0 then begin
          (*Calculate ends and size of each fragment and assign it *)
          (* a place in the ORDER array, to be sorted later.    *)
            SITE:=1;
            while CURRENTFRAG^.NEXT <> TAIL do begin
              with CURRENTFRAG^ do begin
                FINISH:= NEXT^.START-1; SIZE:= FINISH-START+1
                end;
              if (SITE <= MAXFRAGS) then ORDER[SITE]:= CURRENTFRAG
              else TOOMANY_FRAGS:=true;               
              SITE:= SITE+1;
              CURRENTFRAG:= CURRENTFRAG^.NEXT
              end; (* while *)
            (*Last fragment in list is a special case*)
            with CURRENTFRAG^ do
              if CIRCULAR then begin
                FINISH:= HEAD^.NEXT^.START-1;
                SIZE:= (SEQLEN-START)+(FINISH+1);
                if FINISH=0 then FINISH:= SEQLEN
                end
              else begin
                   FINISH:= SEQLEN; SIZE:= FINISH-START+1
                   end;  
              if (SITE <= MAXFRAGS) then ORDER[SITE]:= CURRENTFRAG
              else TOOMANY_FRAGS:=true
            end; (* if LNUM > 0 *)
          end (* with FOUND *)
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
    procedure PRINTLIST(var ORDER:SORTEDLIST);
      var THISFRAG: FRAG;
          SITE:integer;
      begin
        with FOUND do begin
          writeln(OUTFILE,NUMSITES:6);
	  if NUMSITES <= FRAGPRINT then begin
             if NUMSITES = LNUM then THISFRAG:= HEAD
             else THISFRAG:= HEAD^.NEXT;
             SITE:=1;
             while THISFRAG^.NEXT <> TAIL do begin
               THISFRAG:= THISFRAG^.NEXT;
               write(OUTFILE,THISFRAG^.START:MAXPAT+33);
               with ORDER[SITE]^ do writeln(OUTFILE,SIZE:8,START:8,FINISH:8);
               SITE:= SITE+1
               end;
             (* Print last fragment. If no sites are found, a linear *)
             (* sequence will still yeild one fragment.              *)
             if NUMSITES < LNUM then with ORDER[SITE]^ do
                writeln(OUTFILE,SIZE:MAXPAT+41,START:8,FINISH:8)
	     end
	  else begin
	       writeln(OUTFILE,' ':MAXPAT+33,'(Fragments not shown)');
	       
	       (* This now looks unnecessary, but we'll keep it as a comment.	       
	       if TOOMANY_FRAGS then begin
		  writeln('>>> Number of fragments exceeds ',MAXFRAGS);
		  writeln('>>> Increase MAXFRAGS and recompile');
		  writeln(OUTFILE);
		  writeln(OUTFILE, '>>> Number of fragments exceeds ',MAXFRAGS);
		  writeln(OUTFILE, '>>> Increase MAXFRAGS and recompile');
	          end
		  *)
	       end	   
	end;
        writeln(OUTFILE)
      end; (* PRINTLIST *)
 
  begin (* REPORT *)
    TOOMANY_FRAGS:=false;
    CALCULATE(ORDER);
    if (FOUND.LNUM >= FRAGLEAST) and ((FOUND.LNUM <= FRAGMOST) or (FRAGMOST=MAXFRAGS) )then begin
        (* Print a header for each the enzyme *)
	with ENZ do begin
	     WRITEWORD(OUTFILE,NAME,15); write(OUTFILE,' ');
	     WRITEWORD(OUTFILE,RECSTR,MAXPAT);
	     if REBASE then begin 
                write(OUTFILE,' ');	  
	        end
	     else begin 
                  write(OUTFILE,CUT:3);
		  if SYMMETRIC then write(OUTFILE,' ':6)
		  else write(OUTFILE,' (',CUTOPP:3,')')	  
	        end
	     end;     
        if not TOOMANY_FRAGS then BUBBLESORT(FOUND.LNUM,1); 
        PRINTLIST(ORDER)
       end;
    with FOUND do RIDOF(HEAD,TAIL)
  end; (* REPORT *)
 
 
 
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

      (* Initialize horizontal output line *)
      with HLINE do begin
        for I:= 1 to MAXLINE do STR[I]:='_';
        LEN:=MAXLINE
        end;
  
      writeln('Enter filename of sequence to search:');
      GETFILE(INFILE,'I',IFN);
      NAME.LEN:=0;
      READSEQ(INFILE,SEQ,SEQLEN,NAME,CIRCULAR);
      if NAME.LEN=0 then begin
        writeln('Type name to appear on output:');
        INPWORD(NAME);readln
        end;     
      INITIALIZE;
      INITENDS;
 
      (* Open the restriction site file, skip the first two title lines, *)
      (*  and advance to first non-blank character.                      *)
      writeln('Enter restriction site filename:');
      GETFILE(RESTFILE,'I',RFN);
      CHECKFORMAT(RESTFILE,REBASE);

      (* Open an output file. *)
      writeln('Enter output filename:');
      GETFILE(OUTFILE,'O',OFN);
 
      (* MAIN MENU *)
      repeat
        writeln;
        WRITELINE(output,HLINE,80);writeln;
        writeln('BACHREST','MAIN MENU':30);
        WRITELINE(output,HLINE,80);writeln;
        write('Input file:        ');WRITELINE(output,IFN,IFN.LEN);writeln;
        write('Rest. enz. file:   ');WRITELINE(output,RFN,RFN.LEN);writeln;
        write('Output file:       ');WRITELINE(output,OFN,OFN.LEN);writeln;
        WRITELINE(output,HLINE,80);writeln;
        writeln(' ':20,'1) Read in a new sequence');
        writeln(' ':20,'2) Open a new rest. enz. file');
        writeln(' ':20,'3) Open a new output file');
        writeln(' ':20,'4) Set search parameters');
        writeln(' ':20,'5) Search for sites and write output to file');
        WRITELINE(output,HLINE,80);writeln;
        writeln('Type the number of your choice  (0 to quit program)');
        GETINTEGER(CHOICE,0,5);readln;

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
(*!!!         CLOSE(RESTFILE); *)
              writeln('Type Rest. Enz. filename:');
              GETFILE(RESTFILE,'I',RFN);
              CHECKFORMAT(RESTFILE,REBASE)
            end;
          3:begin
(*!!!         CLOSE(OUTFILE); *)
              writeln('Type output filename:');
              GETFILE(OUTFILE,'O',OFN);
            end;
          4:PARAMETERS;
          5:if OFN.LEN>0 then begin
        
               HEADER(OUTFILE);
               
               (* Reset file pointer to beginning of RESTFILE *)
               reset(RESTFILE);
               if REBASE then begin
                  J:=1;
                  while (RESTFILE^ <> '<') and (not eof(RESTFILE)) do begin
                        readln(RESTFILE);
                        J:=J+1
                        end;
                  if not (eof(RESTFILE)) then begin

                     (* The .lst files distributed with REBASE begin with header information
                        and an entry formula for each of six fields. If the header has
                        not been deleted, we need to read past it to the first true entry.*)
                     (* Set FIRSTITEM to # of line on which first data item begins *)
                     repeat
                       READLINE(RESTFILE,TEMPLINE);J:=J+1
                     until (TEMPLINE.STR[2]='1') and (TEMPLINE.STR[4]<> '<');
                     FIRSTITEM:=J-1;
                    (* Reset the file and advance to the first data item. *) 
                    reset(RESTFILE);
                    J:=1;
                    while J<FIRSTITEM do begin
                      readln(RESTFILE);
                      J:= J+1
                      end (* while *)
                    end  (* if not eof *)
                  end (* REBASE *)
               else begin
                    if not eof(RESTFILE) then readln(RESTFILE);
                    if not eof(RESTFILE) then readln(RESTFILE)
                  end;
                  
               (* Search loop *)
               while not eof(RESTFILE) do begin
                  FITSPARAMS:=true;
                  if REBASE then READREBASE(RESTFILE,ENZ,FITSPARAMS)
                  else READSITE(RESTFILE,ENZ);
		  
                  (* Check minimim and maximum lengths of cut sites
                  if ENZ.PLEN < MINSITE then FITSPARAMS:=false;
                  if ENZ.PLEN > MAXSITE then FITSPARAMS:=false;

		  (* Check for symmetry *)
		  TESTSYMMETRY(ENZ);
		  if (SYMM = 'S') and (not ENZ.SYMMETRIC) then FITSPARAMS:= false;
		  if (SYMM = 'A') and (ENZ.SYMMETRIC) then FITSPARAMS:= false;

		  (* Check for cutting position *)
		  TESTCUT(ENZ);
                  case ENZ.FRAGENDS of
		       '5':if PROT5 = 'N' then FITSPARAMS:= false;
		       '3':if PROT3 = 'N' then FITSPARAMS:= false;
		       'B':if BLUNT = 'N' then FITSPARAMS:= false
		       end;
		  
		  (* Run the search. *)  
                  if FITSPARAMS then begin
                     KMP(ENZ);
                     REPORT
                     end
                  end (* while not eof *)
               end 
            else begin writeln('>>> NO OUTPUT FILE CURRENTLY OPEN');
                       readln end
          end (* case *)
        until CHOICE = 0;

(*!!!*)CLOSE(RESTFILE); 
(*!!!*) CLOSE(OUTFILE) 
    end. (* BACHREST *)
 
