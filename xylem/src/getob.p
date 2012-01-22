  (* Modification history:
  14  Aug 2001 As of GenBank 126.0, the LOCUS line is scheduled to
               allow names of 18 char. and sequence lengths up to
               99,000,000,000. This rearrangement also affects the 
               position of the word 'circular'. GETOB has been updated
               to read both old and new formats.
     8  Oct 99 As of Dec. 1999 the NID field will be deleted from 
               the GenBank flatfile format. This code will still
               read the NID field from old files, but WRITEGB will
               only write an NID field if it exists.
     18 Mar 96 procedure RGB has been modified to be more forgiving of
               errors in GenBank entries that might otherwise cause
               GETOB to loop infinitely. Now, if the first column of
               a line doesn't match a legal field identifier (eg. L for 
               LOCUS) it will simply keep reading until if finds one.
               This fix is primarily aimed at GDE, which really messes up
               a number of fields, but may help with GenBank entries that
               have been filtered through other programs.
     27 Jan 96 NCBI_GI changed to NID.
     30 May 95 Changed CRUNCHOFFSET to ASCII 33 = "!". VT100 emulators
               didn't like ASCII 29. 
     21 Dec 94 Read past leading blanks in .ano file that were compressed
               by splitdb -c.
     13 Oct 94 Initialize NUMN. Failure to initialize causes n's to be written
               to expression file for non-translated features.               
     13 Oct 94 MAXSEQ increased to 750,000
     11 Sep 94 Added ability to handle single base_position as a location. 
      9 Sep 94 If codon_start=2,3, write n's to expression file, as is
               already done for outfile.
     13 Apr 94 Literature citation in .msg file
      2 Jan 94 add NCBI_GI field
     25 Oct 93 Fixed bug in COMPLEMENT. Only complemented 3' most location
               in an expression of the form complement(join(L1,L2..Ln))
      7 Jul 93 Added and deleted feature keys to comply with DDBJ/EMBL/GenBank
               FeatureTable:Definition Version 1.04
     18 Jun 92 bug fix: In external (indirect) references, nested expressions 
               were not properly copied to the output file.
     16 Jun 92 added feature keys 'contig' and 'chromosome'
     16 Jun 92 allow arbitrarily long NAMEFILE
     16 Jun 92 poly(<location>,x)
     16 May 92 For each feature evaluated, write a feature expression to
               EXPFILE, which upon evaluation, would return that feature.
               Also, a title field has been added to type OBJECT. This title
               is printed in all output files, replacing the old object
               identifier.  
     13 Sep 91 fixed MAKEFN; used to truncate last letter in name  
     12 Sep 91 -r automatically sets -c
                    FINDDATA automatically converts names or acc.#'s to
                    upper case before searching an index file. *)

  (***********************************************************)
  (*                                                         *)
  (*  GETOB     VERSION  1.2.9    8/14/2001  Standard Pascal *)
  (*            Brian Fristensky                             *)
  (*            Dept. of Plant Science                       *)
  (*            University of Manitoba                       *)
  (*            Winnipeg, MB R3T 2N2  CANADA                 *)
  (*                                                         *)
  (* SYNOPSIS                                                *)
  (* getob infile namefile anofile seqfile indfile message   *)
  (*       outfile [expfile]                                 *)
  (*                                                         *)
  (* DESCRIPTION                                             *)
  (*  Gets objects from GenBank file and writes to OUTFILE   *)
  (*  Conforms to "The DDBJ/EMBL/GenBank Feature Table:      *)
  (*     Definition" version 1.06 Feb. 1, 1994.              *)
  (*                                                         *)
  (*    -f        write each entry to a separate file.       *)
  (*              the filename consists of the locus name,   *)
  (*             followed by a 3-letter file extesion .obj   *)
  (*    -r        resolve indirect references from NAMEFILE  *)
  (*              into objects. Indirect references take the *)
  (*              form:                                      *)
  (*            @[<database>::]<accession>|<locus>:<location>*)
  (*              -r automatically sets -c.                  *)
  (*              -r and -e are mutually exclusive.          *)
  (*    -c        NAMEFILE contains ACCESSION numbers, rather*)
  (*              than locus names                           *)
  (*    -n        By default, the qualifier 'codon_start'    *)
  (*              is used to determine how many n's, if nec- *)
  (*              essary, must be added to the 5' end of CDS,*)
  (*              mat_peptide, or sig_peptide, to preserve   *)
  (*              the reading frame. To turn OFF this feature*)
  (*              -n must be set. -n must be set for GenBank *)
  (*              Release 67.0 or earlier.                   *)
  (*                                                         *)
  (*    INFILE  - instructions for what data to pull         *)
  (*    NAMEFILE- contains names of entries to get           *)
  (*    ANOFILE - contains annotation parts of entries       *)
  (*    SEQFILE - contains sequence parts of entries         *)
  (*    INDFILE - contains linenumbers for beginning of each *)
  (*              entry in ANOFILE and SEQFILE               *)
  (*    MESSAGE - lists objects processed in each entry      *)
  (*    OUTFILE - objects retrieved, if not -f               *)
  (*    EXPFILE - contains feature expressions evaluated     *)
  (*                                                         *)
  (*  Copyright (c) 1989-1996  by Brian Fristensky           *)
  (*  !!! in comment indicates feature which may need change *)
  (***********************************************************)
  program GETOB(INFILE,NAMEFILE,ANOFILE,SEQFILE,INDFILE,MESSAGE,OUTFILE,EXPFILE);
(*!!!  Some Pascals require file parameters in program heading *)

  (******************************************************************)
  (*                     const  definition                          *)
  (******************************************************************)

  const MAXWORD = 35;     (*Max. length of WORD.STR *)
        MAXFK   = 15;     (*Length of FK, Feature Key string *)
        MAXLINE = 132;    (*Max. length of GenBank entry line *)
        MAXLLEN = 1000;   (*Max. chars. in a location *)
        MAXSEQ  = 750000; (*Max. sequence in an entry *)
        MAXOBJ  = 1000;    (*Max. objects in an entry *)
        CRUNCHOFFSET = 33; (* ASCII char. used by splitdb for -c option. *)
        DEBUG   = false;  (*=true, debugging messages are printed to output *)

(* BEGIN MODULE STARTARGNUM *)
	STARTARGNUM=1;    (* SUN Pascal: ARG(1) is 1st command line argument*)
      (*STARTARGNUM=2;       ATT Pascal: ARG(2) is 1st command line argument*)
(* END MODULE STARTARGNUM         VERSION= 'SUNMODS     Version  8/ 9/94'; *)

        VERSION = 'GETOB          Version 1.2.9  14 Aug 2001';

  (******************************************************************)
  (*                     type   definition                          *)
  (******************************************************************)

  type        
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

       (* -  - - - Locations of entries in database files - - - - - - *)
       LOC = record
             NAME:WORD;
             ANOLINE,
             SEQLINE:integer 
             end;
    
      (* - - - - - - - - -  NUCLEOTIDE types - - - - - - - - -  *)
      NUCLEOTIDE  = (T,C,A,G,R,Y,M,W,S,K,D,H,V,B,N);
      SEQUENCE = array[1..MAXSEQ] of NUCLEOTIDE;
      CAR      = array[NUCLEOTIDE] of char;
      NAR      = array[NUCLEOTIDE] of NUCLEOTIDE;

      NP    = ^NODE;


      (* - - - - - - - - - - - Fragment types - - - - - - - - - *)
      (*- used as building blocks to create a model of an object. *)

      LOCLINE = array[1..MAXLLEN] of char; (* string that can be parsed
                                              into a location *)

      LOCTYPE  = (position,baserange,betweenposition,featurename);

      BASEPOSITION = record
        LNUM,RNUM:integer; (* if LNUM<>RNUM, this is a two_base_bound *)
        BOUNDFLAG:char
        end; (* BASEPOSITION *)
        
      FRAGMENT = record
        REPEATS:integer; (* used by POLY(), #times to repeat an expr. *)
        case LT:LOCTYPE of
          position,baserange,betweenposition:
                       (START,FINISH:BASEPOSITION;STRAND:char;SIZE:integer);
          featurename:(EXLOC:LOCLINE;EXLOCLEN:integer)
        end; (* FRAGMENT *)   
           
      LIST = record
                HEAD,TAIL:NP
             end;

      (* - - - - - - - - - - - FEATURES table. - - - - - - - - - -*)

      (* KEY IDentifiers for FEATURES table *)
      FEATUREKEY = (BLANK,
           ALLELE,ATTENUATOR,
           BINDING,
           CAATSIGNAL,CDS,CHROMOSOME, CONFLICT, CONTIG,CREGION,
           DLOOP,DREGION,DSEGMENT,
           ENHANCER,EXON,
           GCSIGNAL,
           HYPHEN,
           IDNA,INTRON,
           JREGION,JSEGMENT,
           LTR,
           MATPEPTIDE,MISCBINDING,MISCDIFFERENCE,
           MISCFEATURE, MISCRECOMB, MISCRNA,MISCSIGNAL,MISCSTRUCTURE,
           MODIFIEDBASE,MRNA, MUTATION,
           NREGION,
           OLDSEQUENCE,
           POLYASIGNAL,POLYASITE,
           PRECURSORRNA, PRIMTRANSCRIPT, PRIMERBIND,PROMOTER,PROTEINBIND,
           RBS, REPEATREGION, REPEATUNIT,REPORIGIN,RRNA,
           SATELLITE,SCRNA,SIGPEPTIDE, SNRNA,SOURCE,STEMLOOP,STS,SREGION,
           TATASIGNAL,TERMINATOR,TRANSITPEPTIDE,TRNA,
           UNSURE,
           VARIATION,VIRION,VREGION,VSEGMENT,
           TPCLIP,TPUTR,
           FPCLIP,
           FPUTR,
           M10SIGNAL,M35SIGNAL,
           OTHER);

      KEYSET = set of FEATUREKEY;
             
         FK      = packed array[1..MAXFK ] of char; (* FEATURE key string *)
         FEATYPE = record
                 KEY:FEATUREKEY;
                 KEYSTRING:FK; 
                 LOCQUAL:LINE (* LOCATION OR QUALIFIER *)
                 end;

         (* - - - - - - - - Linked Lists - - - - - - - - - - - - *)
         (* General purpose linked list structure for GenBank entries and
            also for sequence models. *)
         NODETYPES = (WNODE,LNODE,FEANODE,FRANODE);
         NODE  = record
                 PREV,NEXT:NP;
                 case NTYPE:NODETYPES of
                   WNODE:(WOR:WORD);
                   LNODE:(LIN:LINE);
                   FEANODE:(FEA:FEATYPE);
                   FRANODE:(FRA:FRAGMENT)
                 end;

         (* - - - - - - - - - - GenBank Entry - - - - - - - - - -  *)
         ENTRY = record
                 NAME:WORD;
                 STARTNUM,SEQLEN:integer;
                 STRANDEDNESS:(ss,ds,ms);
                 SEQTYPE:WORD;
                 CIRCULAR:boolean;
                 ENTRYSTATUS:WORD;
                 ENTRYDATE:WORD;
                 DEFINITION:LIST;
                 ACCESSION:LIST;
                 NID:WORD;
                 SEGNUMBER,SEGTOTAL:integer;
                 FEATURETABLE:LIST;
                 ACOMP,CCOMP,GCOMP,TCOMP:integer
                 end; (* ENTRY *)

      (* - - - - - - - - - - - - OBJECTS - - - - - - - - - - - - - *)
      (* Model of sequence, consisting of a linked list of fragments which 
         together, make up a feature. Up to MAXOBJ objects can be evaluated
         for a given entry. *)
      OBJECT = record
                 TITLE:WORD;
                 NOLABEL:boolean; (*=true if feature didn't have a label *)
                 FEALAB:WORD;
                 HEAD,TAIL:NP;
                 end;
      OBJECTS = array[0..MAXOBJ] of OBJECT;

  (******************************************************************)
  (*                       var  definition                          *)
  (******************************************************************)

  var 
      (* File variables *)
      INFILE,NAMEFILE,ANOFILE,SEQFILE,INDFILE,MESSAGE,OUTFILE,EXPFILE:text;
      FILENAME:LINE;
      ARGNUM:integer;

      (* Option variables *)
      FILES,        (* -f: =true, each entry written to separate file *)
      RESOLVEXP,    (* -r: =true, references in NAMEFILE resolved into objects*)
      ACNO,         (* -c: =true, NAMEFILE contains ACCESSION numbers *)
      PAD5PRIME:boolean;    (* -n: =true, if -n is not set *)

      (* Variables associated with NAMEFILE and INDFILE *)
      SEQNAME:WORD;
      INDIRREF:LINE; (* Holds indirect reference from NAMEFILE *)
      LOCATIONS:LOC; (* locations of a locus in .ANO and .SEQ files *)
      FOUND,         (* =true if entry has been found in INDFILE *)
      OKAY:boolean;  (* =true if INDIRREF parsed into legal <location> *)
 
      (* Global type conversion arrays. *)
      NUCHAR:CAR;
      INUC,COMP:NAR;
      KEYSTR:array[FEATUREKEY] of FK;

      (* Global data structures *)
      E:ENTRY;
      SEQ:SEQUENCE;
      SEQLIM, (* SEQ[SEQLIM+1..MAXSEQ] is unused as of yet *)
      NUMOBJ:integer; (* number of objects in OBJ *)
      OBJ:OBJECTS;    (* set of objects found in an entry *)
      OBJECTTYPE:FK;
      ACCESSION:WORD;
      LL:LOCLINE;
      LLEN:integer;
      TRANSCRIPTS,SITESET:KEYSET; (*tell which feature keys are to be
                                     searched for in each entry *)
      FREENODE:NP; (* beginning of free node list *)

      (* Miscellaneous variables for indexing, reading etc. *) 
      CIRC:WORD; (* holds 'circular' as a WORD *)
      CA,CS,     (* current lines in annotation and sequence files*)
      I:integer;
      CH,CRUNCHFLAG:char;

  (******************************************************************)
  (*              procedure and function  definition                *)
  (******************************************************************)

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
(* END MODULE READLINE         VERSION= 'SUNMODS     Version  8/ 9/94'; *)

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
(* END MODULE WRITELINE         VERSION= 'SUNMODS     Version  8/ 9/94'; *)

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
(* END MODULE SAMEWORD         VERSION= 'SUNMODS     Version  8/ 9/94'; *)

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
(* END MODULE TOUPPER         VERSION= 'SUNMODS     Version  8/ 9/94'; *)

  (*****************************************************)
  (* Keyword procedures.                               *)
  (*****************************************************)
  procedure READFK(var F:text; var KEYWORD:FK);
    var I:integer;
    begin
      for I:= 1 to MAXFK do
        if not eoln(F) then read(F,KEYWORD[I])
        else KEYWORD[I]:=' '
    end; (* READFK *)

  (* Convert an input string to a keyword from FEATURETABLE.  If
     string can't be identified, assign it to blank and write message
     to MESSAGE file. *)
  procedure STR2FK(KEYSTRING:FK; var NEWKEY:FEATUREKEY);
    begin
      NEWKEY:=BLANK;
      while (KEYSTRING <> KEYSTR[NEWKEY]) and (NEWKEY < OTHER) do
        NEWKEY:= succ(NEWKEY)
    end; (* STR2FK *)
   
  (*****************************************************)
  (* Initialization procedures.                        *)
  (*****************************************************)

  (* Initialize the variables.                         *)
  procedure INITIALIZE;
    var I:integer;

    (* Create the dummy nodes HEAD & TAIL, which delimit a list *)
    procedure INITLIST(var HEAD,TAIL:NP);
      begin
        new(HEAD);new(TAIL);
        with HEAD^ do begin NEXT:=TAIL;PREV:=nil end;
        with TAIL^ do begin NEXT:=nil;PREV:=HEAD end
      end; (* INITLIST *)

    begin
     (*   Initialize INUC, NUCHAR and COMP, the arrays which *)
     (*   hold the  values of the NUCLEOTIDEs.      *)
        INUC[T]:=T;NUCHAR[T]:='t';COMP[T]:=A;
        INUC[C]:=C;NUCHAR[C]:='c';COMP[C]:=G;
        INUC[A]:=A;NUCHAR[A]:='a';COMP[A]:=T;
        INUC[G]:=G;NUCHAR[G]:='g';COMP[G]:=C;
        INUC[R]:=R;NUCHAR[R]:='r';COMP[R]:=Y;
        INUC[D]:=D;NUCHAR[D]:='d';COMP[D]:=H;
        INUC[V]:=V;NUCHAR[V]:='v';COMP[V]:=B;
        INUC[M]:=M;NUCHAR[M]:='m';COMP[M]:=K;
        INUC[K]:=K;NUCHAR[K]:='k';COMP[K]:=M;
        INUC[B]:=B;NUCHAR[B]:='b';COMP[B]:=V;
        INUC[H]:=H;NUCHAR[H]:='h';COMP[H]:=D;
        INUC[Y]:=Y;NUCHAR[Y]:='y';COMP[Y]:=R;
        INUC[W]:=W;NUCHAR[W]:='w';COMP[W]:=W;
        INUC[S]:=S;NUCHAR[S]:='s';COMP[S]:=S;
        INUC[N]:=N;NUCHAR[N]:='n';COMP[N]:=N;

    (* Initialize keyword array *)
    KEYSTR[BLANK          ]:='               ';
    KEYSTR[ALLELE         ]:='allele         ';
    KEYSTR[ATTENUATOR     ]:='attenuator     ';
    KEYSTR[BINDING        ]:='binding        ';
    KEYSTR[CAATSIGNAL     ]:='CAAT_signal    ';
    KEYSTR[CDS            ]:='CDS            ';
    KEYSTR[CHROMOSOME     ]:='chromosome     ';
    KEYSTR[CONFLICT       ]:='conflict       ';
    KEYSTR[CONTIG         ]:='contig         ';
    KEYSTR[CREGION        ]:='C_region       ';
    KEYSTR[DLOOP          ]:='D-loop         ';
    KEYSTR[DREGION        ]:='D_region       ';
    KEYSTR[DSEGMENT       ]:='D_segment      ';
    KEYSTR[ENHANCER       ]:='enhancer       ';
    KEYSTR[EXON           ]:='exon           ';
    KEYSTR[GCSIGNAL       ]:='GC_signal      ';
    KEYSTR[IDNA           ]:='iDNA           ';
    KEYSTR[INTRON         ]:='intron         ';
    KEYSTR[JREGION        ]:='J_region       ';
    KEYSTR[JSEGMENT       ]:='J_segment      ';
    KEYSTR[LTR            ]:='LTR            ';
    KEYSTR[MATPEPTIDE     ]:='mat_peptide    ';
    KEYSTR[MISCBINDING    ]:='misc_binding   ';
    KEYSTR[MISCDIFFERENCE ]:='misc_difference';
    KEYSTR[MISCFEATURE    ]:='misc_feature   ';
    KEYSTR[MISCRECOMB     ]:='misc_recomb    ';
    KEYSTR[MISCRNA        ]:='misc_RNA       ';
    KEYSTR[MISCSIGNAL     ]:='misc_signal    ';
    KEYSTR[MISCSTRUCTURE  ]:='misc_structure ';
    KEYSTR[MODIFIEDBASE   ]:='modified_base  ';
    KEYSTR[MRNA           ]:='mRNA           ';
    KEYSTR[MUTATION       ]:='mutation       ';
    KEYSTR[NREGION        ]:='N_region       ';
    KEYSTR[OLDSEQUENCE    ]:='old_sequence   ';
    KEYSTR[POLYASIGNAL    ]:='polyA_signal   ';
    KEYSTR[POLYASITE      ]:='polyA_site     ';
    KEYSTR[PRECURSORRNA   ]:='precursor_RNA  ';
    KEYSTR[PRIMERBIND     ]:='primer_bind    ';
    KEYSTR[PRIMTRANSCRIPT ]:='prim_transcript';
    KEYSTR[PROMOTER       ]:='promoter       ';
    KEYSTR[PROTEINBIND    ]:='protein_bind   ';
    KEYSTR[RBS            ]:='RBS            ';
    KEYSTR[REPEATREGION   ]:='repeat_region  ';
    KEYSTR[REPEATUNIT     ]:='repeat_unit    ';
    KEYSTR[REPORIGIN      ]:='rep_origin     '; 
    KEYSTR[RRNA           ]:='rRNA           ';
    KEYSTR[SATELLITE      ]:='satellite      ';
    KEYSTR[SCRNA          ]:='scRNA          ';
    KEYSTR[SIGPEPTIDE     ]:='sig_peptide    ';
    KEYSTR[SNRNA          ]:='snRNA          ';
    KEYSTR[SOURCE         ]:='source         ';
    KEYSTR[SREGION        ]:='S_region       ';
    KEYSTR[STEMLOOP       ]:='stem_loop      ';
    KEYSTR[STS            ]:='STS            ';
    KEYSTR[TATASIGNAL     ]:='TATA_signal    ';
    KEYSTR[TERMINATOR     ]:='terminator     ';
    KEYSTR[TRANSITPEPTIDE ]:='transit_peptide';
    KEYSTR[TRNA           ]:='tRNA           ';
    KEYSTR[UNSURE         ]:='unsure         ';
    KEYSTR[VARIATION      ]:='variation      ';
    KEYSTR[VIRION         ]:='virion         ';
    KEYSTR[VREGION        ]:='V_region       ';
    KEYSTR[VSEGMENT       ]:='V_segment      ';
    KEYSTR[TPCLIP         ]:='3''clip         ';
    KEYSTR[TPUTR          ]:='3''UTR          ';
    KEYSTR[FPCLIP         ]:='5''clip         ';
    KEYSTR[FPUTR          ]:='5''UTR          ';
    KEYSTR[M10SIGNAL      ]:='-10_signal     ';
    KEYSTR[M35SIGNAL      ]:='-35_signal     ';
    KEYSTR[HYPHEN         ]:='-              '; 
    KEYSTR[OTHER          ]:='               '; 

    with E do begin
      CRUNCHFLAG:=chr(CRUNCHOFFSET); (*used in RGB *)
      NAME.LEN:=0; SEQLEN:=0; SEQTYPE.LEN:=0; ENTRYSTATUS.LEN:=0;
      ENTRYDATE.LEN:=0;
      INITLIST(DEFINITION.HEAD,DEFINITION.TAIL);
      INITLIST(ACCESSION.HEAD,ACCESSION.TAIL);
      SEGNUMBER:=0; SEGTOTAL:=0; NID.LEN:=0;
      INITLIST(FEATURETABLE.HEAD,FEATURETABLE.TAIL);
      with FEATURETABLE.TAIL^ do begin
        NTYPE:=FEANODE; FEA.KEY:=BLANK
        end; (* BLANK indicates end of features table *)
      ACOMP:=0; CCOMP:=0; GCOMP:=0; TCOMP:=0;
      for I:= 0 to MAXOBJ do INITLIST(OBJ[I].HEAD,OBJ[I].TAIL);
      FREENODE:=nil;
      
      (* WORDS used to parse input *)
      with CIRC do begin
           STR[1]:='c';STR[2]:='i';STR[3]:='r';STR[4]:='c';STR[5]:='u';
           STR[6]:='l';STR[7]:='a';STR[8]:='r';LEN:=8
        end
        
      end (* E *)
  end; (* INITIALIZE *)
 
    (* - - - - - - - - - - - - - - - - - - - - - - - - - - - *)
    (* Read options from command line and set FILES & RESOLVEXP. *)
    procedure READOPTIONS(var ARGNUM:integer; var FILES,RESOLVEXP,ACNO:boolean);
      var ARGUMENT:packed array[1..132] of char;
     begin (* READOPTIONS *)
       (* Set defaults *)
       FILES:=false; RESOLVEXP:=false; ACNO:=false; PAD5PRIME:=true;
       (* Read options *)
       argv(ARGNUM,ARGUMENT);
       while ARGUMENT[1]='-' do begin
         if ARGUMENT[2] in ['f','r','c','n'] then 
           case ARGUMENT[2] of 
            'f':FILES:=true;
                (* indirect references always contain accession numbers *)
            'r':begin RESOLVEXP:=true; ACNO:= true end;
            'c':ACNO:=true;
            'n':PAD5PRIME:=false
             end;
         ARGNUM:=ARGNUM+1;
         argv(ARGNUM,ARGUMENT)
         end
       end; (* READOPTIONS *)      
 
   (* - - - - - - - - - - - - - - - - - - - - - - - - - - - *)
   (* Open files.*)
   procedure OPENFILES(var ARGNUM:integer; var INFILE,NAMEFILE,ANOFILE,SEQFILE,
			   INDFILE,MESSAGE, OUTFILE, EXPFILE:text);

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
(* END MODULE FILEARGS         VERSION= 'SUNMODS     Version  8/ 9/94'; *)

     begin
       FILEARGS(INFILE,'I',ARGNUM);
       FILEARGS(NAMEFILE,'I',ARGNUM);
       FILEARGS(ANOFILE,'I',ARGNUM);
       FILEARGS(SEQFILE,'I',ARGNUM);
       FILEARGS(INDFILE,'I',ARGNUM);
       FILEARGS(MESSAGE,'O',ARGNUM);
       if not FILES then FILEARGS(OUTFILE,'O',ARGNUM);
       if not RESOLVEXP then FILEARGS(EXPFILE,'O',ARGNUM)
     end; (* OPENFILES *)

   (* - - - - - - - - - - - - - - - - - - - - - - - - - - - *)
   (* Read input file containing instructions for processing each entry.*)
    procedure READINSTRUCTIONS(var INFILE:text; 
                               var TRANSCRIPTS,SITESET:KEYSET;
                               var OBJECTTYPE:FK);
      var OBJID:FK;
          OBJKEY:FEATUREKEY;
      begin
        TRANSCRIPTS:=[]; SITESET:=[];
        READFK(INFILE,OBJECTTYPE);readln(INFILE);
        if (OBJECTTYPE<>'GENBANK        ') then begin
          READFK(INFILE,OBJID);readln(INFILE);
          while (OBJID<>'SITES          ') and (not eof(INFILE)) do begin 
            STR2FK(OBJID,OBJKEY);
            if OBJKEY<>OTHER then TRANSCRIPTS:= TRANSCRIPTS + [OBJKEY];
            if not eof(INFILE) then begin
               READFK(INFILE,OBJID);readln(INFILE) end
            end; (* OBJID <> 'SITES         ' *)
          while not eof(INFILE) do begin
            STR2FK(OBJID,OBJKEY);
            if OBJKEY<>OTHER then SITESET:= SITESET + [OBJKEY];
            if not eof(INFILE) then begin
               READFK(INFILE,OBJID);readln(INFILE) end
            end (* not eof(INFILE) *)
(*!!!*)   end; (* OBJECTTYPE<> 'GENBANK        ' *)
(*!!!   CLOSE(INFILE); *)
      end; (* READINSTRUCTIONS *)
 
  (**************************************************************)
  (* Read a LOCUS name for next sequence to pull from database. *)
  (**************************************************************)
  procedure READNAME(var NAMEFILE:text; var SEQNAME:WORD);
    begin
      (* Skip comment lines *)
      while (NAMEFILE^=';') and (not eof(NAMEFILE)) do readln(NAMEFILE);

      (* Read the first token in the next non-comment line *)
      SEQNAME.LEN:=0;
      if not eof(NAMEFILE) then READWORD(NAMEFILE,SEQNAME);
      if not eof(NAMEFILE) then readln(NAMEFILE)
    end; (* READNAME *)

  (**********************************************************************)
  (* Look for a given LOCUS in the index file.  Since names in NAMEFILE *)
  (* are supposed to be sorted, reaching the end of file means that the *)
  (* locus in question is not in the database.                          *)
  (**********************************************************************)
  procedure FINDDATA(var INDFILE:text; ID:WORD; var LOCATIONS:LOC;
                     var FOUND:boolean);
    var FLAG:char;
        DUMMY:WORD;

   (* Advance to first line with first char of ID *)
   procedure ZOOM;
     begin
       while (not eof(INDFILE)) and (INDFILE^<>FLAG) do readln(INDFILE)
     end; (* ZOOM *)

   (* Search for LOCUS corresponding to ID *) 
   procedure LCRAWL;
     begin
       while (INDFILE^=FLAG) and (not FOUND) do
         with LOCATIONS do begin
           READWORD(INDFILE,NAME);
           if SAMEWORD(ID,NAME) then begin
             FOUND:=true;
             READWORD(INDFILE,DUMMY); (* ignore ACCESSION number *)
             read(INDFILE,ANOLINE,SEQLINE)
             end;
           if not eof(INDFILE) then readln(INDFILE)
           end (* with *)
      end; (* LCRAWL *)

   (* Search for ACCESSION number corresponding to ID *) 
   procedure ACRAWL;
     begin
       while (not FOUND) and (not eof(INDFILE)) do
         with LOCATIONS do begin
           READWORD(INDFILE,DUMMY); (* skip over LOCUS name *)
           READWORD(INDFILE,NAME);
           if SAMEWORD(ID,NAME) then begin
             FOUND:=true;
             read(INDFILE,ANOLINE,SEQLINE)
             end;
           if not eof(INDFILE) then readln(INDFILE)
           end (* with *)
      end; (* ACRAWL *)

   begin (* FINDDATA *)
     if eof(INDFILE) then reset(INDFILE);
     FOUND:=false;
     (* Make sure that the name or accession number is capitalized. *)
     with ID do for I:= 1 to LEN do STR[I]:= TOUPPER(STR[I]);
     FLAG:=ID.STR[1];

     (* Search to end of file or until FOUND *)
     if ACNO then ACRAWL
     else while (not FOUND) and (not eof(INDFILE)) do begin ZOOM; LCRAWL end;

     (* If list isn't sorted, go back to the beginning and search from top *)
     if not FOUND then begin 
       reset(INDFILE);
       if ACNO then ACRAWL
       else while (not FOUND) and (not eof(INDFILE)) do begin ZOOM; LCRAWL end;
       end
  end; (* FINDDATA *)

  (******************************************)
  (* Parse a reference into its components. *)
  (******************************************)
  procedure PARSEREF(var INDIRREF:LINE; var ACCESSION:WORD; 
                     var LL:LOCLINE; var LLEN:integer; var OKAY:boolean);
    var I,J,STARTACC,FINISHACC:integer;
    begin (* PARSREF *)
      with INDIRREF do begin
        (* Ignore DATABASE component. It is assumed that the necessary entries
          have already been obtained from the appropriate database and 
          included in ANOFILE, SEQFILE, INDFILE *)
        I:=2; OKAY:=false;
        while (I<LEN) and (STR[I]<>':') do I:=I+1;
        if I<LEN then begin (* 1*)
          if STR[I+1]=':' then begin
             STARTACC:=I+2;
             I:=STARTACC+1;
             while (I<LEN) and (STR[I]<>':') do I:=I+1
             end
          else STARTACC:=2;
          FINISHACC:=I-1;

          (* Read the accession number *)
          if I<LEN then begin (* 2 *)
            J:=1;
            for I:= STARTACC to FINISHACC do begin
                ACCESSION.STR[J]:=STR[I]; J:=J+1 end;
            ACCESSION.LEN:= FINISHACC-STARTACC+1;
  
            (* Read the <location> expression *)
            J:=1;
            for I:= FINISHACC+2 to LEN do begin
                LL[J]:=STR[I]; J:=J+1 end;
            LLEN:= LEN-FINISHACC-1;
            OKAY:=true
            end (* 2 *)
          end (* 1 *)
        end (* with INDIRREF *)
    end; (* PARSREF *)

  (*****************************************************)
  (* Create a filename using the locus.                *)
  (*****************************************************)
  procedure MAKEFN(NAME:WORD; var FILENAME:LINE);
    var I:integer;
    begin
      with FILENAME do begin
        I:=1;
        while I <= NAME.LEN do begin
          STR[I]:= NAME.STR[I];
          I:=I+1
          end;
        STR[I]:='.';
        STR[I+1]:='o'; STR[I+2]:='b'; STR[I+3]:='j'; 
        LEN:=I+3;
        for I:= LEN+1 to MAXLINE do STR[I]:=' '
        end
    end; (* MAKEFN *)

  (*****************************************)
  (*  Linked-list operations for node list.*)
  (*****************************************)
 
  (*Add a node after DEST*)
  procedure ADDNODE(DEST:NP);
    var TEMP,NEWNODE:NP;
    begin
      (*Get a new node from freelist.*)
      if FREENODE = nil then new(NEWNODE)
      else begin
        NEWNODE:= FREENODE;
        FREENODE:= FREENODE^.NEXT
        end;
      (* Add the node *)
      TEMP:=DEST^.NEXT;
      NEWNODE^.NEXT:= TEMP;
      NEWNODE^.PREV:=DEST;
      DEST^.NEXT:= NEWNODE;
      TEMP^.PREV:= NEWNODE
    end; (* ADDNODE *)
 
  (*Return a list to the top of freelist*)
  procedure RIDOF(var HEAD,TAIL:NP);
    var FIRST,LAST:NP;
    begin
      if HEAD^.NEXT <> TAIL then begin
        FIRST:=HEAD^.NEXT;LAST:=TAIL^.PREV;
        FIRST^.PREV:= nil;LAST^.NEXT:= FREENODE;
        FREENODE:= FIRST;
        HEAD^.NEXT:=TAIL;TAIL^.PREV:=HEAD
        end
    end; (* RIDOF *)

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
(* END MODULE NUC         VERSION= 'SUNMODS     Version  8/ 9/94'; *)
 
    (*******************************************************)
    (* Read a GenBank entry into a structure.              *)
    (*******************************************************)
    procedure RGB(var ANOFILE,SEQFILE:text; LOCATIONS:LOC; var E:ENTRY);
      var DONE:boolean;
          BLANKSET:set of char;
          CH:char;

      procedure SKIP(var F:text; POS:integer);
        var I:integer;
            CH:char;
        begin
          I:= 0;
          (* Read past compressed leading blanks, if any.*)
          if (F^=CRUNCHFLAG) then begin
             read(F,CH);read(F,CH) end
          (* Otherwise, read past the specified number of characters.*)
          else
            while (I < POS) and (not eoln(F)) do begin
              read(F,CH);
              I:=I+1
              end               
        end; (* SKIP *)

      (* - - - - - - - - - - - - - - - - - - - - - - - - - *)
      procedure READLOCUS;
        var TOKEN:WORD;
        begin
          with E do begin
            SKIP(ANOFILE,12);
            READWORD(ANOFILE,NAME);
            read(ANOFILE,SEQLEN);
            READWORD(ANOFILE,TOKEN); (* bp *)
            (* After 'bp', there is an optional field telling
            the type of molecule (ss-RNA, ds-DNA etc.)
            Since this field is optional, we must test the
            next two tokens to see if they are 'circular' *)
            CIRCULAR:=false;
            if not eof(ANOFILE) then READWORD(ANOFILE,TOKEN);
            if SAMEWORD(TOKEN,CIRC) then CIRCULAR:=true
            else begin
                 READWORD(ANOFILE,TOKEN);
                 if SAMEWORD(TOKEN,CIRC) then CIRCULAR:=true
                 end;
            readln(ANOFILE);
(*!!!*)     if DEBUG then begin WRITEWORD(output,NAME,18);writeln(SEQLEN) end
            end; (* with E *)
          CA:= CA + 1          
        end; (* READLOCUS *)

      (* - - - - - - - - - - - - - - - - - - - - - - - - - *)
      procedure READDEFINITION;
        begin
          with E.DEFINITION do
            repeat 
              SKIP(ANOFILE,12);
              ADDNODE(TAIL^.PREV);
              with TAIL^.PREV^ do begin
                NTYPE:=LNODE;
                READLINE(ANOFILE,LIN);
                CA:= CA + 1;
(*!!!*) if DEBUG then begin WRITELINE(output,LIN,LIN.LEN);writeln end
                end (* with TAIL^ *)
            until not(ANOFILE^ in BLANKSET)
        end; (* READDEFINITION *)

      (* - - - - - - - - - - - - - - - - - - - - - - - - - *)
      procedure READACCESSION;
        begin
          with E.ACCESSION do 
            repeat
              SKIP(ANOFILE,12);
              while not eoln(ANOFILE) do begin
                ADDNODE(TAIL^.PREV);
                TAIL^.PREV^.NTYPE:=WNODE;
                READWORD(ANOFILE,TAIL^.PREV^.WOR);
(*!!!*) if DEBUG then begin WRITEWORD(output,TAIL^.PREV^.WOR,10) end
                end;      
             readln(ANOFILE);
(*!!!*) if DEBUG then writeln;
             CA:= CA + 1
           until not(ANOFILE^ in BLANKSET)
        end; (* READACCESSION *)

      (* - - - - - - - - - - - - - - - - - - - - - - - - - *)
      procedure READNID;
        begin
(*!!!*)if DEBUG then writeln(output,'NID');
          SKIP(ANOFILE,11);
          READWORD(ANOFILE,E.NID);readln(ANOFILE);
          CA:= CA+1
        end; (* READNID *)

      (* - - - - - - - - - - - - - - - - - - - - - - - - - *)
      procedure READSEGMENT;
        var DUMMY:WORD;
        begin
(*!!!*)if DEBUG then writeln(output,'SEGMENT');
          SKIP(ANOFILE,11);
          read(ANOFILE,E.SEGNUMBER);
          READWORD(ANOFILE,DUMMY);
          readln(ANOFILE,E.SEGTOTAL);
          CA:= CA+1
        end; (* READSEGMENT *)

      (* - - - - - - - - - - - - - - - - - - - - - - - - - *)
      (* This procedure skips over parts of an entry that
         aren't used by GETOB. *)
      procedure READDUMMY;
        begin
          repeat
            readln(ANOFILE);
            CA:=CA+1
          until not(ANOFILE^ in BLANKSET)
        end; (* READDUMMY *)

      (* - - - - - - - - - - - - - - - - - - - - - - - - - *)
      procedure READFEATURES;

        procedure FEALINE(var DEST:NP);
          var CH:char;  
              NUMBLANKS:integer;     
          begin
              ADDNODE(DEST);
              DEST^.NTYPE:=FEANODE;
              with DEST^.FEA do begin
                if ANOFILE^=CRUNCHFLAG then begin
                   read(ANOFILE,CH); read(ANOFILE,CH);
                   NUMBLANKS:=ord(CH)-CRUNCHOFFSET;
                   if NUMBLANKS=5 then begin (*next field is feature key *)
                      READFK(ANOFILE,KEYSTRING);
                      STR2FK(KEYSTRING,KEY);
                      read(ANOFILE,CH)
                      end
                   else KEY:=BLANK (*next field is qualifier line *)
                   end (* ANOFILE^=CRUNCHFLAG *)
                else begin
                  SKIP(ANOFILE,5);
                  READFK(ANOFILE,KEYSTRING);
                  STR2FK(KEYSTRING,KEY);
                  read(ANOFILE,CH)
                  end;
                READLINE(ANOFILE,LOCQUAL);
                CA:=CA+1;
(*!!!*) if DEBUG then begin
                 write(output,KEYSTR[KEY],'   '); 
                 WRITELINE(output,LOCQUAL,LOCQUAL.LEN);
                 writeln(output)
                 end          
                end (* with DEST^ *)    
          end; (* FEALINE *)

        begin (* READFEATURES *)
(*!!!*)if DEBUG then writeln(output,'FEATURES');
          readln(ANOFILE); CA:=CA+1;

          (* Read the FEATURES table *)
          with E.FEATURETABLE do
            while (ANOFILE^ in BLANKSET) do FEALINE(TAIL^.PREV)
  
        end; (* READFEATURES *)

      (* - - - - - - - - - - - - - - - - - - - - - - - - - *)
      procedure READBASECOUNT;
        var DUMMY:WORD;
        begin
(*!!!*)if DEBUG then writeln(output,'BASE COUNT');
          with E do begin
            SKIP(ANOFILE,12);
            read(ANOFILE,ACOMP);
            READWORD(ANOFILE,DUMMY);
            read(ANOFILE,CCOMP);
            READWORD(ANOFILE,DUMMY);
            read(ANOFILE,GCOMP);
            READWORD(ANOFILE,DUMMY);
            read(ANOFILE,TCOMP);
            readln(ANOFILE);
            CA:= CA + 1
            end
        end; (* READBASECOUNT *)

      (* - - - - - - - - - - - - - - - - - - - - - - - - - *)
      procedure READORIGIN;
        begin
          readln(ANOFILE);
          CA:= CA + 1
        end; (* READORIGIN *)

      (* - - - - - - - - - - - - - - - - - - - - - - - - - *)
      procedure READSEQUENCE(var SEQFILE:text);
        var CH:char;
            TOOLONG:boolean;
        begin
        with E do begin
          SEQLEN := 0; TOOLONG:=false;
          while (not eof(SEQFILE)) and (not (SEQFILE^='>')) do begin
            while not eoln(SEQFILE) do begin
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
              CS:= CS + 1
              end; (* eof *)
          SEQLIM:=SEQLEN;
          if TOOLONG then begin
             WRITEWORD(MESSAGE,NAME,NAME.LEN+1);
             writeln(MESSAGE,
            '>>> WARNING! Sequence length exceeds MAXSEQ-2. Seq. truncated.')
            end (* TOOLONG *)
          end (* with E *)
        end; (* READSEQUENCE *)

      begin (* RGB --------------------------------------------------------*)
        BLANKSET:=[' ',CRUNCHFLAG];
        (* Advance to the entry specified by ANOLINE in ANOFILE *)
        if CA > LOCATIONS.ANOLINE then begin reset(ANOFILE); CA:=1 end;
        while CA < LOCATIONS.ANOLINE do begin
          readln(ANOFILE);
          CA:= CA + 1
          end;

        (* Advance to the sequence specified by SEQLINE in SEQFILE *)
        if CS > LOCATIONS.SEQLINE then begin reset(SEQFILE); CS:=1 end;
        while CS < LOCATIONS.SEQLINE do begin
          readln(SEQFILE);
          CS:= CS + 1
          end;
        readln(SEQFILE);
        CS:= CS + 1;

        (* Read in the entry *)
        DONE:=false;
        while (not DONE) and (not eof(ANOFILE)) do
          if (ANOFILE^ in ['L','D','A','N','K','S','R','C','F','B','O','/'])
             then
             case ANOFILE^ of
               'L':READLOCUS;
               'D':READDEFINITION;
               'A':READACCESSION;
               'N':READNID;
               'K':READDUMMY;
               'S':begin
                     read(ANOFILE,CH);
                     case ANOFILE^ of
                      'E':READSEGMENT;
                      'O':READDUMMY
                     end
                   end;
               'R':READDUMMY;
               'C':READDUMMY;
               'F':READFEATURES;
               'B':READBASECOUNT;
               'O':READORIGIN;
               '/':DONE:=true
               end
           else READDUMMY;
        READSEQUENCE(SEQFILE);
      end; (* RGB *)

  (****************************************************************)
  (* Test to see if a coordinate is in a given fragment.          *)
  (****************************************************************)
    function WITHIN(FR:NP; COORD:integer):boolean;
      var FIRST,LAST:integer;
      function BETWEEN(X,Y,Z:integer):boolean;
        begin
          if (X<=Y) and (Y<=Z) then BETWEEN:=true
          else BETWEEN:=false
        end; (* BETWEEN *)
      begin (* WITHIN *)
        if FR<>nil then 
          with FR^ do if NTYPE=FRANODE then
            with FRA do begin
              FIRST:=START.LNUM;LAST:=FINISH.RNUM;
              case STRAND of
                 ' ':if FIRST<=LAST
                        then WITHIN:= BETWEEN(FIRST,COORD,LAST)
                        else WITHIN:= BETWEEN(FIRST,COORD,E.SEQLEN) or
                                      BETWEEN(1,COORD,LAST);
               'C','c':if FIRST >= LAST
                        then WITHIN:= BETWEEN(LAST,COORD,FIRST)
                        else WITHIN:= BETWEEN(LAST,COORD,E.SEQLEN) or
                                      BETWEEN(1,COORD,FIRST) (* iff CIRCULAR *)
                     end (* case STRAND *)
              end (* with FRA *)
          else WITHIN:=false
        else WITHIN:=false
      end; (* WITHIN *)

  (****************************************************************)
  (* Make a model of the sequence, based on criteria in datafile. *)
  (****************************************************************) 
  procedure MODEL(var E:ENTRY; TRANSCRIPTS,SITESET:KEYSET; var LL:LOCLINE;
                  var LLEN:integer; var OBJ:OBJECTS);
 
    const TABSIZE = 4;     

    var OKAY,DONE,
        TRANSLATED:boolean; (* =true if feature codes for protein *)
        TEMP:NP;
        INDENT,     (* # columns to indent in MESSAGE (see TAB,UNTAB)  *) 
        NUMN,       (* # of n's to pad CDS with if codon_start=2 or 3 *)
        CURRENT,I:integer;    (* index of the current object *)
        FEAKEY:FEATUREKEY; (* Feature key of current feature *)

   (* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *)
   (* To enhance the readability of the MESSAGE file, the indentation of
      output is incremented by TABSIZE columns for each expression 
      encountered. Thus, at the beginning of an expression, TAB adds TABSIZE
      columns to the variable INDENT, and at the end of an expression, 
      UNTAB subtracts TABSIZE columns from INDENT. *)
   procedure TAB;
     begin
       INDENT:=INDENT+TABSIZE
     end; 
   procedure UNTAB;
     begin
       INDENT:= INDENT-TABSIZE;
       if INDENT < 0 then INDENT:=0
     end;

   (* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *)
    (* Copy the location, which may span several FEATURE table lines, into
       the array LL. When DONE, TEMP will point to the line on which the
       first qualifier is found, or to the next feature, or to the end of 
       the FEATURES table.*)
    procedure ONESTRING(var LL:LOCLINE; var LLEN:integer; var TEMP:NP);
      var DONE:boolean;
          I:integer;
      begin
        DONE:=false;LLEN:=0;
        repeat
          (* Copy characters from the current location line until the end or
             until a slash (/) is reached. *)
          with TEMP^.FEA.LOCQUAL do begin
            I:=1;
            while (I<=LEN) and (STR[I]<>'/') do begin
               LLEN:= LLEN+1; LL[LLEN]:=STR[I]; I:=I+1
               end; (* while *)
               if STR[I]='/' then DONE:=true
            end; (* with *)
          (* If the location continues on another line, continue reading.*)
          if (not DONE) or (I>TEMP^.FEA.LOCQUAL.LEN) then 
            with E.FEATURETABLE do begin
              TEMP:=TEMP^.NEXT;
              if TEMP=TAIL then DONE:=true
              else with TEMP^.FEA do begin
                if LOCQUAL.STR[1]='/' then DONE:=true
                else if KEY<> BLANK then DONE:=true
                end (* TEMP^.FEA *) 
            end (* with *)
        until DONE
      end; (* ONESTRING *)

   (* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *)
   (* Extract characters from the qualifier line until a '='  or
      end of line is reached. *)
   procedure EXTRACTQUAL(var STR:CHARARRAY; var LEN,STRPOSN:integer;
                        var TESTSTR:FK);
    var I:integer;
    begin
      I:=0;
      while not(STR[STRPOSN] in ['(',')',',',':','=']) and 
            (STRPOSN<=LEN) do begin
        I:=I+1;
        if I<=MAXFK then TESTSTR[I]:=STR[STRPOSN]; (* use <= MAXFK chars.*)
        STRPOSN:=STRPOSN+1
        end;
      (* Pad the end of OPSTR with blanks *)
      while I < MAXFK do begin
        I:=I+1;
        TESTSTR[I]:=' '
        end
    end; (* EXTRACTQUAL *)

    (* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *)
    function FRAGSIZE(START,FINISH:integer;STRANDFLAG:char):integer;
      begin
        case STRANDFLAG of
           ' ':if START<FINISH
                  then FRAGSIZE:=FINISH-START+1
                  else FRAGSIZE:=(E.SEQLEN-START)+FINISH+1;(* iff CIRCULAR *)
       'c','C':if START>FINISH
                  then FRAGSIZE:=START-FINISH+1
                  else FRAGSIZE:=(E.SEQLEN-FINISH)+START+1
               end (* case STRANDFLAG *)
      end; (* FRAGSIZE *)

    (* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *)
    (* Read an integer from the location line. *)
    procedure NUMBER(var LL:LOCLINE; var LLEN, POSN,NUM:integer);
      var ORDZERO:integer;
      begin
        ORDZERO:=ord('0');
        NUM:=0;
          while (LL[POSN] in ['0'..'9']) and (POSN<=LLEN) do begin
            NUM:=(NUM*10)+(ord(LL[POSN])-ORDZERO);
            POSN:=POSN+1
            end;
         write(MESSAGE,NUM:1) (*forces use of minimal # of chars *)
       end; (* NUMBER *) 
 
   (* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *)
   (* Forward declaration for recursive procedure EVALUATE *)
    procedure EVALUATE(var LL:LOCLINE; var POSN,LLEN:integer;
                       var RESULT:LOCTYPE; var BP:BASEPOSITION; var DEST:NP); 
    forward;
    
   (* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *)
   (* Read a base position from LL.  A base position composed of a single
      number is represented by setting both  LNUM and RNUM to the same
      number.  A two base bound exists when LNUM and RNUM are different. *)
    procedure GETBP(var LL:LOCLINE; var POSN,LLEN:integer; var BP:BASEPOSITION);

      begin  (* GETBP *)
        with BP do begin
          if LL[POSN] in ['>','<'] then begin
             BOUNDFLAG:=LL[POSN]; POSN:=POSN+1 end
          else BOUNDFLAG:=' ';
          write(MESSAGE,BOUNDFLAG);
          NUMBER(LL,LLEN,POSN,LNUM);
          if (LL[POSN]='.') and (LL[POSN+1] in ['0'..'9']) and (POSN<LLEN)
            then begin
                   write(MESSAGE,'.');
                   POSN:=POSN+1;
                   NUMBER(LL,LLEN,POSN,RNUM)
                 end
          else RNUM:=LNUM;
        end (* with BP *)
      end; (* GETBP *)

   (* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *)
   (* COPYCHAR and COPYEXPRESSION are used by ONEOF and FEANAME to copy
      an expression without evaluating it. *)

     (* Copy a character from LL to EXLOC *)
     procedure COPYCHAR(var POSN,EXLOCLEN:integer; var LL,EXLOC:LOCLINE);
       begin
         EXLOCLEN:=EXLOCLEN+1;
         EXLOC[EXLOCLEN]:=LL[POSN];
         POSN:=POSN+1
       end; (* COPYCHAR *)

     (* Copy a parenthetical expression from LL to EXLOC *)       
     procedure COPYEXPRESSION(var POSN,EXLOCLEN,LLEN:integer;
                              var LL,EXLOC:LOCLINE);
       var DONE:boolean;
       begin
         (*Copy the left parenthesis '(' *)
         COPYCHAR(POSN,EXLOCLEN,LL,EXLOC);
         (* Copy the expression *)
         repeat
           while (not(LL[POSN] in ['(',')'])) and (POSN<=LLEN) do
               COPYCHAR(POSN,EXLOCLEN,LL,EXLOC);
           if LL[POSN]='(' then COPYEXPRESSION(POSN,EXLOCLEN,LLEN,LL,EXLOC)
           else if LL[POSN]=')' then begin
                   COPYCHAR(POSN,EXLOCLEN,LL,EXLOC);
                   DONE:=true
                   end
         until DONE or (POSN>LLEN)
       end; (* COPYEXPRESSION *)


    (* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *)
    (* Create a new fragment to represent a literal sequence, and store
       the literal sequence in SEQ[>SEQLIM] *)
    procedure LITERAL(var LL:LOCLINE; var POSN,LLEN:integer;
                      var RESULT:LOCTYPE; var DEST:NP);
      var I:integer;
      begin
        ADDNODE(DEST);
	DEST:=DEST^.NEXT;
        with DEST^ do begin
          NTYPE:=FRANODE;
          with FRA do begin
            LT:=baserange; RESULT:=LT;
            REPEATS:=1;
            with START do begin
              LNUM:=SEQLIM+1; RNUM:=LNUM;
              BOUNDFLAG:=' '
              end (* with START *);
            STRAND:=' ';
            POSN:=POSN+1; I:=START.LNUM;
            while LL[POSN] <> '"' do begin
              SEQ[I]:=NUC(LL[POSN]);
              write(MESSAGE,LL[POSN]);
              POSN:=POSN+1;I:=I+1
              end; (* while *)
            writeln(MESSAGE);
            SIZE:=I-START.LNUM;
            with FINISH do begin
              LNUM:= I-1;RNUM:=LNUM;
              BOUNDFLAG:=' ';
              SEQLIM:=LNUM
              end; (* with FINISH *)
            POSN:=POSN+1
            end (* with FRA *)
          end (* with DEST^ *)
      end; (* LITERAL *)

    (* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *)
    (* Complement a fragment.  This simply means to swap START and FINISH
       values and toggle the strand flag, beginning with CURRENTFRAG and
       ending with LASTFRAG. *)
    procedure COMPLEMENT(var LL:LOCLINE; var POSN,LLEN:integer;
                         var RESULT:LOCTYPE; var DEST:NP);
      var CURRENTFRAG:NP;
          DUMMY:BASEPOSITION;

      (*  Given a linked list in order 1,2,3,..n,
          invert the order to n..3,2,1 *)
      procedure INVERTSTRAND(CURRENTFRAG,DEST:NP);
        var ANODE,LEFT,RIGHT:NP;

        (* Swap the values for ends; Toggle STRAND flag. *)
        procedure INVERTFRAG(var ANODE:NP);
          var TEMPL,TEMPR:integer;
          begin
            with ANODE^.FRA do begin 
              if LT=baserange then begin  
                 TEMPR:=START.LNUM;TEMPL:=START.RNUM;
                 START.LNUM:=FINISH.RNUM; START.RNUM:=FINISH.LNUM;
                 FINISH.LNUM:=TEMPL; FINISH.RNUM:=TEMPR;
                 end;
              if STRAND=' ' then STRAND:='C' else STRAND:=' '
              end (* with ANODE^.FRA *)
          end; (* INVERTFRAG *)

        (* remove a node from between two adjacent nodes *)
        procedure CUTNODE(var LEFT,ANODE:NP);
          var LBORDER,RBORDER:NP;
          begin
            ANODE:=LEFT;
            with ANODE^ do begin
              LBORDER:= PREV; RBORDER:=NEXT;
              PREV:=nil; NEXT:=nil
              end;
            LBORDER^.NEXT:=RBORDER; RBORDER^.PREV:=LBORDER;
            LEFT:=RBORDER 
          end; (* CUTNODE *)

        (* Place ANODE to the right of RIGHT. (Kind of like Pat Buchanan.) *)
        procedure PASTENODE(var ANODE,RIGHT:NP);
          var RBORDER:NP;
          begin
            RBORDER:= RIGHT^.NEXT;
            with ANODE^ do begin
              PREV:=RIGHT; NEXT:=RBORDER;
              end;
            RIGHT^.NEXT:=ANODE; RBORDER^.PREV:=ANODE
          end; (* PASTENODE *)

        begin (* INVERTSTRAND *)
          (* Cut out leftmost node and reset LEFT to the leftmost remaining
             node. Paste it to the right of RIGHT. Done when LEFT = RIGHT *)
          LEFT:=CURRENTFRAG; RIGHT:=DEST;
          DEST:=LEFT;
          while LEFT <> RIGHT do begin 
            CUTNODE(LEFT,ANODE);
            INVERTFRAG(ANODE);
            PASTENODE(ANODE,RIGHT)
            end;  
          INVERTFRAG(RIGHT)     
        end; (* INVERTSTRAND *)

      begin (* COMPLEMENT *)
        writeln(MESSAGE);
        if LL[POSN]='(' then begin
           TAB;
           writeln(MESSAGE,' ':INDENT,'('); 
           POSN:=POSN+1;
           CURRENTFRAG:=DEST;
           EVALUATE(LL,POSN,LLEN,RESULT,DUMMY,DEST);
           if DEST <> CURRENTFRAG then begin (* one or more NODES were added *)
             CURRENTFRAG:=CURRENTFRAG^.NEXT;
             INVERTSTRAND(CURRENTFRAG,DEST)
             end;

           if LL[POSN]=')' then begin
             writeln(MESSAGE,' ':INDENT,')');
             POSN:=POSN+1
             end
           else begin
             writeln(MESSAGE,'>>> Error in database entry.');
             writeln(MESSAGE,'>>> Only one location allowed in complement',
                     ' expression');
             writeln(MESSAGE,'>>> Correct syntax: complement(location)');
             OKAY:=false
             end;
           UNTAB              
           end 
        else OKAY:=false
      end; (* COMPLEMENT *)

 (* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *) 
   (* Join a series of locations in the order specified. *) 
   procedure JOIN(var LL:LOCLINE; var POSN,LLEN:integer; var RESULT:LOCTYPE;
                  var DEST:NP); 
     var DUMMY:BASEPOSITION; 
         
       begin 
         writeln(MESSAGE);
         if  LL[POSN]='(' then begin 
             TAB;
             writeln(MESSAGE,' ':INDENT,'('); 
             POSN:=POSN+1;
             repeat
               (* Evaluate each local location in the location list *)
               EVALUATE(LL,POSN,LLEN,RESULT,DUMMY,DEST)
             until (LL[POSN]=')') or (POSN>LLEN) or not OKAY;
             if LL[POSN]=')' then begin
               writeln(MESSAGE,' ':INDENT,')');
               POSN:=POSN+1
               end;
             UNTAB
             end
          else OKAY:=false
      end; (* JOIN *)

    (* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *)
    (* Create a fragment for each element, but write a message that
       nothing is implied about the reasonableness of joining them.*)
    procedure ORDER(var LL:LOCLINE; var POSN,LLEN:integer; var RESULT:LOCTYPE;
                    var DEST:NP);
      var DUMMY:BASEPOSITION;
      begin
        writeln(MESSAGE);
        writeln(MESSAGE,'>>> ORDER: Joining fragments.');
        writeln(MESSAGE,'>>> User should verify validity of this action.');
        if LL[POSN]='(' then begin
           TAB;
           writeln(MESSAGE,' ':INDENT,'('); 
           POSN:=POSN+1;
           repeat
             (* Evaluate each local location in the location list *)
             EVALUATE(LL,POSN,LLEN,RESULT,DUMMY,DEST)
           until (LL[POSN]=')') or (POSN>LLEN) or not OKAY;
           if LL[POSN]=')' then begin
             writeln(MESSAGE,' ':INDENT,')');
             POSN:=POSN+1
             end;
           UNTAB
           end
        else OKAY:=false
      end; (* ORDER *)
   
    (* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *)
    (* Create a fragment for each element, but write a message that
       no specific order is implied.*)
    procedure GROUP(var LL:LOCLINE; var POSN,LLEN:integer; var RESULT:LOCTYPE;
                    var DEST:NP);
      var DUMMY:BASEPOSITION;
      begin
(*!!!*) writeln('GROUP');
        writeln(MESSAGE);
        writeln(MESSAGE,'>>> GROUP: Joining fragments.');
        writeln(MESSAGE,'>>> User should verify validity of this action.');
        if LL[POSN]='(' then begin
           TAB;
           writeln(MESSAGE,' ':INDENT,'('); 
           POSN:=POSN+1;
           repeat
             (* Evaluate each local location in the location list *)
             EVALUATE(LL,POSN,LLEN,RESULT,DUMMY,DEST)
           until (LL[POSN]=')') or (POSN>LLEN) or not OKAY;
           if LL[POSN]=')' then begin
             writeln(MESSAGE,' ':INDENT,')');
             POSN:=POSN+1
             end;
           UNTAB
           end
         else OKAY:=false
      end; (* GROUP *)

    (* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *)
    (* Evaluate  only the first expression, printing a message that 
       the others are alternatives. In a future version,
       alternative choices may appear as SITES in output. *)
    procedure ONEOF(var LL:LOCLINE; var POSN,LLEN:integer; var RESULT:LOCTYPE;
                    var BP:BASEPOSITION; var DEST:NP);
      var EXLOC:LOCLINE;
          EXLOCLEN,I :integer;
      begin
        writeln(MESSAGE);
        if LL[POSN]='(' then begin
           TAB;
           writeln(MESSAGE,' ':INDENT,'('); 
           POSN:=POSN+1;
           writeln(MESSAGE,'>>> THE CURRENT IMPLEMENTATION ONLY EVALUATES THE',
                ' FIRST ARGUMENT OF one-of().');
           writeln(MESSAGE,
               '    The result of this expression will appear in the output:');
           EVALUATE(LL,POSN,LLEN,RESULT,BP,DEST);
           writeln(MESSAGE,'    The following alternate choices will NOT appear',
                   ' in the output:');
           repeat
             EXLOCLEN:=0;
             COPYEXPRESSION(POSN,EXLOCLEN,LLEN,LL,EXLOC);
             write(MESSAGE,' ':INDENT);
             for I:= 1 to EXLOCLEN do write(MESSAGE,EXLOC[I]);
             writeln(MESSAGE) 
           until (LL[POSN]=')') or (POSN>LLEN) or not OKAY;
           if LL[POSN]=')' then begin
             writeln(MESSAGE,' ':INDENT,')');
             POSN:=POSN+1
             end;
           UNTAB
           end
         else OKAY:=false
      end; (* ONEOF *)

    (* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *)
    (* Replace the first location with (he contents of the second. For
       illustrative purposes, the paradigm of intron splicing is used, but
       this procedure works for all replace operations. *)
    procedure REPLACE(var LL:LOCLINE; var POSN,LLEN:integer; var DEST:NP);
      var RESULT:LOCTYPE;
          BP:BASEPOSITION;
          DONORFRAG,ACFRAG,DUMMY:NP;
          INTRON5P,INTRON3P,DONOR,ACCEPTOR:integer;
          BOTH_ENDS_FOUND:boolean;

      (* Determine the number of the base 5' to LOC *)
      function FIVEPRIME(LOC:integer;FRA:FRAGMENT;STRAND:char):integer;
        begin
          with FRA do
            if LOC=START.LNUM then FIVEPRIME:=0
            else case STRAND of
              ' ':FIVEPRIME:=LOC-1; (* input strand *)
              'c','C':FIVEPRIME:=LOC+1(* complementary strand *)
              end (* case *)
        end; (* FIVEPRIME *)
  
      (* Determine the number of the base 3' to LOC *)
      function THREEPRIME(LOC:integer;FRA:FRAGMENT;STRAND:char):integer;
         begin
          with FRA do
            if LOC=FINISH.RNUM then THREEPRIME:=0
            else case STRAND of
              ' ':THREEPRIME:= LOC+1; (* input strand *)
              'c','C':THREEPRIME:=LOC-1 (* complementary strand *)
              end; (* case *)
        end; (* THREEPRIME *)

      (* Splice an intron out of the object *)
      procedure SPLICE(var DONOR,ACCEPTOR:integer; var DONORFRAG,ACFRAG:NP);
        begin (*SPLICE *)
            if (DONOR=0) and (ACCEPTOR=0) then (* intron spans all of one
                                                  or more fragments *)
              RIDOF(DONORFRAG^.PREV,ACFRAG^.NEXT)
            else if DONOR=0 then (* intron at 5' end of object *)
              with ACFRAG^.FRA do begin
                with START do begin
                  LNUM:=ACCEPTOR;RNUM:=ACCEPTOR;
                  BOUNDFLAG:=' '
                  end; (* with START *)
                SIZE:=FRAGSIZE(ACCEPTOR,FINISH.RNUM,STRAND)
                end (* ACFRAG^.FRA *)

            else if ACCEPTOR=0 then (* intron at 3'end of object *)
              with DONORFRAG^.FRA do begin
                with FINISH do begin
                  LNUM:=DONOR; RNUM:=DONOR;
                  BOUNDFLAG:=' '
                  end; (* with FINISH *)
                SIZE:=FRAGSIZE(START.LNUM,DONOR,STRAND)  
                end (* DONORFRAG^.FRA *)

            else begin
              (* If DONOR and ACCEPTOR are in the same fragment, then 
                 insert a new fragment after DONORFRAG. *)
              if DONORFRAG=ACFRAG then begin
                ADDNODE(DONORFRAG);
                ACFRAG:=DONORFRAG^.NEXT;
		DEST:=ACFRAG; (* make sure DEST still points to last frag. *)
                with ACFRAG^ do begin
                  NTYPE:=FRANODE;
                  with FRA do begin
                    REPEATS:=1;
                    LT:=baserange;
                    FINISH:=DONORFRAG^.FRA.FINISH;
                    STRAND:=DONORFRAG^.FRA.STRAND
                    end (* with FRA *)
                  end (* with ACFRAG^ *)
                 end; (* DONORFRAG=ACFRAG *)

              (* Truncate DONORFRAG's 3' end at DONOR *)
              with DONORFRAG^.FRA.FINISH do begin
                LNUM:=DONOR; RNUM:=DONOR;
                BOUNDFLAG:=' '
                end; (* with DONORFRAG^.FRA.FINISH *)
 
              (* Truncate ACFRAG's 5' end at ACCEPTOR *) 
              with ACFRAG^.FRA.START do begin
                LNUM:=ACCEPTOR; RNUM:=ACCEPTOR;
                BOUNDFLAG:=' '
                end; (* with ACFRAG^.FRA.START *)
 
              (* Recalculate the sizes of DONORFRAG and ACFRAG *)
              with DONORFRAG^.FRA do 
                SIZE:=FRAGSIZE(START.LNUM,FINISH.RNUM,STRAND);
              if DONORFRAG<>ACFRAG then with ACFRAG^.FRA do
                SIZE:=FRAGSIZE(START.LNUM,FINISH.RNUM,STRAND)  
              end
        end; (* SPLICE *)
 
      begin (* REPLACE *)
   (*!!!*)writeln('REPLACE'); 
        writeln(MESSAGE);
        if LL[POSN]='(' then begin
          TAB;
          writeln(MESSAGE,' ':INDENT,'('); 
          POSN:=POSN+1;

          (* Evaluate the target sequence. A simple position is returned as a BP.
             A base range or between position is returned as a temporary node,
             inserted into OBJ[0] by EVALUATE.  The 5' and 3' ends of the intron
             (INTRON5P and INTRON3P) are read from this node, and RIDOF removes 
             it. *)
          DUMMY:=OBJ[0].HEAD;
          EVALUATE(LL,POSN,LLEN,RESULT,BP,DUMMY);
          (* Derive DONOR and ACCEPTOR positions *)
          if RESULT=position then begin
            INTRON5P:=BP.LNUM; INTRON3P:=BP.LNUM end
          else begin
            with DUMMY^.FRA do begin
              INTRON5P:=START.LNUM; INTRON3P:=FINISH.RNUM
              end; (* with DUMMY^.FRA *)
            RIDOF(OBJ[0].HEAD,OBJ[0].TAIL)
            end;
   
          (* Find fragments in which INTRON5P and INTRON3P occur *)
          with OBJ[CURRENT] do begin
            BOTH_ENDS_FOUND:=true;
            (* Set DONORFRAG and ACFRAG to point to the fragments
               containing positions INTRON5P and INTRON3P *) 
              DONORFRAG:=HEAD^.NEXT;
              while not WITHIN(DONORFRAG,INTRON5P) and BOTH_ENDS_FOUND do
                if DONORFRAG=TAIL then begin
                   BOTH_ENDS_FOUND:=false;
                   writeln(MESSAGE,'>>> replace target ',INTRON5P,' not found')
                   end
                else DONORFRAG:=DONORFRAG^.NEXT;

              ACFRAG:=DONORFRAG;
              while not WITHIN(ACFRAG,INTRON3P) and BOTH_ENDS_FOUND do
                if ACFRAG=TAIL then begin
                   BOTH_ENDS_FOUND:=false;
                   writeln(MESSAGE,'>>> replace target ',INTRON3P,' not found')
                   end
                else ACFRAG:=ACFRAG^.NEXT
            end; (* with OBJ[CURRENT] *)

          if BOTH_ENDS_FOUND then begin
            (* Calculate DONOR and ACCEPTOR positions *)
            if RESULT=betweenposition then begin
              DONOR:=INTRON5P;
              ACCEPTOR:=INTRON3P
              end
            else begin
              DONOR:=FIVEPRIME(INTRON5P,DONORFRAG^.FRA,' ');
              ACCEPTOR:=THREEPRIME(INTRON3P,ACFRAG^.FRA,' ')
              end;

            (* Remove the target *)
            SPLICE(DONOR,ACCEPTOR,DONORFRAG,ACFRAG);
        
            (* Insert the resultant fragment(s) *)
            if (LL[POSN]='"') and (LL[POSN+1]='"') then (*deletion,ie.null frag.*)
               POSN:=POSN+2
            else EVALUATE(LL,POSN,LLEN,RESULT,BP,DONORFRAG);
            if (LL[POSN]=')') then begin
               writeln(MESSAGE,' ':INDENT,')');
               POSN:=POSN+1
               end;
             end (* BOTH_ENDS_FOUND *)
          else begin OKAY:=false; POSN:=LLEN+1 end;
          UNTAB
          end
        else OKAY:=false
      end; (* REPLACE *)

    (* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *)
    (*  THIS IS AN EXTENSION TO THE FEATURES TABLE LANGUAGE!!!
        poly(<absolute_location>|<literal>|<feature_name>,n)
        is evaluated to mean that the location is repeated n times.
        Note: in the current implementation, complex expressions involving
        functional_operators are not supported. 

        Essentially, what this procedure does is to evaluate the first
        term in the expression, creating a fragment (.FRA). Then, the 
        second term is evaluated to an integer, telling how many times
        the fragment is to be written by WRITEOBJ.*)
    procedure POLY(var LL:LOCLINE; var POSN,LLEN:integer; var RESULT:LOCTYPE;
                   var DEST:NP);
      var DUMMY:BASEPOSITION;
      begin
        writeln(MESSAGE);
        if LL[POSN]='(' then begin
           TAB;
           writeln(MESSAGE,' ':INDENT,'('); 
           POSN:=POSN+1;

           (* Evaluate the expression *)
           EVALUATE(LL,POSN,LLEN,RESULT,DUMMY,DEST);
           TAB;
           (* Read the and set the repeat factor, REPEATS *)
           write(MESSAGE,' ':INDENT);
           NUMBER(LL,LLEN,POSN,DEST^.FRA.REPEATS);
           writeln(MESSAGE);
           UNTAB;
           
           if LL[POSN]=')' then begin
             writeln(MESSAGE,' ':INDENT,')');
             POSN:=POSN+1
             end
           else begin
              writeln(MESSAGE,'>>> TOO MANY ARGUMENTS TO POLY');
              OKAY:=false
              end;
           UNTAB
           end
         else OKAY:=false
      end; (* POLY *)

  
   (* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *)
   (* Resolve a reference to a labeled feature within the current entry.*)
   procedure INTRAENTRY(var OPSTR,QUALSTR:FK; var RESULT:LOCTYPE; var DEST:NP);
     var LPTR:NP;
         LL:LOCLINE;
         POSN, (*note: this refers to position in labeled feature, NOT
                        the calling feature. *)
         LLEN:integer;
         TESTSTR:FK;
         DUMMY:BASEPOSITION;
         
     (* Given a label, find the corresponding feature within the 
        current entry. *)
     procedure FINDQUALIFIER(var OPSTR,QUALSTR:FK; var LPTR:NP);
       var FOUND:boolean;
           STRPOSN:integer;

        begin (* FINDQUALIFIER *)
(*!!!*)if DEBUG then writeln('FINDQUALIFIER');
         with E.FEATURETABLE do begin 
           LPTR:= HEAD^.NEXT;
           FOUND:=false;
           while (not FOUND and (LPTR<>TAIL)) do begin 
             if LPTR^.FEA.KEY = BLANK then  (*ie. qualifier line *)
               (* OPSTR may occur anywhere on LOCQUAL line *)
               with LPTR^.FEA.LOCQUAL do begin
(*!!!*)if DEBUG then writeln('look for /',OPSTR,'=');
                 STRPOSN:=1;
                 while (STRPOSN<LEN) and (not FOUND) do begin
                   while (STR[STRPOSN]<>'/') and (STRPOSN<LEN) 
                         do STRPOSN:=STRPOSN+1;
                   if STR[STRPOSN]='/' then begin
                     STRPOSN:=STRPOSN+1;
                     EXTRACTQUAL(STR,LEN,STRPOSN,TESTSTR);
                     STRPOSN:=STRPOSN+1; (* read past '=' *)
                     if TESTSTR=OPSTR then begin (* look for QUALSTR *)
(*!!!*)if DEBUG then writeln('look for ',QUALSTR);
                       EXTRACTQUAL(STR,LEN,STRPOSN,TESTSTR);
                       if TESTSTR=QUALSTR then FOUND:=true
                       end (* TESTSTR=OPSTR *)
                     end (* STR[STRPOSN]='/' *)
                   end (* STRPOSN<LEN and not FOUND *)
               end; (* with LPTR^.FEA.LOCQUAL *)
             if not FOUND then LPTR:=LPTR^.NEXT
             end (* while not FOUND and LPTR <> TAIL *)
         end (* with E.FEATURETABLE *)
       end; (* FINDQUALIFIER *)

     begin (* INTRAENTRY *)
(*!!!*)if DEBUG then writeln('INTRAENTRY');
       writeln(MESSAGE);
         (* Find the feature by its label. *)
         FINDQUALIFIER(OPSTR,QUALSTR,LPTR);

         if LPTR <> E.FEATURETABLE.TAIL then begin
           (* Back up to location part of the feature. *)
           while LPTR^.FEA.KEY = BLANK do LPTR:=LPTR^.PREV;

           (* Evaluate the feature *)
           ONESTRING(LL,LLEN,LPTR);
           POSN:=1;
           EVALUATE(LL,POSN,LLEN,RESULT,DUMMY,DEST);

           (* write qualifier lines to MESSAGE file *)
           LPTR:=LPTR^.NEXT;
           while (LPTR^.FEA.KEY = BLANK) and (LPTR <> E.FEATURETABLE.TAIL)
                  do begin
             write(MESSAGE,' ':INDENT);
             with LPTR^.FEA do WRITELINE(MESSAGE,LOCQUAL,LOCQUAL.LEN);
             writeln(MESSAGE);
             LPTR:=LPTR^.NEXT
             end
           end
         else begin
           OKAY:=false;
           writeln(MESSAGE,'>>> ',QUALSTR,' not found')
           end
     end; (* INTRAENTRY *)

   (* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *)
   procedure FEANAME(var LL:LOCLINE; var POSN,LLEN:integer; OPSTR:FK;
                     var DEST:NP);
     var I,STARTLOCATION:integer;
     begin
       ADDNODE(DEST);
       DEST:=DEST^.NEXT;
       with DEST^ do begin
         NTYPE:=FRANODE;
         with FRA do begin
           LT:=featurename;
           REPEATS:=1;

           (* Copy OPSTR into EXLOC *)
           I:=1;EXLOCLEN:=0;
           repeat
             EXLOCLEN:=EXLOCLEN+1;
             EXLOC[EXLOCLEN]:=OPSTR[I];
             I:=I+1
           until OPSTR[I]=' ';
           STARTLOCATION:=I;
        
           (* Copy characters from LL to EXLOC until an expression is found
              or the end of the line is reached *)
           while not(LL[POSN] in [',','(',')']) and (POSN <=  LLEN) do
             COPYCHAR(POSN,EXLOCLEN,LL,EXLOC);
           
           (* Copy an expression from LL to EXLOC, if there is one *)
           if LL[POSN]='(' then COPYEXPRESSION(POSN,EXLOCLEN,LLEN,LL,EXLOC);
           for I:= STARTLOCATION to EXLOCLEN do write(MESSAGE,EXLOC[I]);
           writeln(MESSAGE)
           end (* with FRA *)
         end (* with DEST^ *)
     end; (* FEANAME *)

    (* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *)
    (*  Evaluate a location, returning a linked list after DEST.
        The procedure that calls evaluate decides what to do with the list.*)
    procedure EVALUATE;
      var OPSTR,QUALSTR:FK;
          DUMMY,TAIL:NP; (* DUMMY is a temporary node added when evaluating 2nd 
                            number of base_range or between_position. TAIL
                            is the next node.*)

      (* This procedure checks for the special case in which a label begins
         with numeric characters. *)
      function NOTLABEL(POSN:integer):boolean;
        begin
          if LL[POSN] in ['<','>'] then POSN:=POSN+1;
          while (LL[POSN] in ['0'..'9']) and (POSN<=LLEN) do POSN:=POSN+1;
          if not (LL[POSN] in ['.',',','(',')','^']) and (POSN<=LLEN) 
             then NOTLABEL:=false
             else NOTLABEL:=true
        end; (* NOTLABEL *) 

      (* Extract characters from the location line until a comma, left or
         right parenthesis, colon  or end of line is reached. *)
      procedure EXTRACT(var POSN:integer; var OPSTR:FK);
        var I:integer;
        begin
          I:=0;
          while not(LL[POSN] in ['(',')',',',':','=']) and (POSN<=LLEN) do begin
            I:=I+1;
            if I<=MAXFK then OPSTR[I]:=LL[POSN]; (* use <= MAXFK chars.*)
            POSN:=POSN+1
            end;

          (* Pad the end of OPSTR with blanks *)
          while I < MAXFK do begin
            I:=I+1;
            OPSTR[I]:=' '
            end
        end; (* EXTRACT *)
            
      begin (* EVALUATE *) 
      TAB;
      write(MESSAGE,' ':INDENT);

      if POSN<=LLEN then 
        (* Nested expression, probably base_position *)
        if LL[POSN]='(' then begin
          writeln(MESSAGE,'(');
          POSN:=POSN+1;
          EVALUATE(LL,POSN,LLEN,RESULT,BP,DEST);
          if LL[POSN]=')' then begin
             UNTAB;
             write(MESSAGE,' ':INDENT,')');
             POSN:=POSN+1
             end
           else OKAY:=false
          end
          
        (* base_position *)      
        else if (LL[POSN] in ['0'..'9','>','<']) and NOTLABEL(POSN) then begin
            RESULT:=position;
            GETBP(LL,POSN,LLEN,BP);
            ADDNODE(DEST);
            DEST:=DEST^.NEXT;
            DEST^.NTYPE:=FRANODE;
            with DEST^.FRA do begin
               LT:=position;
               REPEATS:=1;
               START:=BP;
               FINISH:=BP;
               STRAND:=' ';
               SIZE:=1
               end (* with DEST^.FRA *)
            end (* local location *)

        (* literal *)
        else if LL[POSN]='"' then LITERAL(LL,POSN,LLEN,RESULT,DEST)
  
        (* location operator or featurename *)
        else begin
           if LL[POSN]='/' then begin
              write(MESSAGE,'/');
              POSN:=POSN+1
              end;
           EXTRACT(POSN,OPSTR);
           write(MESSAGE,OPSTR);
           if OPSTR='complement     ' then COMPLEMENT(LL,POSN,LLEN,RESULT,DEST)
           else if OPSTR='join           ' then JOIN(LL,POSN,LLEN,RESULT,DEST)
           else if OPSTR='order          ' then ORDER(LL,POSN,LLEN,RESULT,DEST)
           else if OPSTR='group          ' then GROUP(LL,POSN,LLEN,RESULT,DEST)
           else if OPSTR='one-of         ' then ONEOF(LL,POSN,LLEN,RESULT,
                                                      BP,DEST)
           else if OPSTR='replace        ' then REPLACE(LL,POSN,LLEN,DEST)
           else if OPSTR='poly           ' then POLY(LL,POSN,LLEN,RESULT,DEST)
           else if LL[POSN] =':' then  FEANAME(LL,POSN,LLEN,OPSTR,DEST)
           else begin
                 if LL[POSN]='=' then begin
                   write(MESSAGE,'=');
                   POSN:=POSN+1;
                   EXTRACT(POSN,QUALSTR);
                   write(MESSAGE,QUALSTR)
                   end
                 else begin
                      QUALSTR:=OPSTR;
                      OPSTR:='label          '
                      end;
                INTRAENTRY(OPSTR,QUALSTR,RESULT,DEST)
                end  (* INTRAENTRY *)
           end; (* location operator or featurename *)

        (* Evaluate the second part of a base range or between position,
           if any. *)
        if POSN<LLEN then
          if LL[POSN] in ['.','^'] then begin
             POSN:=POSN+1;
             with DEST^.FRA do begin
                if LL[POSN]='.' then begin
                   LT:=baserange; POSN:=POSN+1 end
                else LT:=betweenposition; 

                (* Get the next number *)
                EVALUATE(LL,POSN,LLEN,RESULT,BP,DEST);
               end; (* with DEST^.FRA *)

             (* Assign pointers so that DUMMY points to the temporary
                node and TAIL to the next one after that. This makes
                it possible to call RIDOF to remove the DUMMY node.*)
             DUMMY:=DEST; TAIL:=DUMMY^.NEXT; DEST:=DUMMY^.PREV;

             (* Update DEST with information from the second part of
                base_range or between_position. *)
             with DEST^.FRA do begin
                FINISH:=BP;
                SIZE:=FRAGSIZE(START.LNUM,FINISH.RNUM,STRAND);
                RESULT:=LT
                end; (* with DEST^.FRA *)
             RIDOF(DEST,TAIL); (* remove dummy fragment *)
             writeln(MESSAGE);
             UNTAB;
             end (* base range or between position *)
            else begin
                 UNTAB;
                 writeln(MESSAGE)
                 end
          else begin
               UNTAB;
               writeln(MESSAGE)
               end;

       (* Read past delimiters to next expression or end of LL array *)
       if POSN<=LLEN then
          if  LL[POSN] = ',' then POSN:=POSN+1
      end; (* EVALUATE *)

   (* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *)
   (* Add a fragment to the begining of a CDS, sig_peptide, or mat_peptide,
      containing the result of the literal expression "n" or "nn", which will
      result in adding n's to the sequence, to put it in the correct reading
      frame. *)
   procedure PADCODON(var QUAL:LINE; var NUMN:integer);
     var POSN,LLEN:integer;
         READINGFRAME:char;
         TESTSTR:FK;
         LL:LOCLINE;
         HEADFRAG:NP;
         RESULT:LOCTYPE;  (* dummy *)
         BP:BASEPOSITION; (* dummy *)
     begin 
       with QUAL do 
         if STR[1]='/' then begin
           (* Extract the qualifier TESTSTR from QUAL, and if it is 'codon_
              start', set READINGFRAME to the position indicated (1,2 or 3). If 
              READINGFRAME is 1, do nothing. If READINGFRAME = 2 or 3, create a
              LOCLINE containing a literal expression, to be evaluated by
              EVALUATE, and inserted after the HEAD of OBJ[CURRENT] *)
           POSN:=2;
           EXTRACTQUAL(STR,LEN,POSN,TESTSTR);
           if TESTSTR='codon_start    ' then begin
             READINGFRAME:= STR[POSN+1];
             if READINGFRAME in ['2','3'] then begin
               write(MESSAGE,'n''s added to 5'' end to preserve reading frame: ');
               (* Create a LOCLINE to be evaluated: "n" or "nn"  *)
               LL[1]:= '"';
               LL[2]:= 'n';
               if READINGFRAME='3' then begin
                  LLEN:=3; 
                  LL[3]:= '"';
                  NUMN:=1
                  end
               else begin
                  LLEN:=4;
                  LL[3]:='n';
                  LL[4]:='"';
                  NUMN:=2
                 end;
               (* Evaluate the expression, and add it to the head of the 
                  object *)
               HEADFRAG:=OBJ[CURRENT].HEAD;
               POSN:=1;
               EVALUATE(LL,POSN,LLEN,RESULT,BP,HEADFRAG);
               end (* READINGFRAME in ['2','3'] *)
             else NUMN:=0
             end (* TESTSTR='codon_start' *)
            end (* STR[1]='/' *)
     end; (* PADCODON *)

   (* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *)
       (*  <title>::=  <locus>:<label>|
                       <locus>:<feature key>  *)
  procedure MAKETITLE(LPTR:NP; FEAKEY:FEATUREKEY; var OBJ:OBJECT);
     var  I,J,POSN:integer;
          CH:char;

     (* Read the feature label, if any, from the qualifier line. *)
     procedure FINDLABEL(QUALIFIER:LINE; var FEALAB:WORD; var NOLABEL:boolean);
       var TESTSTR:FK;
       begin (* FINDLABEL *)
         POSN:=2; (* read past '/' *)
         with QUALIFIER do EXTRACTQUAL(STR,LEN,POSN,TESTSTR);
         if TESTSTR='label          ' then begin
            NOLABEL:=false;
            POSN:=POSN+1; (* read past '=' *)
            (* Copy the label to FEALAB *)
            with FEALAB do begin
                 LEN:=0;
                 for I:= POSN to QUALIFIER.LEN do begin
                     LEN:= LEN + 1; 
                     STR[LEN]:= QUALIFIER.STR[I]
                     end; (* for *)
                 for I:= LEN+1 to MAXWORD do STR[I]:=' ' 
                 end (* with FEALAB *)         
            end (* if *)
       end; (* FINDLABEL *)

     begin (* MAKETITLE *)
       with OBJ do begin
         (* Title begins with the locus name *)
         TITLE:=E.NAME;
        
         (* Find the feature label, if any. *)
         NOLABEL:=true;
         while (LPTR^.FEA.KEY = BLANK) and (LPTR <> E.FEATURETABLE.TAIL) and
                NOLABEL do begin
           with LPTR^.FEA do FINDLABEL(LOCQUAL,FEALAB,NOLABEL);
           LPTR:=LPTR^.NEXT
           end;

          (* If there isn't a label, add featurekey to the title, otherwise,
             add the label. *)
          with TITLE do begin
            LEN:=LEN+1;
            STR[LEN]:=':';
            if NOLABEL then begin
               CH:= KEYSTR[FEAKEY][1]; J:=1;
               while (CH <> ' ') and (J < MAXFK) do begin
                   LEN:=LEN+1;
                   STR[LEN]:= CH;
                   J:= J + 1;
                   if J<= MAXFK then CH:= KEYSTR[FEAKEY][J]
                   end (* while *)
               end (* NOLABEL *)
             else begin
               for J:= 1 to FEALAB.LEN do begin
                  LEN:=LEN+1;
                  STR[LEN]:= FEALAB.STR[J]
                  end;(* for *)
               for J:= LEN+1 to MAXWORD do STR[J]:=' '
               end (*else *)
             end (* with TITLE *)
        end (* with OBJ *)
     end; (* MAKETITLE *)

   (* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *)
   procedure BUILD(var LL:LOCLINE; var LLEN:integer; var OBJ:OBJECT);
     var POSN:integer;
         RESULT:LOCTYPE;
         LASTFRAG:NP;
         BP:BASEPOSITION;
     begin
       with OBJ do begin
         if RESOLVEXP then begin
            WRITEWORD(MESSAGE,ACCESSION,ACCESSION.LEN);
            writeln(MESSAGE,':',CURRENT:1);
            end
         else begin
           WRITEWORD(MESSAGE,TITLE,TITLE.LEN);
           if NOLABEL then write(MESSAGE,CURRENT:1);
           writeln(MESSAGE)
           end;
         OKAY:=true;
         POSN:=1;
         INDENT:=0;
         LASTFRAG:=HEAD;
         EVALUATE(LL,POSN,LLEN,RESULT,BP,LASTFRAG);

         (* Although not permitted by the Features Table Definition,
            GenBank has been known to allow locations that evaluate to a
            between position. Since this is biologically meaningless, it
            must be detected and eliminated. *)
         if RESULT = betweenposition then begin
            writeln(MESSAGE,'>>> FEATURE EVALUATES TO between_position.');
            writeln(MESSAGE,'>>> NO SEQUENCE RESULTS FROM THIS EXPRESSION.');
            OKAY:=false
            end;
         writeln(MESSAGE);
         if not OKAY then begin
            RIDOF(HEAD,TAIL);
            TITLE.LEN:=0;
            FEALAB.LEN:=0;
            NOLABEL:=true;
            NUMOBJ:=NUMOBJ-1
            end (* if not OKAY *)
         end (* with *)
     end; (* BUILD *)

  (* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - *)
  (* this is a non-functional skeleton procedure, at present *)
  procedure FINDSITE(var LL:LOCLINE; var LLEN:integer; var OBJ:OBJECT);
     var POSN:integer;
         RESULT:LOCTYPE;
         LASTFRAG:NP;
         BP:BASEPOSITION;
     begin
       OKAY:=true;
       POSN:=1;
       INDENT:=0;
       LASTFRAG:=OBJ.HEAD;
       EVALUATE(LL,POSN,LLEN,RESULT,BP,LASTFRAG);
       writeln(MESSAGE);
       if not OKAY then begin
          RIDOF(OBJ.HEAD,OBJ.TAIL);
          NUMOBJ:=NUMOBJ-1
          end (* if not OKAY *)
     end; (* FINDSITE *)
 
    begin  (* MODEL -------------------------------------------------------*)

      (* If RESOLVEXP is true, then simply evaluate the location found in LL.
         Otherwise, make a pass through the FEATURES table, creating one or
         more objects corresponding to TESTSET.  Finally, make a second pass
         through the data, searching for each site in SITESET and assigning 
         it to a particular object. *)
      with E.FEATURETABLE do begin

          (* 1st Pass - Create objects *)
          if RESOLVEXP then begin
             NUMOBJ:=1;
             CURRENT:=NUMOBJ; (* any procedure can refer to CURRENT to find
                                 out the current OBJECT.  This means that
                                 OBJ doesn't have to be passed to each proc.*)
             BUILD(LL,LLEN,OBJ[CURRENT]);
             writeln(MESSAGE,'//----------------------------------------------')
             end (* RESOLVEXP *)
          else begin
           NUMOBJ:=0;
           TEMP:=HEAD^.NEXT;
            while TEMP <> TAIL do begin
              FEAKEY:= TEMP^.FEA.KEY;
              if FEAKEY in TRANSCRIPTS then begin (* an object *)
                if FEAKEY in [CDS,MATPEPTIDE,SIGPEPTIDE] then
                   TRANSLATED:=true
                else TRANSLATED:=false; 
                NUMN:=0;
                NUMOBJ:=NUMOBJ+1;
                CURRENT:=NUMOBJ; (* any procedure can refer to CURRENT to find
                                    out the current OBJECT.  This means that
                                  OBJ doesn't have to be passed to each proc.*)

                (* Re-code the location to a single line for easy parsing *)
                ONESTRING(LL,LLEN,TEMP);

                (* Each object carries a unique title, to identify the object
                   in all output files.*)
                MAKETITLE(TEMP, FEAKEY, OBJ[CURRENT]);

                (* Build the object *)
                BUILD(LL,LLEN,OBJ[CURRENT]);

               (* Write qualifier lines to message file *)
                DONE:=false;
                while not DONE do
                  if TEMP=TAIL then DONE:=true 
                  else
                    if TEMP^.FEA.KEY in [BLANK,OTHER] then begin
                      (* Write the qualifier *)
                      WRITELINE(MESSAGE,TEMP^.FEA.LOCQUAL,TEMP^.FEA.LOCQUAL.LEN); 
                      writeln(MESSAGE);

                      (* If the qualifier is /codon_start, then add a fragment to
                         the head of the object, containing n's to put the object
                         in the correct reading frame. *)
                      if TRANSLATED and PAD5PRIME then
                            PADCODON(TEMP^.FEA.LOCQUAL,NUMN);
                      TEMP:=TEMP^.NEXT
                     end
                   else DONE:=true;
                   writeln(MESSAGE,'//------------------------------',
                                '----------------');
                (* output to EXPFILE *)
                if OKAY then with OBJ[CURRENT] do begin 
                   (*  ><title> *)
                   write(EXPFILE,'>');
                   WRITEWORD(EXPFILE,TITLE,TITLE.LEN);
                   if NOLABEL then write(EXPFILE,CURRENT:1);
                   writeln(EXPFILE);

                   (* if codon_start<>1, write n's before expression *)
                   if NUMN in [1,2] then begin
                      for I:= 1 to NUMN do write(EXPFILE,'n');
                      writeln(EXPFILE)
                      end;

                   (* Write primary accession number etc.: @<accession>: *)
                   write(EXPFILE,'@');
                   with E.ACCESSION.HEAD^.NEXT^ do
                        WRITEWORD(EXPFILE,WOR,WOR.LEN);
                   write(EXPFILE,':');
                   (* <expression> *)
                   if NOLABEL then for I:= 1 to LLEN do write(EXPFILE,LL[I])
                   else WRITEWORD(EXPFILE,FEALAB,FEALAB.LEN);      
                   writeln(EXPFILE)
                   end (* if OKAY *)
                end (*  an object *)
              else TEMP:= TEMP^.NEXT         
              end (* TEMP <> TAIL *)
            end; (* 1st Pass *)

          (* 2nd Pass -  Assign sites to objects *)
          if (NUMOBJ > 0) and (SITESET <> []) then begin
            writeln(MESSAGE,'SITES:');
            TEMP:=HEAD^.NEXT;
            while TEMP <> TAIL do begin
              if TEMP^.FEA.KEY in SITESET then begin (* a site *)
                FINDSITE(LL,LLEN,OBJ[0]);
                end (*  a site *)
              else TEMP:= TEMP^.NEXT         
              end (* TEMP <> TAIL *)
            end (* NUMOBJ > 0 *)
        end (* with E.FEATURETABLE *)
    end; (* MODEL *)

  (****************************************)
  (* Write the object to the output file. *)
  (****************************************)
  procedure WRITEOBJECT(var F:text; OBJ:OBJECTS);
    var FIRST,LAST,I,J,K,WIDTH:integer;
        TEMP:NP;

    procedure WRITEFRAG(var F:text; START,FINISH,SENSE,REPEATS:integer;
                        NOOKIE:NAR; var WIDTH:integer);
      var I,J:integer;
      begin
        for J:= 1 to REPEATS do begin
          I:=START;
          repeat
            write(F,NUCHAR[NOOKIE[SEQ[I]]]);
            WIDTH:=WIDTH+1;
            if WIDTH=50 then begin
               writeln(F);
               WIDTH:=0
               end;
            I:=I+SENSE
          until I=FINISH+SENSE
          end (* for J *)
      end; (* WRITEFRAG *)

    begin

      for I:= 1 to NUMOBJ do 
        with OBJ[I] do begin
          (* Write the object to F *)
          if not RESOLVEXP then begin
            write(F,'>');WRITEWORD(F,TITLE,TITLE.LEN);
            if NOLABEL then write(F,I:1);
            writeln(F)
            end;
          TEMP:=HEAD^.NEXT;
          repeat
            with TEMP^.FRA do if LT in [baserange,position] then begin
              WIDTH:=0; FIRST:=START.LNUM; LAST:=FINISH.RNUM;
              case STRAND of
                      ' ':if FIRST<=LAST 
                            then WRITEFRAG(F,FIRST,LAST,1,REPEATS,INUC,WIDTH)
                            else begin
                                 WRITEFRAG(F,FIRST,E.SEQLEN,1,REPEATS,INUC,WIDTH);
                                 WRITEFRAG(F,1,LAST,1,REPEATS,INUC,WIDTH)
                                 end;
                  'C','c':if FIRST>=LAST 
                            then WRITEFRAG(F,FIRST,LAST,-1,REPEATS,COMP,WIDTH)
                            else begin
                                 WRITEFRAG(F,FIRST,1,-1,REPEATS,COMP,WIDTH);
                                 WRITEFRAG(F,E.SEQLEN,LAST,-1,REPEATS,COMP,WIDTH)
                                 end
                  end (* case STRAND *)
                end (* baserange *)
              else begin (* featurename *)
                for K:= 1 to REPEATS do begin
                  write(F,'@');
                  for J:= 1 to EXLOCLEN do write(F,EXLOC[J])
                  end (* for K *)
                end; (* featurename *)  
        
            writeln(F);
            TEMP:=TEMP^.NEXT
          until TEMP=TAIL
          end
    end; (* WRITEOBJECT *)


  (****************************************)
  (* Write the GenBank entry to the output*)
  (****************************************)
  procedure WRITEGB(var F:text; E:ENTRY);

    type KEYWORD = packed array[1..12] of char;

    procedure WRITELIST(HEAD,TAIL:NP;KW:KEYWORD);
      var TEMP:NP;
          LINEWIDTH:integer;
      begin
      if HEAD^.NEXT <> TAIL then begin
        write(F,KW);
        TEMP:=HEAD^.NEXT;
        case TEMP^.NTYPE of
          WNODE:begin
                  LINEWIDTH:=12;
                  repeat
                    with TEMP^ do begin
                      LINEWIDTH:= LINEWIDTH+WOR.LEN+1;
                      if LINEWIDTH > 80 then begin
                        writeln(F);write(F,' ':12);LINEWIDTH:=12+WOR.LEN+1 end;
                      WRITEWORD(F,WOR,WOR.LEN); write(F,' ')
                      end;
                    TEMP:= TEMP^.NEXT
                  until TEMP = TAIL;
                  writeln(F)
                end;
          LNODE:begin
                  WRITELINE(F,TEMP^.LIN,TEMP^.LIN.LEN);writeln(F);
                  TEMP:=TEMP^.NEXT;
                  while TEMP <> TAIL do begin
                    write(F,' ':12);
                    WRITELINE(F,TEMP^.LIN,TEMP^.LIN.LEN);writeln(F);
                    TEMP:= TEMP^.NEXT
                    end
                end;
          FRANODE:;
          FEANODE:begin
                    if KW='FEATURES    ' then 
                      writeln(F,'        Location/Qualifiers')
                    else writeln(F);
                    repeat
                      with TEMP^.FEA do begin
                        if KEY=OTHER then write(F,'    ',KEYSTRING,' ')
                        else write(F,'    ',KEYSTR[KEY],' ');
                        WRITELINE(F,LOCQUAL,LOCQUAL.LEN);writeln(F) 
                        end; (* with TEMP.FEA *)
                      TEMP:= TEMP^.NEXT
                    until TEMP = TAIL
                end (* FEANODE *)
          end(* case *)
        end (* HEAD *)
      end; (* WRITELIST *)

    begin (* WRITEGB *)
      with E do begin
      (* Write the annotation data *)
      write(F,'LOCUS       ');WRITEWORD(F,NAME,18);writeln(F,SEQLEN:11,' bp');
      WRITELIST(DEFINITION.HEAD,DEFINITION.TAIL,'DEFINITION  ');
      WRITELIST(ACCESSION.HEAD,ACCESSION.TAIL,'ACCESSION   ');
      (* NID field to be removed from GenBank by Dec. 1999. However,
         this statement will write NID fields from old files, if they
         exist. *)
      if NID.LEN > 0 then begin
         write(F,'NID         ');WRITEWORD(F,NID,NID.LEN);writeln(F) end;
      if SEGTOTAL > 0 then
         writeln(F,'SEGMENT     ',SEGNUMBER:1,' of ',SEGTOTAL:1);
      WRITELIST(FEATURETABLE.HEAD,FEATURETABLE.TAIL,'FEATURES    ');
      writeln(F,'BASECOUNT   ',ACOMP:7,' a',CCOMP:7,' c',GCOMP:7,' g',
              TCOMP:7,' t');
     
      writeln(F,'//')
      end (* with E *)
    end; (* WRITEGB *)

  (*****************************************************)
  (* Reinitialize the entry.                           *)
  (*****************************************************)
  procedure REINIT(var E:ENTRY);
    var I:integer;
    begin
      with E do begin
        NAME.LEN:=0; SEQLEN:=0; ENTRYSTATUS.LEN:=0;
        RIDOF(DEFINITION.HEAD,DEFINITION.TAIL);
        RIDOF(ACCESSION.HEAD,ACCESSION.TAIL);
        NID.LEN:=0;
        SEGTOTAL:=0;
        RIDOF(FEATURETABLE.HEAD,FEATURETABLE.TAIL);
        for I:= 1 to NUMOBJ do
           with OBJ[I] do begin
             RIDOF(HEAD,TAIL);
             TITLE.LEN:=0;
             FEALAB.LEN:=0;
             NOLABEL:=true
             end;
        NUMOBJ:=0
        end
    end; (* REINIT *)

   
  (* -----------------------------------------------------------  *)
  (* ----------------- MAIN  PROCEDURE  -----------------------  *)
    begin
       (* Initialize parameters. *)
       INITIALIZE;
       
       (* Read options from command line and set FILES,RESOLVEXP & ACNO *)
       ARGNUM:=STARTARGNUM; 
       READOPTIONS(ARGNUM,FILES,RESOLVEXP,ACNO);
      
       (* Open files.*)
       OPENFILES(ARGNUM,INFILE,NAMEFILE,ANOFILE,SEQFILE,INDFILE,MESSAGE,OUTFILE,
                 EXPFILE);

       (* Read input file containing instructions for processing each
             entry.*)
       READINSTRUCTIONS(INFILE,TRANSCRIPTS,SITESET,OBJECTTYPE);
 
       writeln(MESSAGE,VERSION);
       writeln(MESSAGE,'Please cite: Fristensky B. (1993)',
                       ' Feature expressions:');
       writeln(MESSAGE,'creating and manipulating sequence datasets.');
       writeln(MESSAGE,'Nucl. Acids Res. 21:5997-6003');
       writeln(MESSAGE);
       CA:=1; CS:=1;

       if RESOLVEXP then  (*ie. -r *) 
         (*- - - - Resolve objects from external references  - - - - - -*)
          while not eof(NAMEFILE) do begin
          (* Copy lines from NAMEFILE to OUTFILE. These may be objects
             resolved in previous runs of GETOB *)
             while (NAMEFILE^ <> '@') and (not eof(NAMEFILE)) do begin
                while not eoln(NAMEFILE) do begin
                  read(NAMEFILE,CH); write(OUTFILE,CH)
                  end;
                if not eof(NAMEFILE) then readln(NAMEFILE);
                writeln(OUTFILE)
                end; (* NAMEFILE^<>'@' *)

          (* Process external reference. These are references left un-
              resolved in a previous run of GETOB. *)
             if NAMEFILE^='@' then begin
               READLINE(NAMEFILE,INDIRREF);
               PARSEREF(INDIRREF,ACCESSION,LL,LLEN,OKAY);
               FOUND:=false;
               if OKAY then FINDDATA(INDFILE,ACCESSION,LOCATIONS,FOUND);
               if FOUND and OKAY then begin
                  RGB(ANOFILE,SEQFILE,LOCATIONS,E);
                  MODEL(E,TRANSCRIPTS,SITESET,LL,LLEN,OBJ);
                  WRITEOBJECT(OUTFILE,OBJ);
                  REINIT(E)
                  end (* FOUND *)
               else begin (* print message; comment out expression in OUTFILE*)
                      write(MESSAGE,'>>> Can''t resolve indirect reference: ');
                      WRITELINE(MESSAGE,INDIRREF,INDIRREF.LEN);
                      writeln(MESSAGE);
                      write(OUTFILE,';');
                      WRITELINE(OUTFILE,INDIRREF,INDIRREF.LEN);
                      writeln(OUTFILE)
                    end
               end (* NAMEFILE^='@' *)
            end (* eof(NAMEFILE) *)

       else begin 
         (* - - - - Search for specified features in each entry - - - - *) 

         (* For each locus in NAMEFILE, find the locations of the annotation 
            and sequence parts in ANOFILE and SEQFILE, and then get the data
            and write to OUTFILE format specified by INFILE. *)
         while not eof(NAMEFILE) do begin
           READNAME(NAMEFILE,SEQNAME);
           FOUND:= false;
           if SEQNAME.LEN > 0 then FINDDATA(INDFILE,SEQNAME,LOCATIONS,FOUND);
           if FOUND then begin
             if FILES then begin
                MAKEFN(SEQNAME,FILENAME);         
                rewrite(OUTFILE,FILENAME.STR)
                end; (* FILES *)
             RGB(ANOFILE,SEQFILE,LOCATIONS,E);
             if OBJECTTYPE='GENBANK        ' then WRITEGB(OUTFILE,E)
             else begin
                  MODEL(E,TRANSCRIPTS,SITESET,LL,LLEN,OBJ);
                  WRITEOBJECT(OUTFILE,OBJ)
                  end;
             REINIT(E);
(*!!!        if FILES then CLOSE(OUTFILE) *)
             end (* if FOUND *)
           end (* not eof(NAMEFILE) *)
         end; (* else *)

(*!!!    if not FILES then CLOSE(OUTFILE) *)
    end. (* GETOB *)
