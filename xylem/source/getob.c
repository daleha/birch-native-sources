/* Output from p2c 1.21alpha-07.Dec.93, the Pascal-to-C translator */
/* From input file "getob.p" */


/* Modification history:
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
                    upper case before searching an index file. */

/***********************************************************/
/*                                                         */
/*  GETOB     VERSION  1.2.9    8/14/2001  Standard Pascal */
/*            Brian Fristensky                             */
/*            Dept. of Plant Science                       */
/*            University of Manitoba                       */
/*            Winnipeg, MB R3T 2N2  CANADA                 */
/*                                                         */
/* SYNOPSIS                                                */
/* getob infile namefile anofile seqfile indfile message   */
/*       outfile [expfile]                                 */
/*                                                         */
/* DESCRIPTION                                             */
/*  Gets objects from GenBank file and writes to OUTFILE   */
/*  Conforms to "The DDBJ/EMBL/GenBank Feature Table:      */
/*     Definition" version 1.06 Feb. 1, 1994.              */
/*                                                         */
/*    -f        write each entry to a separate file.       */
/*              the filename consists of the locus name,   */
/*             followed by a 3-letter file extesion .obj   */
/*    -r        resolve indirect references from NAMEFILE  */
/*              into objects. Indirect references take the */
/*              form:                                      */
/*            @[<database>::]<accession>|<locus>:<location>*/
/*              -r automatically sets -c.                  */
/*              -r and -e are mutually exclusive.          */
/*    -c        NAMEFILE contains ACCESSION numbers, rather*/
/*              than locus names                           */
/*    -n        By default, the qualifier 'codon_start'    */
/*              is used to determine how many n's, if nec- */
/*              essary, must be added to the 5' end of CDS,*/
/*              mat_peptide, or sig_peptide, to preserve   */
/*              the reading frame. To turn OFF this feature*/
/*              -n must be set. -n must be set for GenBank */
/*              Release 67.0 or earlier.                   */
/*                                                         */
/*    INFILE  - instructions for what data to pull         */
/*    NAMEFILE- contains names of entries to get           */
/*    ANOFILE - contains annotation parts of entries       */
/*    SEQFILE - contains sequence parts of entries         */
/*    INDFILE - contains linenumbers for beginning of each */
/*              entry in ANOFILE and SEQFILE               */
/*    MESSAGE - lists objects processed in each entry      */
/*    OUTFILE - objects retrieved, if not -f               */
/*    EXPFILE - contains feature expressions evaluated     */
/*                                                         */
/*  Copyright (c) 1989-1996  by Brian Fristensky           */
/*  !!! in comment indicates feature which may need change */
/***********************************************************/

#include <p2c.h>


/*!!!  Some Pascals require file parameters in program heading */

/******************************************************************/
/*                     const  definition                          */
/******************************************************************/

#define MAXWORD         35   /*Max. length of WORD.STR */
#define MAXFK           15   /*Length of FK, Feature Key string */
#define MAXLINE         132   /*Max. length of GenBank entry line */
#define MAXLLEN         1000   /*Max. chars. in a location */

#define MAXSEQ          750000L   /*Max. sequence in an entry */

#define MAXOBJ          1000   /*Max. objects in an entry */
#define CRUNCHOFFSET    33   /* ASCII char. used by splitdb for -c option. */

#define DEBUG           false
    /*=true, debugging messages are printed to output */

/* BEGIN MODULE STARTARGNUM */

#define STARTARGNUM     1
    /* SUN Pascal: ARG(1) is 1st command line argument*/
/*STARTARGNUM=2;       ATT Pascal: ARG(2) is 1st command line argument*/
/* END MODULE STARTARGNUM         VERSION= 'SUNMODS     Version  8/ 9/94'; */

#define VERSION         "GETOB          Version 1.2.9  14 Aug 2001"


/******************************************************************/
/*                     type   definition                          */
/******************************************************************/

/* BEGIN MODULE TYPE.WORD */
/*   <word>::= <non-blank char>[<non-blank char>] */

typedef struct WORD {
  long LEN;
  Char STR[MAXWORD];
} WORD;

/* END MODULE TYPE.WORD         VERSION= 'SUNMODS     Version  8/ 9/94'; */

/* BEGIN MODULE TYPE.LINE */
typedef Char CHARARRAY[MAXLINE];

typedef struct LINE {
  CHARARRAY STR;
  long LEN;
} LINE;

/* END MODULE TYPE.LINE         VERSION= 'SUNMODS     Version  8/ 9/94'; */

/* -  - - - Locations of entries in database files - - - - - - */

typedef struct LOC {
  WORD NAME;
  long ANOLINE, SEQLINE;
} LOC;

/* - - - - - - - - -  NUCLEOTIDE types - - - - - - - - -  */
typedef enum {
  T, C, A, G, R, Y, M, W, S, K, D, H, V, B, N
} NUCLEOTIDE;
typedef NUCLEOTIDE SEQUENCE[MAXSEQ];
typedef Char CAR[15];
typedef NUCLEOTIDE NAR[15];



/* - - - - - - - - - - - Fragment types - - - - - - - - - */
/*- used as building blocks to create a model of an object. */

typedef Char LOCLINE[MAXLLEN];
/* string that can be parsed
                                        into a location */

typedef enum {
  position, baserange, betweenposition, featurename
} LOCTYPE;

typedef struct BASEPOSITION {
  long LNUM, RNUM;   /* if LNUM<>RNUM, this is a two_base_bound */
  Char BOUNDFLAG;   /* BASEPOSITION */
} BASEPOSITION;

typedef struct FRAGMENT {
  long REPEATS;   /* used by POLY(), #times to repeat an expr. */
  LOCTYPE LT;
  union {
    struct {
      BASEPOSITION START, FINISH;
      Char STRAND;
      long SIZE;
    } U0;
    struct {
      LOCLINE EXLOC;
      long EXLOCLEN;
    } U3;
  } UU;
} FRAGMENT;   /* FRAGMENT */

typedef struct LIST {
  struct NODE *HEAD, *TAIL;
} LIST;

/* - - - - - - - - - - - FEATURES table. - - - - - - - - - -*/

/* KEY IDentifiers for FEATURES table */
typedef enum {
  BLANK, ALLELE, ATTENUATOR, BINDING, CAATSIGNAL, CDS, CHROMOSOME, CONFLICT,
  CONTIG, CREGION, DLOOP, DREGION, DSEGMENT, ENHANCER, EXON, GCSIGNAL, HYPHEN,
  IDNA, INTRON, JREGION, JSEGMENT, LTR, MATPEPTIDE, MISCBINDING,
  MISCDIFFERENCE, MISCFEATURE, MISCRECOMB, MISCRNA, MISCSIGNAL, MISCSTRUCTURE,
  MODIFIEDBASE, MRNA, MUTATION, NREGION, OLDSEQUENCE, POLYASIGNAL, POLYASITE,
  PRECURSORRNA, PRIMTRANSCRIPT, PRIMERBIND, PROMOTER, PROTEINBIND, RBS,
  REPEATREGION, REPEATUNIT, REPORIGIN, RRNA, SATELLITE, SCRNA, SIGPEPTIDE,
  SNRNA, SOURCE, STEMLOOP, STS, SREGION, TATASIGNAL, TERMINATOR,
  TRANSITPEPTIDE, TRNA, UNSURE, VARIATION, VIRION, VREGION, VSEGMENT, TPCLIP,
  TPUTR, FPCLIP, FPUTR, M10SIGNAL, M35SIGNAL, OTHER
/* p2c: getob.p, line 272: 
 * Note: Line breaker spent 0.0 seconds, 5000 tries on line 227 [251] */
} FEATUREKEY;

typedef long KEYSET[4];


typedef Char FK[MAXFK];   /* FEATURE key string */

typedef struct FEATYPE {
  FEATUREKEY KEY;
  FK KEYSTRING;
  /* LOCATION OR QUALIFIER */
  LINE LOCQUAL;
} FEATYPE;

/* - - - - - - - - Linked Lists - - - - - - - - - - - - */
/* General purpose linked list structure for GenBank entries and
   also for sequence models. */
typedef enum {
  WNODE, LNODE, FEANODE, FRANODE
} NODETYPES;

typedef struct NODE {
  struct NODE *PREV, *NEXT;
  NODETYPES NTYPE;
  union {
    WORD WOR;
    LINE LIN;
    FEATYPE FEA;
    FRAGMENT FRA;
  } UU;
} NODE;

/* - - - - - - - - - - GenBank Entry - - - - - - - - - -  */

typedef struct ENTRY {
  WORD NAME;
  long STARTNUM, SEQLEN;
  enum {
    ss, ds, ms
  } STRANDEDNESS;
  WORD SEQTYPE;
  boolean CIRCULAR;
  WORD ENTRYSTATUS, ENTRYDATE;
  LIST DEFINITION, ACCESSION;
  WORD NID;
  long SEGNUMBER, SEGTOTAL;
  LIST FEATURETABLE;
  long ACOMP, CCOMP, GCOMP, TCOMP;   /* ENTRY */
} ENTRY;

/* - - - - - - - - - - - - OBJECTS - - - - - - - - - - - - - */
/* Model of sequence, consisting of a linked list of fragments which
   together, make up a feature. Up to MAXOBJ objects can be evaluated
   for a given entry. */

typedef struct OBJECT {
  WORD TITLE;
  boolean NOLABEL;   /*=true if feature didn't have a label */
  WORD FEALAB;
  NODE *HEAD, *TAIL;
} OBJECT;

typedef OBJECT OBJECTS[MAXOBJ + 1];

/******************************************************************/
/*                       var  definition                          */
/******************************************************************/


/* File variables */
Static _TEXT INFILE, NAMEFILE, ANOFILE, SEQFILE, INDFILE, MESSAGE, OUTFILE,
	     EXPFILE;
Static LINE FILENAME;
Static long ARGNUM;

/* Option variables */
Static boolean FILES;   /* -f: =true, each entry written to separate file */
Static boolean RESOLVEXP;
    /* -r: =true, references in NAMEFILE resolved into objects*/
Static boolean ACNO;   /* -c: =true, NAMEFILE contains ACCESSION numbers */
Static boolean PAD5PRIME;   /* -n: =true, if -n is not set */

/* Variables associated with NAMEFILE and INDFILE */
Static WORD SEQNAME;
Static LINE INDIRREF;   /* Holds indirect reference from NAMEFILE */
Static LOC LOCATIONS;   /* locations of a locus in .ANO and .SEQ files */
Static boolean FOUND;   /* =true if entry has been found in INDFILE */
Static boolean OKAY;   /* =true if INDIRREF parsed into legal <location> */

/* Global type conversion arrays. */
Static CAR NUCHAR;
Static NAR INUC, COMP;
Static FK KEYSTR[71];

/* Global data structures */
Static ENTRY E;
Static SEQUENCE SEQ;
Static long SEQLIM;   /* SEQ[SEQLIM+1..MAXSEQ] is unused as of yet */
Static long NUMOBJ;   /* number of objects in OBJ */
Static OBJECTS OBJ;   /* set of objects found in an entry */
Static FK OBJECTTYPE;
Static WORD ACCESSION;
Static LOCLINE LL;
Static long LLEN;
Static KEYSET TRANSCRIPTS, SITESET;
/*tell which feature keys are to be
                               searched for in each entry */
Static NODE *FREENODE;   /* beginning of free node list */

/* Miscellaneous variables for indexing, reading etc. */
Static WORD CIRC;   /* holds 'circular' as a WORD */
Static long CA, CS;   /* current lines in annotation and sequence files*/
Static long I;
Static Char CH, CRUNCHFLAG;


/******************************************************************/
/*              procedure and function  definition                */
/******************************************************************/

/* BEGIN MODULE READLINE */
/* Read a line from a file, omitting trailing blanks */
Static Void READLINE(F, L)
_TEXT *F;
LINE *L;
{
  long LASTNONBLANK = 0;
  long I;
  Char CH;

  /* with L*/
  L->LEN = 0;
  while (!P_eoln(F->f)) {
    CH = getc(F->f);
    if (CH == '\n')
      CH = ' ';
    if (L->LEN < MAXLINE) {
      L->LEN++;
      L->STR[L->LEN - 1] = CH;
      if (CH != ' ')
	LASTNONBLANK = L->LEN;
    }
  }
  if (!BUFEOF(F->f)) {
    fscanf(F->f, "%*[^\n]");
    getc(F->f);
  }
  L->LEN = LASTNONBLANK;
  for (I = L->LEN; I < MAXLINE; I++)
    L->STR[I] = ' ';
}  /* READLINE */


/* END MODULE READLINE         VERSION= 'SUNMODS     Version  8/ 9/94'; */

/* BEGIN MODULE WRITELINE */
/*  Write a line to a file using L char, left-justified.  */
Static Void WRITELINE(F, W_, L)
_TEXT *F;
LINE W_;
long L;
{
  long I;

  for (I = 1; I <= L; I++) {
    if (I <= W_.LEN)
      putc(W_.STR[I-1], F->f);
    else
      putc(' ', F->f);
  }
}  /* WRITELINE */


/* END MODULE WRITELINE         VERSION= 'SUNMODS     Version  8/ 9/94'; */

/* BEGIN MODULE READWORD */
/*  Read a word from a textfile           */
Static Void READWORD(F, W_)
_TEXT *F;
WORD *W_;
{
  long I;
  Char CH;

  W_->LEN = 0;
  while (P_peek(F->f) == ' ') {
    if (!P_eoln(F->f)) {
      CH = getc(F->f);
      if (CH == '\n')
	CH = ' ';
    } else if (!BUFEOF(F->f)) {
      fscanf(F->f, "%*[^\n]");
      getc(F->f);
    }
  }
  while (P_peek(F->f) != ' ') {
    if (W_->LEN < MAXWORD) {
      W_->LEN++;
      W_->STR[W_->LEN - 1] = getc(F->f);
      if (W_->STR[W_->LEN - 1] == '\n')
	W_->STR[W_->LEN - 1] = ' ';
    } else {
      CH = getc(F->f);
      if (CH == '\n')
	CH = ' ';
    }
  }
  for (I = W_->LEN; I < MAXWORD; I++)
    W_->STR[I] = ' ';
}  /* READWORD */


/* END MODULE READWORD         VERSION= 'SUNMODS     Version  8/ 9/94'; */

/* BEGIN MODULE SAMEWORD */
/* Compare two WORDS for equality */
Static boolean SAMEWORD(W1, W2)
WORD *W1, *W2;
{
  long I;
  boolean T_;

  if (W1->LEN == W2->LEN) {
    T_ = true;
    I = 1;
    while (I <= W1->LEN && T_) {
      if (W1->STR[I-1] == W2->STR[I-1])
	I++;
      else
	T_ = false;
    }
    return T_;
  } else
    return false;
}  /* SAMEWORD */


/* END MODULE SAMEWORD         VERSION= 'SUNMODS     Version  8/ 9/94'; */

/* BEGIN MODULE WRITEWORD */
/*  Write a word to a file using L char, left-justified.  */
Static Void WRITEWORD(F, W_, L)
_TEXT *F;
WORD W_;
long L;
{
  long I;

  for (I = 1; I <= L; I++) {
    if (I <= W_.LEN)
      putc(W_.STR[I-1], F->f);
    else
      putc(' ', F->f);
  }
}


/* END MODULE WRITEWORD         VERSION= 'SUNMODS     Version  8/ 9/94'; */

/* BEGIN MODULE TOUPPER */
/* Change a character from lower to uppercase */
Static Char TOUPPER(CH)
Char CH;
{
  Char Result;

  if (!islower(CH))
    return CH;
  switch (CH) {

  case 'a':
    Result = 'A';
    break;

  case 'b':
    Result = 'B';
    break;

  case 'c':
    Result = 'C';
    break;

  case 'd':
    Result = 'D';
    break;

  case 'e':
    Result = 'E';
    break;

  case 'f':
    Result = 'F';
    break;

  case 'g':
    Result = 'G';
    break;

  case 'h':
    Result = 'H';
    break;

  case 'i':
    Result = 'I';
    break;

  case 'j':
    Result = 'J';
    break;

  case 'k':
    Result = 'K';
    break;

  case 'l':
    Result = 'L';
    break;

  case 'm':
    Result = 'M';
    break;

  case 'n':
    Result = 'N';
    break;

  case 'o':
    Result = 'O';
    break;

  case 'p':
    Result = 'P';
    break;

  case 'q':
    Result = 'Q';
    break;

  case 'r':
    Result = 'R';
    break;

  case 's':
    Result = 'S';
    break;

  case 't':
    Result = 'T';
    break;

  case 'u':
    Result = 'U';
    break;

  case 'v':
    Result = 'V';
    break;

  case 'w':
    Result = 'W';
    break;

  case 'x':
    Result = 'X';
    break;

  case 'y':
    Result = 'Y';
    break;

  case 'z':
    Result = 'Z';
    break;
  }
  return Result;
}  /* TOUPPER */


/* END MODULE TOUPPER         VERSION= 'SUNMODS     Version  8/ 9/94'; */

/*****************************************************/
/* Keyword procedures.                               */
/*****************************************************/
Static Void READFK(F, KEYWORD)
_TEXT *F;
Char *KEYWORD;
{
  long I;

  for (I = 0; I < MAXFK; I++) {
    if (!P_eoln(F->f)) {
      KEYWORD[I] = getc(F->f);
      if (KEYWORD[I] == '\n')
	KEYWORD[I] = ' ';
    } else
      KEYWORD[I] = ' ';
  }
}  /* READFK */


/* Convert an input string to a keyword from FEATURETABLE.  If
   string can't be identified, assign it to blank and write message
   to MESSAGE file. */
Static Void STR2FK(KEYSTRING, NEWKEY)
Char *KEYSTRING;
FEATUREKEY *NEWKEY;
{
  *NEWKEY = BLANK;
  while (strncmp(KEYSTRING, KEYSTR[(long)(*NEWKEY)], sizeof(FK)) &&
	 (long)(*NEWKEY) < (long)OTHER)
    *NEWKEY = (FEATUREKEY)((long)(*NEWKEY) + 1);
}  /* STR2FK */


/* Create the dummy nodes HEAD & TAIL, which delimit a list */
Local Void INITLIST(HEAD, TAIL)
NODE **HEAD, **TAIL;
{
  NODE *WITH;

  *HEAD = (NODE *)Malloc(sizeof(NODE));
  *TAIL = (NODE *)Malloc(sizeof(NODE));
  WITH = *HEAD;
  WITH->NEXT = *TAIL;
  WITH->PREV = NULL;
  WITH = *TAIL;
  WITH->NEXT = NULL;
  WITH->PREV = *HEAD;
}  /* INITLIST */


/*****************************************************/
/* Initialization procedures.                        */
/*****************************************************/

/* Initialize the variables.                         */
Static Void INITIALIZE()
{
  long I;
  NODE *WITH;

  /*   Initialize INUC, NUCHAR and COMP, the arrays which */
  /*   hold the  values of the NUCLEOTIDEs.      */
  INUC[(long)T] = T;
  NUCHAR[(long)T] = 't';
  COMP[(long)T] = A;
  INUC[(long)C] = C;
  NUCHAR[(long)C] = 'c';
  COMP[(long)C] = G;
  INUC[(long)A] = A;
  NUCHAR[(long)A] = 'a';
  COMP[(long)A] = T;
  INUC[(long)G] = G;
  NUCHAR[(long)G] = 'g';
  COMP[(long)G] = C;
  INUC[(long)R] = R;
  NUCHAR[(long)R] = 'r';
  COMP[(long)R] = Y;
  INUC[(long)D] = D;
  NUCHAR[(long)D] = 'd';
  COMP[(long)D] = H;
  INUC[(long)V] = V;
  NUCHAR[(long)V] = 'v';
  COMP[(long)V] = B;
  INUC[(long)M] = M;
  NUCHAR[(long)M] = 'm';
  COMP[(long)M] = K;
  INUC[(long)K] = K;
  NUCHAR[(long)K] = 'k';
  COMP[(long)K] = M;
  INUC[(long)B] = B;
  NUCHAR[(long)B] = 'b';
  COMP[(long)B] = V;
  INUC[(long)H] = H;
  NUCHAR[(long)H] = 'h';
  COMP[(long)H] = D;
  INUC[(long)Y] = Y;
  NUCHAR[(long)Y] = 'y';
  COMP[(long)Y] = R;
  INUC[(long)W] = W;
  NUCHAR[(long)W] = 'w';
  COMP[(long)W] = W;
  INUC[(long)S] = S;
  NUCHAR[(long)S] = 's';
  COMP[(long)S] = S;
  INUC[(long)N] = N;
  NUCHAR[(long)N] = 'n';
  COMP[(long)N] = N;

  /* Initialize keyword array */
  memcpy(KEYSTR[(long)BLANK], "               ", sizeof(FK));
  memcpy(KEYSTR[(long)ALLELE], "allele         ", sizeof(FK));
  memcpy(KEYSTR[(long)ATTENUATOR], "attenuator     ", sizeof(FK));
  memcpy(KEYSTR[(long)BINDING], "binding        ", sizeof(FK));
  memcpy(KEYSTR[(long)CAATSIGNAL], "CAAT_signal    ", sizeof(FK));
  memcpy(KEYSTR[(long)CDS], "CDS            ", sizeof(FK));
  memcpy(KEYSTR[(long)CHROMOSOME], "chromosome     ", sizeof(FK));
  memcpy(KEYSTR[(long)CONFLICT], "conflict       ", sizeof(FK));
  memcpy(KEYSTR[(long)CONTIG], "contig         ", sizeof(FK));
  memcpy(KEYSTR[(long)CREGION], "C_region       ", sizeof(FK));
  memcpy(KEYSTR[(long)DLOOP], "D-loop         ", sizeof(FK));
  memcpy(KEYSTR[(long)DREGION], "D_region       ", sizeof(FK));
  memcpy(KEYSTR[(long)DSEGMENT], "D_segment      ", sizeof(FK));
  memcpy(KEYSTR[(long)ENHANCER], "enhancer       ", sizeof(FK));
  memcpy(KEYSTR[(long)EXON], "exon           ", sizeof(FK));
  memcpy(KEYSTR[(long)GCSIGNAL], "GC_signal      ", sizeof(FK));
  memcpy(KEYSTR[(long)IDNA], "iDNA           ", sizeof(FK));
  memcpy(KEYSTR[(long)INTRON], "intron         ", sizeof(FK));
  memcpy(KEYSTR[(long)JREGION], "J_region       ", sizeof(FK));
  memcpy(KEYSTR[(long)JSEGMENT], "J_segment      ", sizeof(FK));
  memcpy(KEYSTR[(long)LTR], "LTR            ", sizeof(FK));
  memcpy(KEYSTR[(long)MATPEPTIDE], "mat_peptide    ", sizeof(FK));
  memcpy(KEYSTR[(long)MISCBINDING], "misc_binding   ", sizeof(FK));
  memcpy(KEYSTR[(long)MISCDIFFERENCE], "misc_difference", sizeof(FK));
  memcpy(KEYSTR[(long)MISCFEATURE], "misc_feature   ", sizeof(FK));
  memcpy(KEYSTR[(long)MISCRECOMB], "misc_recomb    ", sizeof(FK));
  memcpy(KEYSTR[(long)MISCRNA], "misc_RNA       ", sizeof(FK));
  memcpy(KEYSTR[(long)MISCSIGNAL], "misc_signal    ", sizeof(FK));
  memcpy(KEYSTR[(long)MISCSTRUCTURE], "misc_structure ", sizeof(FK));
  memcpy(KEYSTR[(long)MODIFIEDBASE], "modified_base  ", sizeof(FK));
  memcpy(KEYSTR[(long)MRNA], "mRNA           ", sizeof(FK));
  memcpy(KEYSTR[(long)MUTATION], "mutation       ", sizeof(FK));
  memcpy(KEYSTR[(long)NREGION], "N_region       ", sizeof(FK));
  memcpy(KEYSTR[(long)OLDSEQUENCE], "old_sequence   ", sizeof(FK));
  memcpy(KEYSTR[(long)POLYASIGNAL], "polyA_signal   ", sizeof(FK));
  memcpy(KEYSTR[(long)POLYASITE], "polyA_site     ", sizeof(FK));
  memcpy(KEYSTR[(long)PRECURSORRNA], "precursor_RNA  ", sizeof(FK));
  memcpy(KEYSTR[(long)PRIMERBIND], "primer_bind    ", sizeof(FK));
  memcpy(KEYSTR[(long)PRIMTRANSCRIPT], "prim_transcript", sizeof(FK));
  memcpy(KEYSTR[(long)PROMOTER], "promoter       ", sizeof(FK));
  memcpy(KEYSTR[(long)PROTEINBIND], "protein_bind   ", sizeof(FK));
  memcpy(KEYSTR[(long)RBS], "RBS            ", sizeof(FK));
  memcpy(KEYSTR[(long)REPEATREGION], "repeat_region  ", sizeof(FK));
  memcpy(KEYSTR[(long)REPEATUNIT], "repeat_unit    ", sizeof(FK));
  memcpy(KEYSTR[(long)REPORIGIN], "rep_origin     ", sizeof(FK));
  memcpy(KEYSTR[(long)RRNA], "rRNA           ", sizeof(FK));
  memcpy(KEYSTR[(long)SATELLITE], "satellite      ", sizeof(FK));
  memcpy(KEYSTR[(long)SCRNA], "scRNA          ", sizeof(FK));
  memcpy(KEYSTR[(long)SIGPEPTIDE], "sig_peptide    ", sizeof(FK));
  memcpy(KEYSTR[(long)SNRNA], "snRNA          ", sizeof(FK));
  memcpy(KEYSTR[(long)SOURCE], "source         ", sizeof(FK));
  memcpy(KEYSTR[(long)SREGION], "S_region       ", sizeof(FK));
  memcpy(KEYSTR[(long)STEMLOOP], "stem_loop      ", sizeof(FK));
  memcpy(KEYSTR[(long)STS], "STS            ", sizeof(FK));
  memcpy(KEYSTR[(long)TATASIGNAL], "TATA_signal    ", sizeof(FK));
  memcpy(KEYSTR[(long)TERMINATOR], "terminator     ", sizeof(FK));
  memcpy(KEYSTR[(long)TRANSITPEPTIDE], "transit_peptide", sizeof(FK));
  memcpy(KEYSTR[(long)TRNA], "tRNA           ", sizeof(FK));
  memcpy(KEYSTR[(long)UNSURE], "unsure         ", sizeof(FK));
  memcpy(KEYSTR[(long)VARIATION], "variation      ", sizeof(FK));
  memcpy(KEYSTR[(long)VIRION], "virion         ", sizeof(FK));
  memcpy(KEYSTR[(long)VREGION], "V_region       ", sizeof(FK));
  memcpy(KEYSTR[(long)VSEGMENT], "V_segment      ", sizeof(FK));
  memcpy(KEYSTR[(long)TPCLIP], "3'clip         ", sizeof(FK));
  memcpy(KEYSTR[(long)TPUTR], "3'UTR          ", sizeof(FK));
  memcpy(KEYSTR[(long)FPCLIP], "5'clip         ", sizeof(FK));
  memcpy(KEYSTR[(long)FPUTR], "5'UTR          ", sizeof(FK));
  memcpy(KEYSTR[(long)M10SIGNAL], "-10_signal     ", sizeof(FK));
  memcpy(KEYSTR[(long)M35SIGNAL], "-35_signal     ", sizeof(FK));
  memcpy(KEYSTR[(long)HYPHEN], "-              ", sizeof(FK));
  memcpy(KEYSTR[(long)OTHER], "               ", sizeof(FK));

  /* E */
  CRUNCHFLAG = (Char)CRUNCHOFFSET;   /*used in RGB */
  E.NAME.LEN = 0;
  E.SEQLEN = 0;
  E.SEQTYPE.LEN = 0;
  E.ENTRYSTATUS.LEN = 0;
  E.ENTRYDATE.LEN = 0;
  INITLIST(&E.DEFINITION.HEAD, &E.DEFINITION.TAIL);
  INITLIST(&E.ACCESSION.HEAD, &E.ACCESSION.TAIL);
  E.SEGNUMBER = 0;
  E.SEGTOTAL = 0;
  E.NID.LEN = 0;
  INITLIST(&E.FEATURETABLE.HEAD, &E.FEATURETABLE.TAIL);
  WITH = E.FEATURETABLE.TAIL;
  /* BLANK indicates end of features table */
  WITH->NTYPE = FEANODE;
  WITH->UU.FEA.KEY = BLANK;
  E.ACOMP = 0;
  E.CCOMP = 0;
  E.GCOMP = 0;
  E.TCOMP = 0;
  for (I = 0; I <= MAXOBJ; I++)
    INITLIST(&OBJ[I].HEAD, &OBJ[I].TAIL);
  FREENODE = NULL;

  /* WORDS used to parse input */

  CIRC.STR[0] = 'c';
  CIRC.STR[1] = 'i';
  CIRC.STR[2] = 'r';
  CIRC.STR[3] = 'c';
  CIRC.STR[4] = 'u';
  CIRC.STR[5] = 'l';
  CIRC.STR[6] = 'a';
  CIRC.STR[7] = 'r';
  CIRC.LEN = 8;
}  /* INITIALIZE */


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Read options from command line and set FILES & RESOLVEXP. */
Static Void READOPTIONS(ARGNUM, FILES, RESOLVEXP, ACNO)
long *ARGNUM;
boolean *FILES, *RESOLVEXP, *ACNO;
{
  Char ARGUMENT[132];

  /* Set defaults */
  *FILES = false;
  *RESOLVEXP = false;
  *ACNO = false;
  PAD5PRIME = true;
  /* Read options */
  P_sun_argv(ARGUMENT, 132, (int)(*ARGNUM));
  while (ARGUMENT[0] == '-') {
    if (ARGUMENT[1] == 'n' || ARGUMENT[1] == 'c' || ARGUMENT[1] == 'r' ||
	ARGUMENT[1] == 'f') {
      switch (ARGUMENT[1]) {

      case 'f':
	*FILES = true;
	break;

      /* indirect references always contain accession numbers */
      case 'r':
	*RESOLVEXP = true;
	*ACNO = true;
	break;

      case 'c':
	*ACNO = true;
	break;

      case 'n':
	PAD5PRIME = false;
	break;
      }
    }
    (*ARGNUM)++;
    P_sun_argv(ARGUMENT, 132, (int)(*ARGNUM));
  }
}  /* READOPTIONS */


Local Void FILEARGS(F, FTYPE, ARGNUM)
_TEXT *F;
Char FTYPE;
long *ARGNUM;
{
  Char ARGUMENT[132];

  P_sun_argv(ARGUMENT, 132, (int)(*ARGNUM));
  if (FTYPE == 'I') {
    strcpy(F->name, P_trimname(ARGUMENT, 132));
    if (F->f != NULL)
      F->f = freopen(F->name, "r", F->f);
    else
      F->f = fopen(F->name, "r");
    if (F->f == NULL)
      _EscIO2(FileNotFound, F->name);
    RESETBUF(F->f, Char);
  } else {
    strcpy(F->name, P_trimname(ARGUMENT, 132));
    if (F->f != NULL)
      F->f = freopen(F->name, "w", F->f);
    else
      F->f = fopen(F->name, "w");
    if (F->f == NULL)
      _EscIO2(FileNotFound, F->name);
    SETUPBUF(F->f, Char);
  }
  (*ARGNUM)++;
}  /* FILEARGS */


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Open files.*/
Static Void OPENFILES(ARGNUM, INFILE, NAMEFILE, ANOFILE, SEQFILE, INDFILE,
		      MESSAGE, OUTFILE, EXPFILE)
long *ARGNUM;
_TEXT *INFILE, *NAMEFILE, *ANOFILE, *SEQFILE, *INDFILE, *MESSAGE, *OUTFILE,
      *EXPFILE;
{

  /* BEGIN MODULE FILEARGS */
  /* This procedure overcomes one of the stupidest aspects of UNIX Pascal,
     namely the fact that filenames in the program statement are supposed to
     be actual UNIX filenames!  To overcome this, the 2-argument version of
     reset and rewrite must be used with string variables.  This module
     need only contain the reset and rewrite statements in any normal
     implementation of Pascal. */
  /* END MODULE FILEARGS         VERSION= 'SUNMODS     Version  8/ 9/94'; */

  FILEARGS(INFILE, 'I', ARGNUM);
  FILEARGS(NAMEFILE, 'I', ARGNUM);
  FILEARGS(ANOFILE, 'I', ARGNUM);
  FILEARGS(SEQFILE, 'I', ARGNUM);
  FILEARGS(INDFILE, 'I', ARGNUM);
  FILEARGS(MESSAGE, 'O', ARGNUM);
  if (!FILES)
    FILEARGS(OUTFILE, 'O', ARGNUM);
  if (!RESOLVEXP)
    FILEARGS(EXPFILE, 'O', ARGNUM);
}  /* OPENFILES */


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Read input file containing instructions for processing each entry.*/
Static Void READINSTRUCTIONS(INFILE, TRANSCRIPTS, SITESET, OBJECTTYPE)
_TEXT *INFILE;
long *TRANSCRIPTS, *SITESET;
Char *OBJECTTYPE;
{
  FK OBJID;
  FEATUREKEY OBJKEY;
  KEYSET SET;

  P_expset(TRANSCRIPTS, 0L);
  P_expset(SITESET, 0L);
  READFK(INFILE, OBJECTTYPE);
  fscanf(INFILE->f, "%*[^\n]");
  getc(INFILE->f);
  if (!strncmp(OBJECTTYPE, "GENBANK        ", sizeof(FK))) {
    return;
  }  /* OBJECTTYPE<> 'GENBANK        ' */
  /*!!!   CLOSE(INFILE); */
  READFK(INFILE, OBJID);
  fscanf(INFILE->f, "%*[^\n]");
  getc(INFILE->f);
  while ((strncmp(OBJID, "SITES          ", sizeof(FK)) != 0) &
	 (!BUFEOF(INFILE->f))) {
    STR2FK(OBJID, &OBJKEY);
    if (OBJKEY != OTHER)
      P_addset(TRANSCRIPTS, (int)OBJKEY);
    if (!BUFEOF(INFILE->f)) {
      READFK(INFILE, OBJID);
      fscanf(INFILE->f, "%*[^\n]");
      getc(INFILE->f);
    }
  }  /* OBJID <> 'SITES         ' */
  while (!BUFEOF(INFILE->f)) {
    STR2FK(OBJID, &OBJKEY);
    if (OBJKEY != OTHER)
      P_addset(SITESET, (int)OBJKEY);
    if (BUFEOF(INFILE->f))
      break;
    READFK(INFILE, OBJID);
    fscanf(INFILE->f, "%*[^\n]");
    getc(INFILE->f);
  }  /* not eof(INFILE) */
  /*!!!*/
}  /* READINSTRUCTIONS */


/**************************************************************/
/* Read a LOCUS name for next sequence to pull from database. */
/**************************************************************/
Static Void READNAME(NAMEFILE, SEQNAME)
_TEXT *NAMEFILE;
WORD *SEQNAME;
{
  /* Skip comment lines */
  while ((P_peek(NAMEFILE->f) == ';') & (!BUFEOF(NAMEFILE->f))) {
    fscanf(NAMEFILE->f, "%*[^\n]");
    getc(NAMEFILE->f);
  }

  /* Read the first token in the next non-comment line */
  SEQNAME->LEN = 0;
  if (!BUFEOF(NAMEFILE->f))
    READWORD(NAMEFILE, SEQNAME);
  if (!BUFEOF(NAMEFILE->f)) {
    fscanf(NAMEFILE->f, "%*[^\n]");
    getc(NAMEFILE->f);
  }
}  /* READNAME */


/* Local variables for FINDDATA: */
struct LOC_FINDDATA {
  _TEXT *INDFILE;
  WORD ID;
  LOC *LOCATIONS;
  boolean *FOUND;
  Char FLAG;
  WORD DUMMY;
} ;

/* Advance to first line with first char of ID */
Local Void ZOOM(LINK)
struct LOC_FINDDATA *LINK;
{
  while ((!BUFEOF(LINK->INDFILE->f)) & (P_peek(LINK->INDFILE->f) != LINK->FLAG)) {
    fscanf(LINK->INDFILE->f, "%*[^\n]");
    getc(LINK->INDFILE->f);
  }
}  /* ZOOM */

/* Search for LOCUS corresponding to ID */
Local Void LCRAWL(LINK)
struct LOC_FINDDATA *LINK;
{
  LOC *WITH;

  while (P_peek(LINK->INDFILE->f) == LINK->FLAG && !*LINK->FOUND) {
    WITH = LINK->LOCATIONS;
    READWORD(LINK->INDFILE, &WITH->NAME);
    if (SAMEWORD(&LINK->ID, &WITH->NAME)) {
      *LINK->FOUND = true;
      READWORD(LINK->INDFILE, &LINK->DUMMY);   /* ignore ACCESSION number */
      fscanf(LINK->INDFILE->f, "%ld%ld", &WITH->ANOLINE, &WITH->SEQLINE);
    }
    if (!BUFEOF(LINK->INDFILE->f)) {
      fscanf(LINK->INDFILE->f, "%*[^\n]");
      getc(LINK->INDFILE->f);
    }
  }
}  /* LCRAWL */

/* Search for ACCESSION number corresponding to ID */
Local Void ACRAWL(LINK)
struct LOC_FINDDATA *LINK;
{
  LOC *WITH;

  while ((!*LINK->FOUND) & (!BUFEOF(LINK->INDFILE->f))) {
    WITH = LINK->LOCATIONS;
    READWORD(LINK->INDFILE, &LINK->DUMMY);   /* skip over LOCUS name */
    READWORD(LINK->INDFILE, &WITH->NAME);
    if (SAMEWORD(&LINK->ID, &WITH->NAME)) {
      *LINK->FOUND = true;
      fscanf(LINK->INDFILE->f, "%ld%ld", &WITH->ANOLINE, &WITH->SEQLINE);
    }
    if (!BUFEOF(LINK->INDFILE->f)) {
      fscanf(LINK->INDFILE->f, "%*[^\n]");
      getc(LINK->INDFILE->f);
    }
  }
}  /* ACRAWL */


/**********************************************************************/
/* Look for a given LOCUS in the index file.  Since names in NAMEFILE */
/* are supposed to be sorted, reaching the end of file means that the */
/* locus in question is not in the database.                          */
/**********************************************************************/
Static Void FINDDATA(INDFILE_, ID_, LOCATIONS_, FOUND_)
_TEXT *INDFILE_;
WORD ID_;
LOC *LOCATIONS_;
boolean *FOUND_;
{
  struct LOC_FINDDATA V;
  long FORLIM;

  V.INDFILE = INDFILE_;
  V.ID = ID_;
  V.LOCATIONS = LOCATIONS_;
  V.FOUND = FOUND_;
  if (BUFEOF(V.INDFILE->f)) {
    if (*V.INDFILE->name != '\0') {
      if (V.INDFILE->f != NULL)
	V.INDFILE->f = freopen(V.INDFILE->name, "r", V.INDFILE->f);
      else
	V.INDFILE->f = fopen(V.INDFILE->name, "r");
    } else
      rewind(V.INDFILE->f);
    if (V.INDFILE->f == NULL)
      _EscIO2(FileNotFound, V.INDFILE->name);
    RESETBUF(V.INDFILE->f, Char);
  }
  *V.FOUND = false;
  FORLIM = V.ID.LEN;
  /* Make sure that the name or accession number is capitalized. */
  for (I = 1; I <= FORLIM; I++)
    V.ID.STR[I-1] = TOUPPER(V.ID.STR[I-1]);
  V.FLAG = V.ID.STR[0];

  /* Search to end of file or until FOUND */
  if (ACNO)
    ACRAWL(&V);
  else {
    while ((!*V.FOUND) & (!BUFEOF(V.INDFILE->f))) {
      ZOOM(&V);
      LCRAWL(&V);
    }
  }

  /* If list isn't sorted, go back to the beginning and search from top */
  if (*V.FOUND)
    return;
  if (*V.INDFILE->name != '\0') {
    if (V.INDFILE->f != NULL)
      V.INDFILE->f = freopen(V.INDFILE->name, "r", V.INDFILE->f);
    else
      V.INDFILE->f = fopen(V.INDFILE->name, "r");
  } else
    rewind(V.INDFILE->f);
  if (V.INDFILE->f == NULL)
    _EscIO2(FileNotFound, V.INDFILE->name);
  RESETBUF(V.INDFILE->f, Char);
  if (ACNO) {
    ACRAWL(&V);
    return;
  }
  while ((!*V.FOUND) & (!BUFEOF(V.INDFILE->f))) {
    ZOOM(&V);
    LCRAWL(&V);
  }
}  /* FINDDATA */


/******************************************/
/* Parse a reference into its components. */
/******************************************/
Static Void PARSEREF(INDIRREF, ACCESSION, LL, LLEN, OKAY)
LINE *INDIRREF;
WORD *ACCESSION;
Char *LL;
long *LLEN;
boolean *OKAY;
{  /* PARSREF */
  long I = 2, J = 1;
  long STARTACC, FINISHACC, FORLIM;

  /* with INDIRREF */
  /* Ignore DATABASE component. It is assumed that the necessary entries
    have already been obtained from the appropriate database and
    included in ANOFILE, SEQFILE, INDFILE */
  *OKAY = false;
  while (I < INDIRREF->LEN && INDIRREF->STR[I-1] != ':')
    I++;
  if (I >= INDIRREF->LEN) {  /* 1*/
    return;
  }  /* 1 */
  if (INDIRREF->STR[I] == ':') {
    STARTACC = I + 2;
    I = STARTACC + 1;
    while (I < INDIRREF->LEN && INDIRREF->STR[I-1] != ':')
      I++;
  } else
    STARTACC = 2;
  FINISHACC = I - 1;

  /* Read the accession number */
  if (I >= INDIRREF->LEN) {  /* 2 */
    return;
  }  /* 2 */
  for (I = STARTACC - 1; I < FINISHACC; I++) {
    ACCESSION->STR[J-1] = INDIRREF->STR[I];
    J++;
  }
  ACCESSION->LEN = FINISHACC - STARTACC + 1;

  /* Read the <location> expression */
  J = 1;
  FORLIM = INDIRREF->LEN;
  for (I = FINISHACC + 1; I < FORLIM; I++) {
    LL[J-1] = INDIRREF->STR[I];
    J++;
  }
  *LLEN = INDIRREF->LEN - FINISHACC - 1;
  *OKAY = true;
}  /* PARSREF */


/*****************************************************/
/* Create a filename using the locus.                */
/*****************************************************/
Static Void MAKEFN(NAME, FILENAME)
WORD NAME;
LINE *FILENAME;
{
  long I = 1;

  while (I <= NAME.LEN) {
    FILENAME->STR[I-1] = NAME.STR[I-1];
    I++;
  }
  FILENAME->STR[I-1] = '.';
  FILENAME->STR[I] = 'o';
  FILENAME->STR[I+1] = 'b';
  FILENAME->STR[I+2] = 'j';
  FILENAME->LEN = I + 3;
  for (I = FILENAME->LEN; I < MAXLINE; I++)
    FILENAME->STR[I] = ' ';
}  /* MAKEFN */


/*****************************************/
/*  Linked-list operations for node list.*/
/*****************************************/

/*Add a node after DEST*/
Static Void ADDNODE(DEST)
NODE *DEST;
{
  NODE *TEMP, *NEWNODE;

  /*Get a new node from freelist.*/
  if (FREENODE == NULL)
    NEWNODE = (NODE *)Malloc(sizeof(NODE));
  else {
    NEWNODE = FREENODE;
    FREENODE = FREENODE->NEXT;
  }
  /* Add the node */
  TEMP = DEST->NEXT;
  NEWNODE->NEXT = TEMP;
  NEWNODE->PREV = DEST;
  DEST->NEXT = NEWNODE;
  TEMP->PREV = NEWNODE;
}  /* ADDNODE */


/*Return a list to the top of freelist*/
Static Void RIDOF(HEAD, TAIL)
NODE **HEAD, **TAIL;
{
  NODE *FIRST, *LAST;

  if ((*HEAD)->NEXT == *TAIL)
    return;
  FIRST = (*HEAD)->NEXT;
  LAST = (*TAIL)->PREV;
  FIRST->PREV = NULL;
  LAST->NEXT = FREENODE;
  FREENODE = FIRST;
  (*HEAD)->NEXT = *TAIL;
  (*TAIL)->PREV = *HEAD;
}  /* RIDOF */


/* BEGIN MODULE NUC */
/*****************************************************************/
/*  Convert a character to the appropriate nucleotide symbol.    */
/*****************************************************************/
Static NUCLEOTIDE NUC(CH)
Char CH;
{
  NUCLEOTIDE Result;

  switch (CH) {

  case 'A':
  case 'a':
    Result = A;
    break;

  case 'C':
  case 'c':
    Result = C;
    break;

  case 'G':
  case 'g':
    Result = G;
    break;

  case 'T':
  case 't':
  case 'U':
  case 'u':
    Result = T;
    break;

  case 'R':
  case 'r':
    Result = R;
    break;

  case 'M':
  case 'm':
    Result = M;
    break;

  case 'B':
  case 'b':
    Result = B;
    break;

  case 'N':
  case 'n':
    Result = N;
    break;

  case 'Y':
  case 'y':
    Result = Y;
    break;

  case 'K':
  case 'k':
    Result = K;
    break;

  case 'D':
  case 'd':
    Result = D;
    break;

  case 'S':
  case 's':
    Result = S;
    break;

  case 'W':
  case 'w':
    Result = W;
    break;

  case 'H':
  case 'h':
    Result = H;
    break;

  case 'V':
  case 'v':
    Result = V;
    break;
  }
  return Result;
}


/* Local variables for RGB: */
struct LOC_RGB {
  _TEXT *ANOFILE;
  ENTRY *E;
  long BLANKSET[9];
} ;

Local Void SKIP(F, POS, LINK)
_TEXT *F;
long POS;
struct LOC_RGB *LINK;
{
  long I = 0;
  Char CH;

  /* Read past compressed leading blanks, if any.*/
  if (P_peek(F->f) == CRUNCHFLAG) {
    CH = getc(F->f);
    if (CH == '\n')
      CH = ' ';
    CH = getc(F->f);
    if (CH == '\n')
      CH = ' ';
    return;
  }
  /* Otherwise, read past the specified number of characters.*/
  while ((I < POS) & (!P_eoln(F->f))) {
    CH = getc(F->f);
    if (CH == '\n')
      CH = ' ';
    I++;
  }
}  /* SKIP */

/* - - - - - - - - - - - - - - - - - - - - - - - - - */
Local Void READLOCUS(LINK)
struct LOC_RGB *LINK;
{
  WORD TOKEN;
  ENTRY *WITH;
  _TEXT TEMP;

  WITH = LINK->E;
  /* with E */
  SKIP(LINK->ANOFILE, 12L, LINK);
  READWORD(LINK->ANOFILE, &WITH->NAME);
  fscanf(LINK->ANOFILE->f, "%ld", &WITH->SEQLEN);
  READWORD(LINK->ANOFILE, &TOKEN);   /* bp */
  /* After 'bp', there is an optional field telling
  the type of molecule (ss-RNA, ds-DNA etc.)
  Since this field is optional, we must test the
  next two tokens to see if they are 'circular' */
  WITH->CIRCULAR = false;
  if (!BUFEOF(LINK->ANOFILE->f))
    READWORD(LINK->ANOFILE, &TOKEN);
  if (SAMEWORD(&TOKEN, &CIRC))
    WITH->CIRCULAR = true;
  else {
    READWORD(LINK->ANOFILE, &TOKEN);
    if (SAMEWORD(&TOKEN, &CIRC))
      WITH->CIRCULAR = true;
  }
  fscanf(LINK->ANOFILE->f, "%*[^\n]");
  getc(LINK->ANOFILE->f);
  /*!!!*/
  if (DEBUG) {
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITEWORD(&TEMP, WITH->NAME, 18L);
    printf("%12ld\n", WITH->SEQLEN);
  }
  CA++;
}  /* READLOCUS */

/* - - - - - - - - - - - - - - - - - - - - - - - - - */
Local Void READDEFINITION(LINK)
struct LOC_RGB *LINK;
{
  LIST *WITH;
  NODE *WITH1;
  _TEXT TEMP;

  WITH = &LINK->E->DEFINITION;
  do {
    SKIP(LINK->ANOFILE, 12L, LINK);
    ADDNODE(WITH->TAIL->PREV);
    WITH1 = WITH->TAIL->PREV;
    /* with TAIL^ */
    WITH1->NTYPE = LNODE;
    READLINE(LINK->ANOFILE, &WITH1->UU.LIN);
    CA++;
    /*!!!*/
    if (DEBUG) {
      TEMP.f = stdout;
      *TEMP.name = '\0';
      WRITELINE(&TEMP, WITH1->UU.LIN, WITH1->UU.LIN.LEN);
      putchar('\n');
    }
  } while (P_inset(P_peek(LINK->ANOFILE->f), LINK->BLANKSET));
      /* READDEFINITION */
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - */
Local Void READACCESSION(LINK)
struct LOC_RGB *LINK;
{
  LIST *WITH;
  _TEXT TEMP;

  WITH = &LINK->E->ACCESSION;
  do {
    SKIP(LINK->ANOFILE, 12L, LINK);
    while (!P_eoln(LINK->ANOFILE->f)) {
      ADDNODE(WITH->TAIL->PREV);
      WITH->TAIL->PREV->NTYPE = WNODE;
      READWORD(LINK->ANOFILE, &WITH->TAIL->PREV->UU.WOR);
      /*!!!*/
      if (DEBUG) {
	TEMP.f = stdout;
	*TEMP.name = '\0';
	WRITEWORD(&TEMP, WITH->TAIL->PREV->UU.WOR, 10L);
      }
    }
    fscanf(LINK->ANOFILE->f, "%*[^\n]");
    getc(LINK->ANOFILE->f);
    /*!!!*/
    if (DEBUG)
      putchar('\n');
    CA++;
  } while (P_inset(P_peek(LINK->ANOFILE->f), LINK->BLANKSET));
      /* READACCESSION */
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - */
Local Void READNID(LINK)
struct LOC_RGB *LINK;
{
  /*!!!*/
  if (DEBUG)
    printf("NID\n");
  SKIP(LINK->ANOFILE, 11L, LINK);
  READWORD(LINK->ANOFILE, &LINK->E->NID);
  fscanf(LINK->ANOFILE->f, "%*[^\n]");
  getc(LINK->ANOFILE->f);
  CA++;
}  /* READNID */

/* - - - - - - - - - - - - - - - - - - - - - - - - - */
Local Void READSEGMENT(LINK)
struct LOC_RGB *LINK;
{
  WORD DUMMY;

  /*!!!*/
  if (DEBUG)
    printf("SEGMENT\n");
  SKIP(LINK->ANOFILE, 11L, LINK);
  fscanf(LINK->ANOFILE->f, "%ld", &LINK->E->SEGNUMBER);
  READWORD(LINK->ANOFILE, &DUMMY);
  fscanf(LINK->ANOFILE->f, "%ld%*[^\n]", &LINK->E->SEGTOTAL);
  getc(LINK->ANOFILE->f);
  CA++;
}  /* READSEGMENT */

/* - - - - - - - - - - - - - - - - - - - - - - - - - */
/* This procedure skips over parts of an entry that
   aren't used by GETOB. */
Local Void READDUMMY(LINK)
struct LOC_RGB *LINK;
{
  do {
    fscanf(LINK->ANOFILE->f, "%*[^\n]");
    getc(LINK->ANOFILE->f);
    CA++;
  } while (P_inset(P_peek(LINK->ANOFILE->f), LINK->BLANKSET));
      /* READDUMMY */
}

/* Local variables for READFEATURES: */
struct LOC_READFEATURES {
  struct LOC_RGB *LINK;
} ;

Local Void FEALINE(DEST, LINK)
NODE **DEST;
struct LOC_READFEATURES *LINK;
{
  Char CH;
  long NUMBLANKS;
  FEATYPE *WITH;
  _TEXT TEMP;

  ADDNODE(*DEST);
  (*DEST)->NTYPE = FEANODE;
  WITH = &(*DEST)->UU.FEA;
  /* with DEST^ */
  if (P_peek(LINK->LINK->ANOFILE->f) == CRUNCHFLAG) {
    CH = getc(LINK->LINK->ANOFILE->f);
    if (CH == '\n')
      CH = ' ';
    CH = getc(LINK->LINK->ANOFILE->f);
    if (CH == '\n')
      CH = ' ';
    NUMBLANKS = CH - CRUNCHOFFSET;
    if (NUMBLANKS == 5) {  /*next field is feature key */
      READFK(LINK->LINK->ANOFILE, WITH->KEYSTRING);
      STR2FK(WITH->KEYSTRING, &WITH->KEY);
      CH = getc(LINK->LINK->ANOFILE->f);
      if (CH == '\n')
	CH = ' ';
    } else
      WITH->KEY = BLANK;   /*next field is qualifier line */
  }  /* ANOFILE^=CRUNCHFLAG */
  else {
    SKIP(LINK->LINK->ANOFILE, 5L, LINK->LINK);
    READFK(LINK->LINK->ANOFILE, WITH->KEYSTRING);
    STR2FK(WITH->KEYSTRING, &WITH->KEY);
    CH = getc(LINK->LINK->ANOFILE->f);
    if (CH == '\n')
      CH = ' ';
  }
  READLINE(LINK->LINK->ANOFILE, &WITH->LOCQUAL);
  CA++;
  /*!!!*/
  if (!DEBUG)
    return;
  printf("%.*s   ", MAXFK, KEYSTR[(long)WITH->KEY]);
  TEMP.f = stdout;
  *TEMP.name = '\0';
  WRITELINE(&TEMP, WITH->LOCQUAL, WITH->LOCQUAL.LEN);
  putchar('\n');
}  /* FEALINE */

/* - - - - - - - - - - - - - - - - - - - - - - - - - */
Local Void READFEATURES(LINK)
struct LOC_RGB *LINK;
{
  struct LOC_READFEATURES V;
  LIST *WITH;

  V.LINK = LINK;
  /*!!!*/
  if (DEBUG)
    printf("FEATURES\n");
  fscanf(LINK->ANOFILE->f, "%*[^\n]");
  getc(LINK->ANOFILE->f);
  CA++;

  /* Read the FEATURES table */
  WITH = &LINK->E->FEATURETABLE;
  while (P_inset(P_peek(LINK->ANOFILE->f), LINK->BLANKSET))
    FEALINE(&WITH->TAIL->PREV, &V);

}  /* READFEATURES */

/* - - - - - - - - - - - - - - - - - - - - - - - - - */
Local Void READBASECOUNT(LINK)
struct LOC_RGB *LINK;
{
  WORD DUMMY;
  ENTRY *WITH;

  /*!!!*/
  if (DEBUG)
    printf("BASE COUNT\n");
  WITH = LINK->E;
  SKIP(LINK->ANOFILE, 12L, LINK);
  fscanf(LINK->ANOFILE->f, "%ld", &WITH->ACOMP);
  READWORD(LINK->ANOFILE, &DUMMY);
  fscanf(LINK->ANOFILE->f, "%ld", &WITH->CCOMP);
  READWORD(LINK->ANOFILE, &DUMMY);
  fscanf(LINK->ANOFILE->f, "%ld", &WITH->GCOMP);
  READWORD(LINK->ANOFILE, &DUMMY);
  fscanf(LINK->ANOFILE->f, "%ld", &WITH->TCOMP);
  fscanf(LINK->ANOFILE->f, "%*[^\n]");
  getc(LINK->ANOFILE->f);
  CA++;
}  /* READBASECOUNT */

/* - - - - - - - - - - - - - - - - - - - - - - - - - */
Local Void READORIGIN(LINK)
struct LOC_RGB *LINK;
{
  fscanf(LINK->ANOFILE->f, "%*[^\n]");
  getc(LINK->ANOFILE->f);
  CA++;
}  /* READORIGIN */

/* - - - - - - - - - - - - - - - - - - - - - - - - - */
Local Void READSEQUENCE(SEQFILE, LINK)
_TEXT *SEQFILE;
struct LOC_RGB *LINK;
{
  Char CH;
  boolean TOOLONG = false;
  ENTRY *WITH;

  WITH = LINK->E;
  /* with E */
  WITH->SEQLEN = 0;
  while ((!BUFEOF(SEQFILE->f)) & (P_peek(SEQFILE->f) != '>')) {
    while (!P_eoln(SEQFILE->f)) {
      CH = getc(SEQFILE->f);
      if (CH == '\n')
	CH = ' ';
      if (CH == 'b' || CH == 'B' || CH == 'v' || CH == 'V' || CH == 'h' ||
	  CH == 'H' || CH == 'd' || CH == 'D' || CH == 'k' || CH == 'K' ||
	  CH == 's' || CH == 'S' || CH == 'w' || CH == 'W' || CH == 'm' ||
	  CH == 'M' || CH == 'y' || CH == 'Y' || CH == 'r' || CH == 'R' ||
	  CH == 'n' || CH == 'N' || CH == 'u' || CH == 'U' || CH == 't' ||
	  CH == 'T' || CH == 'g' || CH == 'G' || CH == 'c' || CH == 'C' ||
	  CH == 'a' || CH == 'A') {
/* p2c: getob.p, line 1111: 
 * Note: Line breaker spent 0.0 seconds, 5000 tries on line 1668 [251] */
	if (WITH->SEQLEN < MAXSEQ - 2) {
	  WITH->SEQLEN++;
	  SEQ[WITH->SEQLEN - 1] = NUC(CH);
	} else
	  TOOLONG = true;
      }
    }  /* eoln */
    fscanf(SEQFILE->f, "%*[^\n]");
    getc(SEQFILE->f);
    CS++;
  }  /* eof */
  SEQLIM = WITH->SEQLEN;
  if (TOOLONG) {
    WRITEWORD(&MESSAGE, WITH->NAME, WITH->NAME.LEN + 1);
    fprintf(MESSAGE.f,
      ">>> WARNING! Sequence length exceeds MAXSEQ-2. Seq. truncated.\n");
  }  /* TOOLONG */
}  /* READSEQUENCE */


/* END MODULE NUC         VERSION= 'SUNMODS     Version  8/ 9/94'; */

/*******************************************************/
/* Read a GenBank entry into a structure.              */
/*******************************************************/
Static Void RGB(ANOFILE_, SEQFILE, LOCATIONS, E_)
_TEXT *ANOFILE_, *SEQFILE;
LOC LOCATIONS;
ENTRY *E_;
{  /* RGB --------------------------------------------------------*/
  struct LOC_RGB V;
  boolean DONE = false;
  Char CH;

  V.ANOFILE = ANOFILE_;
  V.E = E_;
  P_addset(P_expset(V.BLANKSET, 0L), ' ');
  P_addset(V.BLANKSET, CRUNCHFLAG);
  /* Advance to the entry specified by ANOLINE in ANOFILE */
  if (CA > LOCATIONS.ANOLINE) {
    if (*V.ANOFILE->name != '\0') {
      if (V.ANOFILE->f != NULL)
	V.ANOFILE->f = freopen(V.ANOFILE->name, "r", V.ANOFILE->f);
      else
	V.ANOFILE->f = fopen(V.ANOFILE->name, "r");
    } else
      rewind(V.ANOFILE->f);
    if (V.ANOFILE->f == NULL)
      _EscIO2(FileNotFound, V.ANOFILE->name);
    RESETBUF(V.ANOFILE->f, Char);
    CA = 1;
  }
  while (CA < LOCATIONS.ANOLINE) {
    fscanf(V.ANOFILE->f, "%*[^\n]");
    getc(V.ANOFILE->f);
    CA++;
  }

  /* Advance to the sequence specified by SEQLINE in SEQFILE */
  if (CS > LOCATIONS.SEQLINE) {
    if (*SEQFILE->name != '\0') {
      if (SEQFILE->f != NULL)
	SEQFILE->f = freopen(SEQFILE->name, "r", SEQFILE->f);
      else
	SEQFILE->f = fopen(SEQFILE->name, "r");
    } else
      rewind(SEQFILE->f);
    if (SEQFILE->f == NULL)
      _EscIO2(FileNotFound, SEQFILE->name);
    RESETBUF(SEQFILE->f, Char);
    CS = 1;
  }
  while (CS < LOCATIONS.SEQLINE) {
    fscanf(SEQFILE->f, "%*[^\n]");
    getc(SEQFILE->f);
    CS++;
  }
  fscanf(SEQFILE->f, "%*[^\n]");
  getc(SEQFILE->f);
  CS++;

  /* Read in the entry */
  while ((!DONE) & (!BUFEOF(V.ANOFILE->f))) {
    if (P_peek(V.ANOFILE->f) != '/' && P_peek(V.ANOFILE->f) != 'O' &&
	P_peek(V.ANOFILE->f) != 'B' && P_peek(V.ANOFILE->f) != 'F' &&
	P_peek(V.ANOFILE->f) != 'C' && P_peek(V.ANOFILE->f) != 'R' &&
	P_peek(V.ANOFILE->f) != 'S' && P_peek(V.ANOFILE->f) != 'K' &&
	P_peek(V.ANOFILE->f) != 'N' && P_peek(V.ANOFILE->f) != 'A' &&
	P_peek(V.ANOFILE->f) != 'D' && P_peek(V.ANOFILE->f) != 'L') {
      READDUMMY(&V);
      continue;
    }
    switch (P_peek(V.ANOFILE->f)) {

    case 'L':
      READLOCUS(&V);
      break;

    case 'D':
      READDEFINITION(&V);
      break;

    case 'A':
      READACCESSION(&V);
      break;

    case 'N':
      READNID(&V);
      break;

    case 'K':
      READDUMMY(&V);
      break;

    case 'S':
      CH = getc(V.ANOFILE->f);
      if (CH == '\n')
	CH = ' ';
      switch (P_peek(V.ANOFILE->f)) {

      case 'E':
	READSEGMENT(&V);
	break;

      case 'O':
	READDUMMY(&V);
	break;
      }
      break;

    case 'R':
      READDUMMY(&V);
      break;

    case 'C':
      READDUMMY(&V);
      break;

    case 'F':
      READFEATURES(&V);
      break;

    case 'B':
      READBASECOUNT(&V);
      break;

    case 'O':
      READORIGIN(&V);
      break;

    case '/':
      DONE = true;
      break;
    }
  }
  READSEQUENCE(SEQFILE, &V);
}  /* RGB */


Local boolean BETWEEN(X, Y_, Z)
long X, Y_, Z;
{
  if (X <= Y_ && Y_ <= Z)
    return true;
  else
    return false;
}  /* BETWEEN */


/****************************************************************/
/* Test to see if a coordinate is in a given fragment.          */
/****************************************************************/
Static boolean WITHIN(FR, COORD)
NODE *FR;
long COORD;
{
  boolean Result;
  long FIRST, LAST;
  FRAGMENT *WITH1;

  if (FR == NULL)
    return false;
  if (FR->NTYPE != FRANODE)
    return false;
  WITH1 = &FR->UU.FRA;
  /* with FRA */
  FIRST = WITH1->UU.U0.START.LNUM;
  LAST = WITH1->UU.U0.FINISH.RNUM;
  switch (WITH1->UU.U0.STRAND) {

  case ' ':
    if (FIRST <= LAST)
      Result = BETWEEN(FIRST, COORD, LAST);
    else
      Result = BETWEEN(FIRST, COORD, E.SEQLEN) | BETWEEN(1L, COORD, LAST);
    break;

  case 'C':
  case 'c':
    if (FIRST >= LAST)
      Result = BETWEEN(LAST, COORD, FIRST);
    else
      Result = BETWEEN(LAST, COORD, E.SEQLEN) | BETWEEN(1L, COORD, FIRST);
	  /* iff CIRCULAR */
    break;
  }/* case STRAND */
  return Result;
}  /* WITHIN */


#define TABSIZE         4


/* Local variables for MODEL: */
struct LOC_MODEL {
  ENTRY *E;
  OBJECT *OBJ;

  boolean OKAY;
  long INDENT;   /* # columns to indent in MESSAGE (see TAB,UNTAB)  */
  long CURRENT;
} ;

Local Void EVALUATE PP((Char *LL, long *POSN, long *LLEN, LOCTYPE *RESULT,
			BASEPOSITION *BP, NODE **DEST,
			struct LOC_MODEL *LINK));

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* To enhance the readability of the MESSAGE file, the indentation of
   output is incremented by TABSIZE columns for each expression
   encountered. Thus, at the beginning of an expression, TAB adds TABSIZE
   columns to the variable INDENT, and at the end of an expression,
   UNTAB subtracts TABSIZE columns from INDENT. */
Local Void TAB(LINK)
struct LOC_MODEL *LINK;
{
  LINK->INDENT += TABSIZE;
}

Local Void UNTAB(LINK)
struct LOC_MODEL *LINK;
{
  LINK->INDENT -= TABSIZE;
  if (LINK->INDENT < 0)
    LINK->INDENT = 0;
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Copy the location, which may span several FEATURE table lines, into
   the array LL. When DONE, TEMP will point to the line on which the
   first qualifier is found, or to the next feature, or to the end of
   the FEATURES table.*/
Local Void ONESTRING(LL, LLEN, TEMP, LINK)
Char *LL;
long *LLEN;
NODE **TEMP;
struct LOC_MODEL *LINK;
{
  boolean DONE = false;
  long I;
  LINE *WITH;
  LIST *WITH1;
  FEATYPE *WITH2;

  *LLEN = 0;
  do {
    /* Copy characters from the current location line until the end or
       until a slash (/) is reached. */
    WITH = &(*TEMP)->UU.FEA.LOCQUAL;
    I = 1;
    while (I <= WITH->LEN && WITH->STR[I-1] != '/') {
      (*LLEN)++;
      LL[*LLEN - 1] = WITH->STR[I-1];
      I++;
    }
    if (WITH->STR[I-1] == '/')
      DONE = true;
    /* If the location continues on another line, continue reading.*/
    if (!DONE || I > (*TEMP)->UU.FEA.LOCQUAL.LEN) {
      WITH1 = &LINK->E->FEATURETABLE;
      *TEMP = (*TEMP)->NEXT;
      if (*TEMP == WITH1->TAIL)
	DONE = true;
      else {
	WITH2 = &(*TEMP)->UU.FEA;
	/* TEMP^.FEA */
	if (WITH2->LOCQUAL.STR[0] == '/')
	  DONE = true;
	else if (WITH2->KEY != BLANK)
	  DONE = true;
      }
    }
  } while (!DONE);   /* ONESTRING */
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Extract characters from the qualifier line until a '='  or
   end of line is reached. */
Local Void EXTRACTQUAL(STR, LEN, STRPOSN, TESTSTR, LINK)
Char *STR;
long *LEN, *STRPOSN;
Char *TESTSTR;
struct LOC_MODEL *LINK;
{
  long I = 0;

  while (STR[*STRPOSN - 1] != '=' && STR[*STRPOSN - 1] != ':' &&
	 STR[*STRPOSN - 1] != ',' && STR[*STRPOSN - 1] != ')' &&
	 STR[*STRPOSN - 1] != '(' && *STRPOSN <= *LEN) {
    I++;
    if (I <= MAXFK)   /* use <= MAXFK chars.*/
      TESTSTR[I-1] = STR[*STRPOSN - 1];
    (*STRPOSN)++;
  }
  /* Pad the end of OPSTR with blanks */
  while (I < MAXFK) {
    I++;
    TESTSTR[I-1] = ' ';
  }
}  /* EXTRACTQUAL */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
Local long FRAGSIZE(START, FINISH, STRANDFLAG, LINK)
long START, FINISH;
Char STRANDFLAG;
struct LOC_MODEL *LINK;
{
  long Result;

  switch (STRANDFLAG) {

  case ' ':   /* iff CIRCULAR */
    if (START < FINISH)
      Result = FINISH - START + 1;
    else
      Result = LINK->E->SEQLEN - START + FINISH + 1;
    break;

  case 'c':
  case 'C':
    if (START > FINISH)
      Result = START - FINISH + 1;
    else
      Result = LINK->E->SEQLEN - FINISH + START + 1;
    break;
  }/* case STRANDFLAG */
  return Result;
}  /* FRAGSIZE */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Read an integer from the location line. */
Local Void NUMBER(LL, LLEN, POSN, NUM, LINK)
Char *LL;
long *LLEN, *POSN, *NUM;
struct LOC_MODEL *LINK;
{
  long ORDZERO = '0';

  *NUM = 0;
  while (isdigit(LL[*POSN - 1]) && *POSN <= *LLEN) {
    *NUM = *NUM * 10 + LL[*POSN - 1] - ORDZERO;
    (*POSN)++;
  }
  fprintf(MESSAGE.f, "%ld", *NUM);   /*forces use of minimal # of chars */
}  /* NUMBER */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Read a base position from LL.  A base position composed of a single
   number is represented by setting both  LNUM and RNUM to the same
   number.  A two base bound exists when LNUM and RNUM are different. */
Local Void GETBP(LL, POSN, LLEN, BP, LINK)
Char *LL;
long *POSN, *LLEN;
BASEPOSITION *BP;
struct LOC_MODEL *LINK;
{
  /* with BP */
  if (LL[*POSN - 1] == '<' || LL[*POSN - 1] == '>') {
    BP->BOUNDFLAG = LL[*POSN - 1];
    (*POSN)++;
  } else
    BP->BOUNDFLAG = ' ';
  putc(BP->BOUNDFLAG, MESSAGE.f);
  NUMBER(LL, LLEN, POSN, &BP->LNUM, LINK);
  if (!(LL[*POSN - 1] == '.' && isdigit(LL[*POSN]) && *POSN < *LLEN)) {
    BP->RNUM = BP->LNUM;
    return;
  }
  putc('.', MESSAGE.f);
  (*POSN)++;
  NUMBER(LL, LLEN, POSN, &BP->RNUM, LINK);
}  /* GETBP */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* COPYCHAR and COPYEXPRESSION are used by ONEOF and FEANAME to copy
   an expression without evaluating it. */

/* Copy a character from LL to EXLOC */
Local Void COPYCHAR(POSN, EXLOCLEN, LL, EXLOC, LINK)
long *POSN, *EXLOCLEN;
Char *LL, *EXLOC;
struct LOC_MODEL *LINK;
{
  (*EXLOCLEN)++;
  EXLOC[*EXLOCLEN - 1] = LL[*POSN - 1];
  (*POSN)++;
}  /* COPYCHAR */

/* Copy a parenthetical expression from LL to EXLOC */
Local Void COPYEXPRESSION(POSN, EXLOCLEN, LLEN, LL, EXLOC, LINK)
long *POSN, *EXLOCLEN, *LLEN;
Char *LL, *EXLOC;
struct LOC_MODEL *LINK;
{
  boolean DONE;

  /*Copy the left parenthesis '(' */
  COPYCHAR(POSN, EXLOCLEN, LL, EXLOC, LINK);
  /* Copy the expression */
  do {
    while (LL[*POSN - 1] != ')' && LL[*POSN - 1] != '(' && *POSN <= *LLEN)
      COPYCHAR(POSN, EXLOCLEN, LL, EXLOC, LINK);
    if (LL[*POSN - 1] == '(')
      COPYEXPRESSION(POSN, EXLOCLEN, LLEN, LL, EXLOC, LINK);
    else if (LL[*POSN - 1] == ')') {
      COPYCHAR(POSN, EXLOCLEN, LL, EXLOC, LINK);
      DONE = true;
    }
  } while (!(DONE || *POSN > *LLEN));   /* COPYEXPRESSION */
}


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Create a new fragment to represent a literal sequence, and store
   the literal sequence in SEQ[>SEQLIM] */
Local Void LITERAL(LL, POSN, LLEN, RESULT, DEST, LINK)
Char *LL;
long *POSN, *LLEN;
LOCTYPE *RESULT;
NODE **DEST;
struct LOC_MODEL *LINK;
{
  long I;
  NODE *WITH;
  FRAGMENT *WITH1;
  BASEPOSITION *WITH2;

  ADDNODE(*DEST);
  *DEST = (*DEST)->NEXT;
  WITH = *DEST;
  /* with DEST^ */
  WITH->NTYPE = FRANODE;
  WITH1 = &WITH->UU.FRA;
  /* with FRA */
  WITH1->LT = baserange;
  *RESULT = WITH1->LT;
  WITH1->REPEATS = 1;
  WITH2 = &WITH1->UU.U0.START;
  /* with START */
  WITH2->LNUM = SEQLIM + 1;
  WITH2->RNUM = WITH2->LNUM;
  WITH2->BOUNDFLAG = ' ';
  WITH1->UU.U0.STRAND = ' ';
  (*POSN)++;
  I = WITH1->UU.U0.START.LNUM;
  while (LL[*POSN - 1] != '"') {
    SEQ[I-1] = NUC(LL[*POSN - 1]);
    putc(LL[*POSN - 1], MESSAGE.f);
    (*POSN)++;
    I++;
  }
  putc('\n', MESSAGE.f);
  WITH1->UU.U0.SIZE = I - WITH1->UU.U0.START.LNUM;
  WITH2 = &WITH1->UU.U0.FINISH;
  /* with FINISH */
  WITH2->LNUM = I - 1;
  WITH2->RNUM = WITH2->LNUM;
  WITH2->BOUNDFLAG = ' ';
  SEQLIM = WITH2->LNUM;
  (*POSN)++;
}  /* LITERAL */

/* Swap the values for ends; Toggle STRAND flag. */
Local Void INVERTFRAG(ANODE)
NODE **ANODE;
{
  long TEMPL, TEMPR;
  FRAGMENT *WITH;

  WITH = &(*ANODE)->UU.FRA;
  /* with ANODE^.FRA */
  if (WITH->LT == baserange) {
    TEMPR = WITH->UU.U0.START.LNUM;
    TEMPL = WITH->UU.U0.START.RNUM;
    WITH->UU.U0.START.LNUM = WITH->UU.U0.FINISH.RNUM;
    WITH->UU.U0.START.RNUM = WITH->UU.U0.FINISH.LNUM;
    WITH->UU.U0.FINISH.LNUM = TEMPL;
    WITH->UU.U0.FINISH.RNUM = TEMPR;
  }
  if (WITH->UU.U0.STRAND == ' ')
    WITH->UU.U0.STRAND = 'C';
  else
    WITH->UU.U0.STRAND = ' ';
}  /* INVERTFRAG */

/* remove a node from between two adjacent nodes */
Local Void CUTNODE(LEFT, ANODE)
NODE **LEFT, **ANODE;
{
  NODE *LBORDER, *RBORDER, *WITH;

  *ANODE = *LEFT;
  WITH = *ANODE;
  LBORDER = WITH->PREV;
  RBORDER = WITH->NEXT;
  WITH->PREV = NULL;
  WITH->NEXT = NULL;
  LBORDER->NEXT = RBORDER;
  RBORDER->PREV = LBORDER;
  *LEFT = RBORDER;
}  /* CUTNODE */

/* Place ANODE to the right of RIGHT. (Kind of like Pat Buchanan.) */
Local Void PASTENODE(ANODE, RIGHT)
NODE **ANODE, **RIGHT;
{
  NODE *RBORDER, *WITH;

  RBORDER = (*RIGHT)->NEXT;
  WITH = *ANODE;
  WITH->PREV = *RIGHT;
  WITH->NEXT = RBORDER;
  (*RIGHT)->NEXT = *ANODE;
  RBORDER->PREV = *ANODE;
}  /* PASTENODE */

/*  Given a linked list in order 1,2,3,..n,
    invert the order to n..3,2,1 */
Local Void INVERTSTRAND(CURRENTFRAG, DEST)
NODE *CURRENTFRAG, *DEST;
{
  NODE *ANODE, *LEFT = CURRENTFRAG, *RIGHT = DEST;

  /* Cut out leftmost node and reset LEFT to the leftmost remaining
     node. Paste it to the right of RIGHT. Done when LEFT = RIGHT */
  DEST = LEFT;
  while (LEFT != RIGHT) {
    CUTNODE(&LEFT, &ANODE);
    INVERTFRAG(&ANODE);
    PASTENODE(&ANODE, &RIGHT);
  }
  INVERTFRAG(&RIGHT);
}  /* INVERTSTRAND */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Complement a fragment.  This simply means to swap START and FINISH
   values and toggle the strand flag, beginning with CURRENTFRAG and
   ending with LASTFRAG. */
Local Void COMPLEMENT(LL, POSN, LLEN, RESULT, DEST, LINK)
Char *LL;
long *POSN, *LLEN;
LOCTYPE *RESULT;
NODE **DEST;
struct LOC_MODEL *LINK;
{
  NODE *CURRENTFRAG;
  BASEPOSITION DUMMY;

  putc('\n', MESSAGE.f);
  if (LL[*POSN - 1] != '(') {
    LINK->OKAY = false;
    return;
  }
  TAB(LINK);
  fprintf(MESSAGE.f, "%*c(\n", (int)LINK->INDENT, ' ');
  (*POSN)++;
  CURRENTFRAG = *DEST;
  EVALUATE(LL, POSN, LLEN, RESULT, &DUMMY, DEST, LINK);
  if (*DEST != CURRENTFRAG) {  /* one or more NODES were added */
    CURRENTFRAG = CURRENTFRAG->NEXT;
    INVERTSTRAND(CURRENTFRAG, *DEST);
  }

  if (LL[*POSN - 1] == ')') {
    fprintf(MESSAGE.f, "%*c)\n", (int)LINK->INDENT, ' ');
    (*POSN)++;
  } else {
    fprintf(MESSAGE.f, ">>> Error in database entry.\n");
    fprintf(MESSAGE.f,
	    ">>> Only one location allowed in complement expression\n");
    fprintf(MESSAGE.f, ">>> Correct syntax: complement(location)\n");
    LINK->OKAY = false;
  }
  UNTAB(LINK);
}  /* COMPLEMENT */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Join a series of locations in the order specified. */
Local Void JOIN(LL, POSN, LLEN, RESULT, DEST, LINK)
Char *LL;
long *POSN, *LLEN;
LOCTYPE *RESULT;
NODE **DEST;
struct LOC_MODEL *LINK;
{
  BASEPOSITION DUMMY;

  putc('\n', MESSAGE.f);
  if (LL[*POSN - 1] != '(') {
    LINK->OKAY = false;
    return;
  }
  TAB(LINK);
  fprintf(MESSAGE.f, "%*c(\n", (int)LINK->INDENT, ' ');
  (*POSN)++;
  do {
    /* Evaluate each local location in the location list */
    EVALUATE(LL, POSN, LLEN, RESULT, &DUMMY, DEST, LINK);
  } while (LL[*POSN - 1] != ')' && *POSN <= *LLEN && LINK->OKAY);
  if (LL[*POSN - 1] == ')') {
    fprintf(MESSAGE.f, "%*c)\n", (int)LINK->INDENT, ' ');
    (*POSN)++;
  }
  UNTAB(LINK);
}  /* JOIN */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Create a fragment for each element, but write a message that
   nothing is implied about the reasonableness of joining them.*/
Local Void ORDER(LL, POSN, LLEN, RESULT, DEST, LINK)
Char *LL;
long *POSN, *LLEN;
LOCTYPE *RESULT;
NODE **DEST;
struct LOC_MODEL *LINK;
{
  BASEPOSITION DUMMY;

  fprintf(MESSAGE.f, "\n>>> ORDER: Joining fragments.\n");
  fprintf(MESSAGE.f, ">>> User should verify validity of this action.\n");
  if (LL[*POSN - 1] != '(') {
    LINK->OKAY = false;
    return;
  }
  TAB(LINK);
  fprintf(MESSAGE.f, "%*c(\n", (int)LINK->INDENT, ' ');
  (*POSN)++;
  do {
    /* Evaluate each local location in the location list */
    EVALUATE(LL, POSN, LLEN, RESULT, &DUMMY, DEST, LINK);
  } while (LL[*POSN - 1] != ')' && *POSN <= *LLEN && LINK->OKAY);
  if (LL[*POSN - 1] == ')') {
    fprintf(MESSAGE.f, "%*c)\n", (int)LINK->INDENT, ' ');
    (*POSN)++;
  }
  UNTAB(LINK);
}  /* ORDER */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Create a fragment for each element, but write a message that
   no specific order is implied.*/
Local Void GROUP(LL, POSN, LLEN, RESULT, DEST, LINK)
Char *LL;
long *POSN, *LLEN;
LOCTYPE *RESULT;
NODE **DEST;
struct LOC_MODEL *LINK;
{
  BASEPOSITION DUMMY;

  /*!!!*/
  printf("GROUP\n");
  fprintf(MESSAGE.f, "\n>>> GROUP: Joining fragments.\n");
  fprintf(MESSAGE.f, ">>> User should verify validity of this action.\n");
  if (LL[*POSN - 1] != '(') {
    LINK->OKAY = false;
    return;
  }
  TAB(LINK);
  fprintf(MESSAGE.f, "%*c(\n", (int)LINK->INDENT, ' ');
  (*POSN)++;
  do {
    /* Evaluate each local location in the location list */
    EVALUATE(LL, POSN, LLEN, RESULT, &DUMMY, DEST, LINK);
  } while (LL[*POSN - 1] != ')' && *POSN <= *LLEN && LINK->OKAY);
  if (LL[*POSN - 1] == ')') {
    fprintf(MESSAGE.f, "%*c)\n", (int)LINK->INDENT, ' ');
    (*POSN)++;
  }
  UNTAB(LINK);
}  /* GROUP */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Evaluate  only the first expression, printing a message that
   the others are alternatives. In a future version,
   alternative choices may appear as SITES in output. */
Local Void ONEOF(LL, POSN, LLEN, RESULT, BP, DEST, LINK)
Char *LL;
long *POSN, *LLEN;
LOCTYPE *RESULT;
BASEPOSITION *BP;
NODE **DEST;
struct LOC_MODEL *LINK;
{
  LOCLINE EXLOC;
  long EXLOCLEN, I;

  putc('\n', MESSAGE.f);
  if (LL[*POSN - 1] != '(') {
    LINK->OKAY = false;
    return;
  }
  TAB(LINK);
  fprintf(MESSAGE.f, "%*c(\n", (int)LINK->INDENT, ' ');
  (*POSN)++;
  fprintf(MESSAGE.f,
    ">>> THE CURRENT IMPLEMENTATION ONLY EVALUATES THE FIRST ARGUMENT OF one-of().\n");
  fprintf(MESSAGE.f,
	  "    The result of this expression will appear in the output:\n");
  EVALUATE(LL, POSN, LLEN, RESULT, BP, DEST, LINK);
  fprintf(MESSAGE.f,
    "    The following alternate choices will NOT appear in the output:\n");
  do {
    EXLOCLEN = 0;
    COPYEXPRESSION(POSN, &EXLOCLEN, LLEN, LL, EXLOC, LINK);
    fprintf(MESSAGE.f, "%*c", (int)LINK->INDENT, ' ');
    for (I = 0; I < EXLOCLEN; I++)
      putc(EXLOC[I], MESSAGE.f);
    putc('\n', MESSAGE.f);
  } while (LL[*POSN - 1] != ')' && *POSN <= *LLEN && LINK->OKAY);
  if (LL[*POSN - 1] == ')') {
    fprintf(MESSAGE.f, "%*c)\n", (int)LINK->INDENT, ' ');
    (*POSN)++;
  }
  UNTAB(LINK);
}  /* ONEOF */

/* Local variables for REPLACE: */
struct LOC_REPLACE {
  struct LOC_MODEL *LINK;
  NODE **DEST;
} ;

/* Determine the number of the base 5' to LOC */
Local long FIVEPRIME(LOC_, FRA, STRAND, LINK)
long LOC_;
FRAGMENT FRA;
Char STRAND;
struct LOC_REPLACE *LINK;
{
  long Result;

  if (LOC_ == FRA.UU.U0.START.LNUM)
    return 0;
  switch (FRA.UU.U0.STRAND) {

  case ' ':   /* input strand */
    Result = LOC_ - 1;
    break;

  case 'c':
  case 'C':
    Result = LOC_ + 1;   /* complementary strand */
    break;
  }
  return Result;
}  /* FIVEPRIME */

/* Determine the number of the base 3' to LOC */
Local long THREEPRIME(LOC_, FRA, STRAND, LINK)
long LOC_;
FRAGMENT FRA;
Char STRAND;
struct LOC_REPLACE *LINK;
{
  long Result;

  if (LOC_ == FRA.UU.U0.FINISH.RNUM)
    return 0;
  switch (FRA.UU.U0.STRAND) {

  case ' ':   /* input strand */
    Result = LOC_ + 1;
    break;

  case 'c':
  case 'C':
    Result = LOC_ - 1;   /* complementary strand */
    break;
  }
  return Result;
}  /* THREEPRIME */

/* Splice an intron out of the object */
Local Void SPLICE(DONOR, ACCEPTOR, DONORFRAG, ACFRAG, LINK)
long *DONOR, *ACCEPTOR;
NODE **DONORFRAG, **ACFRAG;
struct LOC_REPLACE *LINK;
{
  FRAGMENT *WITH;
  BASEPOSITION *WITH1;
  NODE *WITH2;

  if (*DONOR == 0 && *ACCEPTOR == 0) {
    /* intron spans all of one
                                          or more fragments */
    RIDOF(&(*DONORFRAG)->PREV, &(*ACFRAG)->NEXT);
    return;
  }
  if (*DONOR == 0) {   /* intron at 5' end of object */
    WITH = &(*ACFRAG)->UU.FRA;
    /* ACFRAG^.FRA */

    WITH1 = &WITH->UU.U0.START;
    /* with START */
    WITH1->LNUM = *ACCEPTOR;
    WITH1->RNUM = *ACCEPTOR;
    WITH1->BOUNDFLAG = ' ';
    WITH->UU.U0.SIZE = FRAGSIZE(*ACCEPTOR, WITH->UU.U0.FINISH.RNUM,
				WITH->UU.U0.STRAND, LINK->LINK);
    return;
  }
  if (*ACCEPTOR == 0) {   /* intron at 3'end of object */
    WITH = &(*DONORFRAG)->UU.FRA;
    /* DONORFRAG^.FRA */

    WITH1 = &WITH->UU.U0.FINISH;
    /* with FINISH */
    WITH1->LNUM = *DONOR;
    WITH1->RNUM = *DONOR;
    WITH1->BOUNDFLAG = ' ';
    WITH->UU.U0.SIZE = FRAGSIZE(WITH->UU.U0.START.LNUM, *DONOR,
				WITH->UU.U0.STRAND, LINK->LINK);
    return;
  }
  /* If DONOR and ACCEPTOR are in the same fragment, then
     insert a new fragment after DONORFRAG. */
  if (*DONORFRAG == *ACFRAG) {
    ADDNODE(*DONORFRAG);
    *ACFRAG = (*DONORFRAG)->NEXT;
    *LINK->DEST = *ACFRAG;   /* make sure DEST still points to last frag. */
    WITH2 = *ACFRAG;
    /* with ACFRAG^ */
    WITH2->NTYPE = FRANODE;
    WITH = &WITH2->UU.FRA;
    /* with FRA */
    WITH->REPEATS = 1;
    WITH->LT = baserange;
    WITH->UU.U0.FINISH = (*DONORFRAG)->UU.FRA.UU.U0.FINISH;
    WITH->UU.U0.STRAND = (*DONORFRAG)->UU.FRA.UU.U0.STRAND;
  }  /* DONORFRAG=ACFRAG */

  /* Truncate DONORFRAG's 3' end at DONOR */
  WITH1 = &(*DONORFRAG)->UU.FRA.UU.U0.FINISH;
  /* with DONORFRAG^.FRA.FINISH */

  WITH1->LNUM = *DONOR;
  WITH1->RNUM = *DONOR;
  WITH1->BOUNDFLAG = ' ';
  /* Truncate ACFRAG's 5' end at ACCEPTOR */
  WITH1 = &(*ACFRAG)->UU.FRA.UU.U0.START;
  /* with ACFRAG^.FRA.START */

  WITH1->LNUM = *ACCEPTOR;
  WITH1->RNUM = *ACCEPTOR;
  WITH1->BOUNDFLAG = ' ';
  /* Recalculate the sizes of DONORFRAG and ACFRAG */
  WITH = &(*DONORFRAG)->UU.FRA;
  WITH->UU.U0.SIZE = FRAGSIZE(WITH->UU.U0.START.LNUM, WITH->UU.U0.FINISH.RNUM,
			      WITH->UU.U0.STRAND, LINK->LINK);
  if (*DONORFRAG != *ACFRAG) {
    WITH = &(*ACFRAG)->UU.FRA;
    WITH->UU.U0.SIZE = FRAGSIZE(WITH->UU.U0.START.LNUM,
	WITH->UU.U0.FINISH.RNUM, WITH->UU.U0.STRAND, LINK->LINK);
  }
}  /* SPLICE */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Replace the first location with (he contents of the second. For
   illustrative purposes, the paradigm of intron splicing is used, but
   this procedure works for all replace operations. */
Local Void REPLACE(LL, POSN, LLEN, DEST_, LINK)
Char *LL;
long *POSN, *LLEN;
NODE **DEST_;
struct LOC_MODEL *LINK;
{
  struct LOC_REPLACE V;
  LOCTYPE RESULT;
  BASEPOSITION BP;
  NODE *DONORFRAG, *ACFRAG, *DUMMY;
  long INTRON5P, INTRON3P, DONOR, ACCEPTOR;
  boolean BOTH_ENDS_FOUND = true;
  FRAGMENT *WITH;
  OBJECT *WITH1;

  /*!!!*/
  V.LINK = LINK;
  V.DEST = DEST_;
  printf("REPLACE\n");
  putc('\n', MESSAGE.f);
  if (LL[*POSN - 1] != '(') {
    LINK->OKAY = false;
    return;
  }
  TAB(LINK);
  fprintf(MESSAGE.f, "%*c(\n", (int)LINK->INDENT, ' ');
  (*POSN)++;

  /* Evaluate the target sequence. A simple position is returned as a BP.
     A base range or between position is returned as a temporary node,
     inserted into OBJ[0] by EVALUATE.  The 5' and 3' ends of the intron
     (INTRON5P and INTRON3P) are read from this node, and RIDOF removes
     it. */
  DUMMY = LINK->OBJ[0].HEAD;
  EVALUATE(LL, POSN, LLEN, &RESULT, &BP, &DUMMY, LINK);
  /* Derive DONOR and ACCEPTOR positions */
  if (RESULT == position) {
    INTRON5P = BP.LNUM;
    INTRON3P = BP.LNUM;
  } else {
    WITH = &DUMMY->UU.FRA;
    /* with DUMMY^.FRA */
    INTRON5P = WITH->UU.U0.START.LNUM;
    INTRON3P = WITH->UU.U0.FINISH.RNUM;
    RIDOF(&LINK->OBJ[0].HEAD, &LINK->OBJ[0].TAIL);
  }

  /* Find fragments in which INTRON5P and INTRON3P occur */
  WITH1 = &LINK->OBJ[LINK->CURRENT];
  /* with OBJ[CURRENT] */

  /* Set DONORFRAG and ACFRAG to point to the fragments
     containing positions INTRON5P and INTRON3P */
  DONORFRAG = WITH1->HEAD->NEXT;
  while (!WITHIN(DONORFRAG, INTRON5P) && BOTH_ENDS_FOUND) {
    if (DONORFRAG == WITH1->TAIL) {
      BOTH_ENDS_FOUND = false;
      fprintf(MESSAGE.f, ">>> replace target %12ld not found\n", INTRON5P);
    } else
      DONORFRAG = DONORFRAG->NEXT;
  }

  ACFRAG = DONORFRAG;
  while (!WITHIN(ACFRAG, INTRON3P) && BOTH_ENDS_FOUND) {
    if (ACFRAG == WITH1->TAIL) {
      BOTH_ENDS_FOUND = false;
      fprintf(MESSAGE.f, ">>> replace target %12ld not found\n", INTRON3P);
    } else
      ACFRAG = ACFRAG->NEXT;
  }
  if (BOTH_ENDS_FOUND) {
    /* Calculate DONOR and ACCEPTOR positions */
    if (RESULT == betweenposition) {
      DONOR = INTRON5P;
      ACCEPTOR = INTRON3P;
    } else {
      DONOR = FIVEPRIME(INTRON5P, DONORFRAG->UU.FRA, ' ', &V);
      ACCEPTOR = THREEPRIME(INTRON3P, ACFRAG->UU.FRA, ' ', &V);
    }

    /* Remove the target */
    SPLICE(&DONOR, &ACCEPTOR, &DONORFRAG, &ACFRAG, &V);

    /* Insert the resultant fragment(s) */
    if (LL[*POSN - 1] == '"' && LL[*POSN] == '"')   /*deletion,ie.null frag.*/
      *POSN += 2;
    else
      EVALUATE(LL, POSN, LLEN, &RESULT, &BP, &DONORFRAG, LINK);
    if (LL[*POSN - 1] == ')') {
      fprintf(MESSAGE.f, "%*c)\n", (int)LINK->INDENT, ' ');
      (*POSN)++;
    }
  }  /* BOTH_ENDS_FOUND */
  else {
    LINK->OKAY = false;
    *POSN = *LLEN + 1;
  }
  UNTAB(LINK);
}  /* REPLACE */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  THIS IS AN EXTENSION TO THE FEATURES TABLE LANGUAGE!!!
    poly(<absolute_location>|<literal>|<feature_name>,n)
    is evaluated to mean that the location is repeated n times.
    Note: in the current implementation, complex expressions involving
    functional_operators are not supported.

    Essentially, what this procedure does is to evaluate the first
    term in the expression, creating a fragment (.FRA). Then, the
    second term is evaluated to an integer, telling how many times
    the fragment is to be written by WRITEOBJ.*/
Local Void POLY(LL, POSN, LLEN, RESULT, DEST, LINK)
Char *LL;
long *POSN, *LLEN;
LOCTYPE *RESULT;
NODE **DEST;
struct LOC_MODEL *LINK;
{
  BASEPOSITION DUMMY;

  putc('\n', MESSAGE.f);
  if (LL[*POSN - 1] != '(') {
    LINK->OKAY = false;
    return;
  }
  TAB(LINK);
  fprintf(MESSAGE.f, "%*c(\n", (int)LINK->INDENT, ' ');
  (*POSN)++;

  /* Evaluate the expression */
  EVALUATE(LL, POSN, LLEN, RESULT, &DUMMY, DEST, LINK);
  TAB(LINK);
  /* Read the and set the repeat factor, REPEATS */
  fprintf(MESSAGE.f, "%*c", (int)LINK->INDENT, ' ');
  NUMBER(LL, LLEN, POSN, &(*DEST)->UU.FRA.REPEATS, LINK);
  putc('\n', MESSAGE.f);
  UNTAB(LINK);

  if (LL[*POSN - 1] == ')') {
    fprintf(MESSAGE.f, "%*c)\n", (int)LINK->INDENT, ' ');
    (*POSN)++;
  } else {
    fprintf(MESSAGE.f, ">>> TOO MANY ARGUMENTS TO POLY\n");
    LINK->OKAY = false;
  }
  UNTAB(LINK);
}  /* POLY */

/* Local variables for INTRAENTRY: */
struct LOC_INTRAENTRY {
  struct LOC_MODEL *LINK;
  FK TESTSTR;
} ;

/* Given a label, find the corresponding feature within the
   current entry. */
Local Void FINDQUALIFIER(OPSTR, QUALSTR, LPTR, LINK)
Char *OPSTR, *QUALSTR;
NODE **LPTR;
struct LOC_INTRAENTRY *LINK;
{
  boolean FOUND = false;
  long STRPOSN;
  LIST *WITH;
  LINE *WITH1;

  /*!!!*/
  if (DEBUG)
    printf("FINDQUALIFIER\n");
  WITH = &LINK->LINK->E->FEATURETABLE;
  /* with E.FEATURETABLE */
  *LPTR = WITH->HEAD->NEXT;
  while (!FOUND && *LPTR != WITH->TAIL) {
    if ((*LPTR)->UU.FEA.KEY == BLANK) {   /*ie. qualifier line */
      /* OPSTR may occur anywhere on LOCQUAL line */
      WITH1 = &(*LPTR)->UU.FEA.LOCQUAL;
      /* with LPTR^.FEA.LOCQUAL */
      /*!!!*/
      if (DEBUG)
	printf("look for /%.*s=\n", MAXFK, OPSTR);
      STRPOSN = 1;
      while (STRPOSN < WITH1->LEN && !FOUND) {
	while (WITH1->STR[STRPOSN-1] != '/' && STRPOSN < WITH1->LEN)
	  STRPOSN++;
	if (WITH1->STR[STRPOSN-1] != '/') {
	  continue;
	}  /* STR[STRPOSN]='/' */
	STRPOSN++;
	EXTRACTQUAL(WITH1->STR, &WITH1->LEN, &STRPOSN, LINK->TESTSTR,
		    LINK->LINK);
	STRPOSN++;   /* read past '=' */
	if (strncmp(LINK->TESTSTR, OPSTR, sizeof(FK))) {
	      /* look for QUALSTR */
		continue;
	}  /* TESTSTR=OPSTR */
	/*!!!*/
	if (DEBUG)
	  printf("look for %.*s\n", MAXFK, QUALSTR);
	EXTRACTQUAL(WITH1->STR, &WITH1->LEN, &STRPOSN, LINK->TESTSTR,
		    LINK->LINK);
	if (!strncmp(LINK->TESTSTR, QUALSTR, sizeof(FK)))
	  FOUND = true;
      }  /* STRPOSN<LEN and not FOUND */
    }
    if (!FOUND)
      *LPTR = (*LPTR)->NEXT;
  }  /* while not FOUND and LPTR <> TAIL */
}  /* FINDQUALIFIER */


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Resolve a reference to a labeled feature within the current entry.*/
Local Void INTRAENTRY(OPSTR, QUALSTR, RESULT, DEST, LINK)
Char *OPSTR, *QUALSTR;
LOCTYPE *RESULT;
NODE **DEST;
struct LOC_MODEL *LINK;
{
  struct LOC_INTRAENTRY V;
  NODE *LPTR;
  LOCLINE LL;
  long POSN = 1;
  /*note: this refers to position in labeled feature, NOT
                 the calling feature. */
  long LLEN;
  BASEPOSITION DUMMY;
  FEATYPE *WITH;

  V.LINK = LINK;
  /*!!!*/
  if (DEBUG)
    printf("INTRAENTRY\n");
  putc('\n', MESSAGE.f);
  /* Find the feature by its label. */
  FINDQUALIFIER(OPSTR, QUALSTR, &LPTR, &V);

  if (LPTR == LINK->E->FEATURETABLE.TAIL) {
    LINK->OKAY = false;
    fprintf(MESSAGE.f, ">>> %.*s not found\n", MAXFK, QUALSTR);
    return;
  }
  /* Back up to location part of the feature. */
  while (LPTR->UU.FEA.KEY == BLANK)
    LPTR = LPTR->PREV;

  /* Evaluate the feature */
  ONESTRING(LL, &LLEN, &LPTR, LINK);
  EVALUATE(LL, &POSN, &LLEN, RESULT, &DUMMY, DEST, LINK);

  /* write qualifier lines to MESSAGE file */
  LPTR = LPTR->NEXT;
  while (LPTR->UU.FEA.KEY == BLANK && LPTR != LINK->E->FEATURETABLE.TAIL) {
    fprintf(MESSAGE.f, "%*c", (int)LINK->INDENT, ' ');
    WITH = &LPTR->UU.FEA;
    WRITELINE(&MESSAGE, WITH->LOCQUAL, WITH->LOCQUAL.LEN);
    putc('\n', MESSAGE.f);
    LPTR = LPTR->NEXT;
  }
}  /* INTRAENTRY */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
Local Void FEANAME(LL, POSN, LLEN, OPSTR, DEST, LINK)
Char *LL;
long *POSN, *LLEN;
Char *OPSTR;
NODE **DEST;
struct LOC_MODEL *LINK;
{
  long I = 1;
  long STARTLOCATION;
  NODE *WITH;
  FRAGMENT *WITH1;
  long FORLIM;

  ADDNODE(*DEST);
  *DEST = (*DEST)->NEXT;
  WITH = *DEST;
  /* with DEST^ */
  WITH->NTYPE = FRANODE;
  WITH1 = &WITH->UU.FRA;
  /* with FRA */
  WITH1->LT = featurename;
  WITH1->REPEATS = 1;

  /* Copy OPSTR into EXLOC */
  WITH1->UU.U3.EXLOCLEN = 0;
  do {
    WITH1->UU.U3.EXLOCLEN++;
    WITH1->UU.U3.EXLOC[WITH1->UU.U3.EXLOCLEN - 1] = OPSTR[I-1];
    I++;
  } while (OPSTR[I-1] != ' ');
  STARTLOCATION = I;

  /* Copy characters from LL to EXLOC until an expression is found
     or the end of the line is reached */
  while (LL[*POSN - 1] != ')' && LL[*POSN - 1] != '(' &&
	 LL[*POSN - 1] != ',' && *POSN <= *LLEN)
    COPYCHAR(POSN, &WITH1->UU.U3.EXLOCLEN, LL, WITH1->UU.U3.EXLOC, LINK);

  /* Copy an expression from LL to EXLOC, if there is one */
  if (LL[*POSN - 1] == '(')
    COPYEXPRESSION(POSN, &WITH1->UU.U3.EXLOCLEN, LLEN, LL, WITH1->UU.U3.EXLOC,
		   LINK);
  FORLIM = WITH1->UU.U3.EXLOCLEN;
  for (I = STARTLOCATION - 1; I < FORLIM; I++)
    putc(WITH1->UU.U3.EXLOC[I], MESSAGE.f);
  putc('\n', MESSAGE.f);
}  /* FEANAME */

/* Local variables for EVALUATE: */
struct LOC_EVALUATE {
  struct LOC_MODEL *LINK;
  Char *LL;
  long *LLEN;
} ;

/* DUMMY is a temporary node added when evaluating 2nd
                  number of base_range or between_position. TAIL
                  is the next node.*/

/* This procedure checks for the special case in which a label begins
   with numeric characters. */
Local boolean NOTLABEL(POSN, LINK)
long POSN;
struct LOC_EVALUATE *LINK;
{
  if (LINK->LL[POSN-1] == '>' || LINK->LL[POSN-1] == '<')
    POSN++;
  while (isdigit(LINK->LL[POSN-1]) && POSN <= *LINK->LLEN)
    POSN++;
  if (LINK->LL[POSN-1] != '^' && LINK->LL[POSN-1] != ')' &&
      LINK->LL[POSN-1] != '(' && LINK->LL[POSN-1] != ',' &&
      LINK->LL[POSN-1] != '.' && POSN <= *LINK->LLEN)
    return false;
  else
    return true;
}  /* NOTLABEL */

/* Extract characters from the location line until a comma, left or
   right parenthesis, colon  or end of line is reached. */
Local Void EXTRACT(POSN, OPSTR, LINK)
long *POSN;
Char *OPSTR;
struct LOC_EVALUATE *LINK;
{
  long I = 0;

  while (LINK->LL[*POSN - 1] != '=' && LINK->LL[*POSN - 1] != ':' &&
	 LINK->LL[*POSN - 1] != ',' && LINK->LL[*POSN - 1] != ')' &&
	 LINK->LL[*POSN - 1] != '(' && *POSN <= *LINK->LLEN) {
    I++;
    if (I <= MAXFK)   /* use <= MAXFK chars.*/
      OPSTR[I-1] = LINK->LL[*POSN - 1];
    (*POSN)++;
  }

  /* Pad the end of OPSTR with blanks */
  while (I < MAXFK) {
    I++;
    OPSTR[I-1] = ' ';
  }
}  /* EXTRACT */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  Evaluate a location, returning a linked list after DEST.
    The procedure that calls evaluate decides what to do with the list.*/
Local Void EVALUATE(LL_, POSN, LLEN_, RESULT, BP, DEST, LINK)
Char *LL_;
long *POSN, *LLEN_;
LOCTYPE *RESULT;
BASEPOSITION *BP;
NODE **DEST;
struct LOC_MODEL *LINK;
{
  struct LOC_EVALUATE V;
  FK OPSTR, QUALSTR;
  NODE *DUMMY, *TAIL;
  FRAGMENT *WITH;

  V.LINK = LINK;
  V.LL = LL_;
  V.LLEN = LLEN_;
  TAB(LINK);
  fprintf(MESSAGE.f, "%*c", (int)LINK->INDENT, ' ');

  if (*POSN <= *V.LLEN) {
    /* Nested expression, probably base_position */
    if (V.LL[*POSN - 1] == '(') {
      fprintf(MESSAGE.f, "(\n");
      (*POSN)++;
      EVALUATE(V.LL, POSN, V.LLEN, RESULT, BP, DEST, LINK);
      if (V.LL[*POSN - 1] == ')') {
	UNTAB(LINK);
	fprintf(MESSAGE.f, "%*c)", (int)LINK->INDENT, ' ');
	(*POSN)++;
      } else
	LINK->OKAY = false;
    }

    /* base_position */
    else if ((V.LL[*POSN - 1] == '<' || V.LL[*POSN - 1] == '>' ||
	      isdigit(V.LL[*POSN - 1])) & NOTLABEL(*POSN, &V)) {
      *RESULT = position;
      GETBP(V.LL, POSN, V.LLEN, BP, LINK);
      ADDNODE(*DEST);
      *DEST = (*DEST)->NEXT;
      (*DEST)->NTYPE = FRANODE;
      WITH = &(*DEST)->UU.FRA;
      /* with DEST^.FRA */
      WITH->LT = position;
      WITH->REPEATS = 1;
      WITH->UU.U0.START = *BP;
      WITH->UU.U0.FINISH = *BP;
      WITH->UU.U0.STRAND = ' ';
      WITH->UU.U0.SIZE = 1;
    }  /* local location */

    /* literal */
    else if (V.LL[*POSN - 1] == '"') {

      /* location operator or featurename */
      LITERAL(V.LL, POSN, V.LLEN, RESULT, DEST, LINK);
    } else {
      if (V.LL[*POSN - 1] == '/') {
	putc('/', MESSAGE.f);
	(*POSN)++;
      }
      EXTRACT(POSN, OPSTR, &V);
      fprintf(MESSAGE.f, "%.*s", MAXFK, OPSTR);
      if (!strncmp(OPSTR, "complement     ", sizeof(FK)))
	COMPLEMENT(V.LL, POSN, V.LLEN, RESULT, DEST, LINK);
      else if (!strncmp(OPSTR, "join           ", sizeof(FK)))
	JOIN(V.LL, POSN, V.LLEN, RESULT, DEST, LINK);
      else if (!strncmp(OPSTR, "order          ", sizeof(FK)))
	ORDER(V.LL, POSN, V.LLEN, RESULT, DEST, LINK);
      else if (!strncmp(OPSTR, "group          ", sizeof(FK)))
	GROUP(V.LL, POSN, V.LLEN, RESULT, DEST, LINK);
      else if (!strncmp(OPSTR, "one-of         ", sizeof(FK)))
	ONEOF(V.LL, POSN, V.LLEN, RESULT, BP, DEST, LINK);
      else if (!strncmp(OPSTR, "replace        ", sizeof(FK)))
	REPLACE(V.LL, POSN, V.LLEN, DEST, LINK);
      else if (!strncmp(OPSTR, "poly           ", sizeof(FK)))
	POLY(V.LL, POSN, V.LLEN, RESULT, DEST, LINK);
      else if (V.LL[*POSN - 1] == ':')
	FEANAME(V.LL, POSN, V.LLEN, OPSTR, DEST, LINK);
      else {
	if (V.LL[*POSN - 1] == '=') {
	  putc('=', MESSAGE.f);
	  (*POSN)++;
	  EXTRACT(POSN, QUALSTR, &V);
	  fprintf(MESSAGE.f, "%.*s", MAXFK, QUALSTR);
	} else {
	  memcpy(QUALSTR, OPSTR, sizeof(FK));
	  memcpy(OPSTR, "label          ", sizeof(FK));
	}
	INTRAENTRY(OPSTR, QUALSTR, RESULT, DEST, LINK);
	/* INTRAENTRY */
      }
    }
  }

  /* location operator or featurename */
  /* Evaluate the second part of a base range or between position,
     if any. */
  if (*POSN < *V.LLEN) {
    if (V.LL[*POSN - 1] == '^' || V.LL[*POSN - 1] == '.') {
      (*POSN)++;
      WITH = &(*DEST)->UU.FRA;
      /* with DEST^.FRA */

      if (V.LL[*POSN - 1] == '.') {
	WITH->LT = baserange;
	(*POSN)++;
      } else
	WITH->LT = betweenposition;

      /* Get the next number */
      EVALUATE(V.LL, POSN, V.LLEN, RESULT, BP, DEST, LINK);
      /* Assign pointers so that DUMMY points to the temporary
         node and TAIL to the next one after that. This makes
         it possible to call RIDOF to remove the DUMMY node.*/
      DUMMY = *DEST;
      TAIL = DUMMY->NEXT;
      *DEST = DUMMY->PREV;

      /* Update DEST with information from the second part of
         base_range or between_position. */
      WITH = &(*DEST)->UU.FRA;
      /* with DEST^.FRA */
      WITH->UU.U0.FINISH = *BP;
      WITH->UU.U0.SIZE = FRAGSIZE(WITH->UU.U0.START.LNUM,
	  WITH->UU.U0.FINISH.RNUM, WITH->UU.U0.STRAND, LINK);
      *RESULT = WITH->LT;
      RIDOF(DEST, &TAIL);   /* remove dummy fragment */
      putc('\n', MESSAGE.f);
      UNTAB(LINK);
    }  /* base range or between position */
    else {
      UNTAB(LINK);
      putc('\n', MESSAGE.f);
    }
  } else {
    UNTAB(LINK);
    putc('\n', MESSAGE.f);
  }

  /* Read past delimiters to next expression or end of LL array */
  if (*POSN <= *V.LLEN) {
    if (V.LL[*POSN - 1] == ',')
      (*POSN)++;
  }
}  /* EVALUATE */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Add a fragment to the begining of a CDS, sig_peptide, or mat_peptide,
   containing the result of the literal expression "n" or "nn", which will
   result in adding n's to the sequence, to put it in the correct reading
   frame. */
Local Void PADCODON(QUAL, NUMN, LINK)
LINE *QUAL;
long *NUMN;
struct LOC_MODEL *LINK;
{
  long POSN = 2;
  long LLEN;
  Char READINGFRAME;
  FK TESTSTR;
  LOCLINE LL;
  NODE *HEADFRAG;
  LOCTYPE RESULT;   /* dummy */
  BASEPOSITION BP;   /* dummy */

  if (QUAL->STR[0] != '/') {
    return;
  }  /* STR[1]='/' */
  /* Extract the qualifier TESTSTR from QUAL, and if it is 'codon_
     start', set READINGFRAME to the position indicated (1,2 or 3). If
     READINGFRAME is 1, do nothing. If READINGFRAME = 2 or 3, create a
     LOCLINE containing a literal expression, to be evaluated by
     EVALUATE, and inserted after the HEAD of OBJ[CURRENT] */
  EXTRACTQUAL(QUAL->STR, &QUAL->LEN, &POSN, TESTSTR, LINK);
  if (strncmp(TESTSTR, "codon_start    ", sizeof(FK))) {
    return;
  }  /* TESTSTR='codon_start' */
  READINGFRAME = QUAL->STR[POSN];
  if (READINGFRAME != '3' && READINGFRAME != '2') {
    *NUMN = 0;
    return;
  }  /* READINGFRAME in ['2','3'] */
  fprintf(MESSAGE.f, "n's added to 5' end to preserve reading frame: ");
  /* Create a LOCLINE to be evaluated: "n" or "nn"  */
  LL[0] = '"';
  LL[1] = 'n';
  if (READINGFRAME == '3') {
    LLEN = 3;
    LL[2] = '"';
    *NUMN = 1;
  } else {
    LLEN = 4;
    LL[2] = 'n';
    LL[3] = '"';
    *NUMN = 2;
  }
  /* Evaluate the expression, and add it to the head of the
     object */
  HEADFRAG = LINK->OBJ[LINK->CURRENT].HEAD;
  POSN = 1;
  EVALUATE(LL, &POSN, &LLEN, &RESULT, &BP, &HEADFRAG, LINK);
}  /* PADCODON */

/* Local variables for MAKETITLE: */
struct LOC_MAKETITLE {
  struct LOC_MODEL *LINK;
  long I, POSN;
} ;

/* Read the feature label, if any, from the qualifier line. */
Local Void FINDLABEL(QUALIFIER, FEALAB, NOLABEL, LINK)
LINE QUALIFIER;
WORD *FEALAB;
boolean *NOLABEL;
struct LOC_MAKETITLE *LINK;
{
  FK TESTSTR;

  LINK->POSN = 2;   /* read past '/' */
  EXTRACTQUAL(QUALIFIER.STR, &QUALIFIER.LEN, &LINK->POSN, TESTSTR, LINK->LINK);
  if (strncmp(TESTSTR, "label          ", sizeof(FK)))
    return;
  *NOLABEL = false;
  LINK->POSN++;   /* read past '=' */
  /* Copy the label to FEALAB */
  /* with FEALAB */
  FEALAB->LEN = 0;
  for (LINK->I = LINK->POSN; LINK->I <= QUALIFIER.LEN; LINK->I++) {
    FEALAB->LEN++;
    FEALAB->STR[FEALAB->LEN - 1] = QUALIFIER.STR[LINK->I-1];
  }
  for (LINK->I = FEALAB->LEN + 1; LINK->I <= MAXWORD; LINK->I++)
    FEALAB->STR[LINK->I-1] = ' ';
}  /* FINDLABEL */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/*  <title>::=  <locus>:<label>|
                <locus>:<feature key>  */
Local Void MAKETITLE(LPTR, FEAKEY, OBJ, LINK)
NODE *LPTR;
FEATUREKEY FEAKEY;
OBJECT *OBJ;
struct LOC_MODEL *LINK;
{
  struct LOC_MAKETITLE V;
  long J;
  Char CH;
  FEATYPE *WITH1;
  WORD *WITH2;
  long FORLIM;

  V.LINK = LINK;
  /* Title begins with the locus name */
  /* with OBJ */
  OBJ->TITLE = LINK->E->NAME;

  /* Find the feature label, if any. */
  OBJ->NOLABEL = true;
  while (LPTR->UU.FEA.KEY == BLANK && LPTR != LINK->E->FEATURETABLE.TAIL &&
	 OBJ->NOLABEL) {
    WITH1 = &LPTR->UU.FEA;
    FINDLABEL(WITH1->LOCQUAL, &OBJ->FEALAB, &OBJ->NOLABEL, &V);
    LPTR = LPTR->NEXT;
  }

  /* If there isn't a label, add featurekey to the title, otherwise,
     add the label. */
  WITH2 = &OBJ->TITLE;
  /* with TITLE */
  WITH2->LEN++;
  WITH2->STR[WITH2->LEN - 1] = ':';
  if (OBJ->NOLABEL) {
    CH = KEYSTR[(long)FEAKEY][0];
    J = 1;
    while (CH != ' ' && J < MAXFK) {
      WITH2->LEN++;
      WITH2->STR[WITH2->LEN - 1] = CH;
      J++;
      if (J <= MAXFK)
	CH = KEYSTR[(long)FEAKEY][J-1];
    }
    return;
  }  /* NOLABEL */
  FORLIM = OBJ->FEALAB.LEN;
  for (J = 0; J < FORLIM; J++) {
    WITH2->LEN++;
    WITH2->STR[WITH2->LEN - 1] = OBJ->FEALAB.STR[J];
  }
  for (J = WITH2->LEN; J < MAXWORD; J++)
    WITH2->STR[J] = ' ';
}  /* MAKETITLE */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
Local Void BUILD(LL, LLEN, OBJ, LINK)
Char *LL;
long *LLEN;
OBJECT *OBJ;
struct LOC_MODEL *LINK;
{
  long POSN = 1;
  LOCTYPE RESULT;
  NODE *LASTFRAG;
  BASEPOSITION BP;

  if (RESOLVEXP) {
    WRITEWORD(&MESSAGE, ACCESSION, ACCESSION.LEN);
    fprintf(MESSAGE.f, ":%ld\n", LINK->CURRENT);
  } else {
    WRITEWORD(&MESSAGE, OBJ->TITLE, OBJ->TITLE.LEN);
    if (OBJ->NOLABEL)
      fprintf(MESSAGE.f, "%ld", LINK->CURRENT);
    putc('\n', MESSAGE.f);
  }
  LINK->OKAY = true;
  LINK->INDENT = 0;
  LASTFRAG = OBJ->HEAD;
  EVALUATE(LL, &POSN, LLEN, &RESULT, &BP, &LASTFRAG, LINK);

  /* Although not permitted by the Features Table Definition,
     GenBank has been known to allow locations that evaluate to a
     between position. Since this is biologically meaningless, it
     must be detected and eliminated. */
  if (RESULT == betweenposition) {
    fprintf(MESSAGE.f, ">>> FEATURE EVALUATES TO between_position.\n");
    fprintf(MESSAGE.f, ">>> NO SEQUENCE RESULTS FROM THIS EXPRESSION.\n");
    LINK->OKAY = false;
  }
  putc('\n', MESSAGE.f);
  if (LINK->OKAY) {
    return;
  }  /* if not OKAY */
  RIDOF(&OBJ->HEAD, &OBJ->TAIL);
  OBJ->TITLE.LEN = 0;
  OBJ->FEALAB.LEN = 0;
  OBJ->NOLABEL = true;
  NUMOBJ--;
}  /* BUILD */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* this is a non-functional skeleton procedure, at present */
Local Void FINDSITE(LL, LLEN, OBJ, LINK)
Char *LL;
long *LLEN;
OBJECT *OBJ;
struct LOC_MODEL *LINK;
{
  long POSN = 1;
  LOCTYPE RESULT;
  NODE *LASTFRAG;
  BASEPOSITION BP;

  LINK->OKAY = true;
  LINK->INDENT = 0;
  LASTFRAG = OBJ->HEAD;
  EVALUATE(LL, &POSN, LLEN, &RESULT, &BP, &LASTFRAG, LINK);
  putc('\n', MESSAGE.f);
  if (!LINK->OKAY) {
    RIDOF(&OBJ->HEAD, &OBJ->TAIL);
    NUMOBJ--;
  }  /* if not OKAY */
}  /* FINDSITE */


/****************************************************************/
/* Make a model of the sequence, based on criteria in datafile. */
/****************************************************************/

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* Forward declaration for recursive procedure EVALUATE */
Static Void MODEL(E_, TRANSCRIPTS, SITESET, LL, LLEN, OBJ_)
ENTRY *E_;
long *TRANSCRIPTS, *SITESET;
Char *LL;
long *LLEN;
OBJECT *OBJ_;
{  /* MODEL -------------------------------------------------------*/
  struct LOC_MODEL V;
  boolean DONE, TRANSLATED;   /* =true if feature codes for protein */
  NODE *TEMP;
  long NUMN;   /* # of n's to pad CDS with if codon_start=2 or 3 */
  long I;   /* index of the current object */
  FEATUREKEY FEAKEY;   /* Feature key of current feature */
  LIST *WITH;
  OBJECT *WITH1;
  long FORLIM;
  NODE *WITH2;


  V.E = E_;
  V.OBJ = OBJ_;
  /* If RESOLVEXP is true, then simply evaluate the location found in LL.
     Otherwise, make a pass through the FEATURES table, creating one or
     more objects corresponding to TESTSET.  Finally, make a second pass
     through the data, searching for each site in SITESET and assigning
     it to a particular object. */
  WITH = &V.E->FEATURETABLE;

  /* with E.FEATURETABLE */
  /* 1st Pass - Create objects */
  if (RESOLVEXP) {
    NUMOBJ = 1;
    V.CURRENT = NUMOBJ;
    /* any procedure can refer to CURRENT to find
                        out the current OBJECT.  This means that
                        OBJ doesn't have to be passed to each proc.*/
    BUILD(LL, LLEN, &V.OBJ[V.CURRENT], &V);
    fprintf(MESSAGE.f, "//----------------------------------------------\n");
  }  /* RESOLVEXP */
  else {
    NUMOBJ = 0;
    TEMP = WITH->HEAD->NEXT;
    while (TEMP != WITH->TAIL) {
      FEAKEY = TEMP->UU.FEA.KEY;
      if (!P_inset(FEAKEY, TRANSCRIPTS)) {  /* an object */
	TEMP = TEMP->NEXT;
	continue;
      }  /*  an object */
      if (FEAKEY == (int)SIGPEPTIDE || FEAKEY == (int)MATPEPTIDE ||
	  FEAKEY == (int)CDS)
	TRANSLATED = true;
      else
	TRANSLATED = false;
      NUMN = 0;
      NUMOBJ++;
      V.CURRENT = NUMOBJ;
      /* any procedure can refer to CURRENT to find
                          out the current OBJECT.  This means that
                        OBJ doesn't have to be passed to each proc.*/

      /* Re-code the location to a single line for easy parsing */
      ONESTRING(LL, LLEN, &TEMP, &V);

      /* Each object carries a unique title, to identify the object
         in all output files.*/
      MAKETITLE(TEMP, FEAKEY, &V.OBJ[V.CURRENT], &V);

      /* Build the object */
      BUILD(LL, LLEN, &V.OBJ[V.CURRENT], &V);

      /* Write qualifier lines to message file */
      DONE = false;
      while (!DONE) {
	if (TEMP == WITH->TAIL) {
	  DONE = true;
	  break;
	}
	if (TEMP->UU.FEA.KEY != (int)OTHER && TEMP->UU.FEA.KEY != (int)BLANK) {
	  DONE = true;
	  break;
	}
	/* Write the qualifier */
	WRITELINE(&MESSAGE, TEMP->UU.FEA.LOCQUAL, TEMP->UU.FEA.LOCQUAL.LEN);
	putc('\n', MESSAGE.f);

	/* If the qualifier is /codon_start, then add a fragment to
	   the head of the object, containing n's to put the object
	   in the correct reading frame. */
	if (TRANSLATED && PAD5PRIME)
	  PADCODON(&TEMP->UU.FEA.LOCQUAL, &NUMN, &V);
	TEMP = TEMP->NEXT;
      }
      fprintf(MESSAGE.f, "//----------------------------------------------\n");
      /* output to EXPFILE */
      if (!V.OKAY)
	continue;
      WITH1 = &V.OBJ[V.CURRENT];
      /*  ><title> */
      /* if OKAY */
      putc('>', EXPFILE.f);
      WRITEWORD(&EXPFILE, WITH1->TITLE, WITH1->TITLE.LEN);
      if (WITH1->NOLABEL)
	fprintf(EXPFILE.f, "%ld", V.CURRENT);
      putc('\n', EXPFILE.f);

      /* if codon_start<>1, write n's before expression */
      if ((unsigned long)NUMN < 32 && ((1L << NUMN) & 0x6) != 0) {
	for (I = 1; I <= NUMN; I++)
	  putc('n', EXPFILE.f);
	putc('\n', EXPFILE.f);
      }

      /* Write primary accession number etc.: @<accession>: */
      putc('@', EXPFILE.f);
      WITH2 = V.E->ACCESSION.HEAD->NEXT;
      WRITEWORD(&EXPFILE, WITH2->UU.WOR, WITH2->UU.WOR.LEN);
      putc(':', EXPFILE.f);
      /* <expression> */
      if (WITH1->NOLABEL) {
	FORLIM = *LLEN;
	for (I = 0; I < FORLIM; I++)
	  putc(LL[I], EXPFILE.f);
      } else
	WRITEWORD(&EXPFILE, WITH1->FEALAB, WITH1->FEALAB.LEN);
      putc('\n', EXPFILE.f);
    }  /* TEMP <> TAIL */
  }
  /* 1st Pass */

  /* 2nd Pass -  Assign sites to objects */
  if (NUMOBJ <= 0 || *SITESET == 0) {
    return;
  }  /* NUMOBJ > 0 */
  fprintf(MESSAGE.f, "SITES:\n");
  TEMP = WITH->HEAD->NEXT;
  while (TEMP != WITH->TAIL) {
    if (P_inset(TEMP->UU.FEA.KEY, SITESET)) {  /* a site */
      FINDSITE(LL, LLEN, V.OBJ, &V);
    }  /*  a site */
    else
      TEMP = TEMP->NEXT;
  }  /* TEMP <> TAIL */
}  /* MODEL */

#undef TABSIZE


Local Void WRITEFRAG(F, START, FINISH, SENSE, REPEATS, NOOKIE, WIDTH)
_TEXT *F;
long START, FINISH, SENSE, REPEATS;
NUCLEOTIDE *NOOKIE;
long *WIDTH;
{
  long I, J;

  for (J = 1; J <= REPEATS; J++) {
    I = START;
    do {
      putc(NUCHAR[(long)NOOKIE[(long)SEQ[I-1]]], F->f);
      (*WIDTH)++;
      if (*WIDTH == 50) {
	putc('\n', F->f);
	*WIDTH = 0;
      }
      I += SENSE;
    } while (I != FINISH + SENSE);   /* for J */
  }
}  /* WRITEFRAG */


/****************************************/
/* Write the object to the output file. */
/****************************************/
Static Void WRITEOBJECT(F, OBJ_)
_TEXT *F;
OBJECT *OBJ_;
{
  OBJECTS OBJ;
  long FIRST, LAST, I, J, K_, WIDTH;
  NODE *TEMP;
  long FORLIM;
  OBJECT *WITH;
  FRAGMENT *WITH1;
  long FORLIM1, FORLIM2;


  memcpy(OBJ, OBJ_, sizeof(OBJECTS));
  FORLIM = NUMOBJ;
  for (I = 1; I <= FORLIM; I++) {
    WITH = &OBJ[I];
    /* Write the object to F */
    if (!RESOLVEXP) {
      putc('>', F->f);
      WRITEWORD(F, WITH->TITLE, WITH->TITLE.LEN);
      if (WITH->NOLABEL)
	fprintf(F->f, "%ld", I);
      putc('\n', F->f);
    }
    TEMP = WITH->HEAD->NEXT;
    do {
      WITH1 = &TEMP->UU.FRA;

      if (((1L << ((long)WITH1->LT)) &
	   ((1L << ((long)baserange)) | (1L << ((long)position)))) != 0) {
	WIDTH = 0;
	FIRST = WITH1->UU.U0.START.LNUM;
	LAST = WITH1->UU.U0.FINISH.RNUM;
	switch (WITH1->UU.U0.STRAND) {

	case ' ':
	  if (FIRST <= LAST)
	    WRITEFRAG(F, FIRST, LAST, 1L, WITH1->REPEATS, INUC, &WIDTH);
	  else {
	    WRITEFRAG(F, FIRST, E.SEQLEN, 1L, WITH1->REPEATS, INUC, &WIDTH);
	    WRITEFRAG(F, 1L, LAST, 1L, WITH1->REPEATS, INUC, &WIDTH);
	  }
	  break;

	case 'C':
	case 'c':
	  if (FIRST >= LAST)
	    WRITEFRAG(F, FIRST, LAST, -1L, WITH1->REPEATS, COMP, &WIDTH);
	  else {
	    WRITEFRAG(F, FIRST, 1L, -1L, WITH1->REPEATS, COMP, &WIDTH);
	    WRITEFRAG(F, E.SEQLEN, LAST, -1L, WITH1->REPEATS, COMP, &WIDTH);
	  }
	  break;
	}/* case STRAND */
      }  /* baserange */
      else {
	FORLIM1 = WITH1->REPEATS;
	for (K_ = 1; K_ <= FORLIM1; K_++) {
	  putc('@', F->f);
	  FORLIM2 = WITH1->UU.U3.EXLOCLEN;
	  for (J = 0; J < FORLIM2; J++)
	    putc(WITH1->UU.U3.EXLOC[J], F->f);
	}  /* for K */
      }
      /* featurename */
      /* featurename */
      putc('\n', F->f);
      TEMP = TEMP->NEXT;
    } while (TEMP != WITH->TAIL);
  }
}  /* WRITEOBJECT */


typedef Char KEYWORD[12];


/* Local variables for WRITEGB: */
struct LOC_WRITEGB {
  _TEXT *F;
} ;

Local Void WRITELIST(HEAD, TAIL, KW, LINK)
NODE *HEAD, *TAIL;
Char *KW;
struct LOC_WRITEGB *LINK;
{
  NODE *TEMP;
  long LINEWIDTH;
  FEATYPE *WITH;

  if (HEAD->NEXT == TAIL) {
    return;
  }  /* HEAD */
  fprintf(LINK->F->f, "%.12s", KW);
  TEMP = HEAD->NEXT;
  switch (TEMP->NTYPE) {

  case WNODE:
    LINEWIDTH = 12;
    do {
      LINEWIDTH += TEMP->UU.WOR.LEN + 1;
      if (LINEWIDTH > 80) {
	fprintf(LINK->F->f, "\n%12c", ' ');
	LINEWIDTH = TEMP->UU.WOR.LEN + 13;
      }
      WRITEWORD(LINK->F, TEMP->UU.WOR, TEMP->UU.WOR.LEN);
      putc(' ', LINK->F->f);
      TEMP = TEMP->NEXT;
    } while (TEMP != TAIL);
    putc('\n', LINK->F->f);
    break;

  case LNODE:
    WRITELINE(LINK->F, TEMP->UU.LIN, TEMP->UU.LIN.LEN);
    putc('\n', LINK->F->f);
    TEMP = TEMP->NEXT;
    while (TEMP != TAIL) {
      fprintf(LINK->F->f, "%12c", ' ');
      WRITELINE(LINK->F, TEMP->UU.LIN, TEMP->UU.LIN.LEN);
      putc('\n', LINK->F->f);
      TEMP = TEMP->NEXT;
    }
    break;

  case FRANODE:
    /* blank case */
    break;

  case FEANODE:
    if (!strncmp(KW, "FEATURES    ", sizeof(KEYWORD)))
      fprintf(LINK->F->f, "        Location/Qualifiers\n");
    else
      putc('\n', LINK->F->f);
    do {
      WITH = &TEMP->UU.FEA;
      /* with TEMP.FEA */
      if (WITH->KEY == OTHER)
	fprintf(LINK->F->f, "    %.*s ", MAXFK, WITH->KEYSTRING);
      else
	fprintf(LINK->F->f, "    %.*s ", MAXFK, KEYSTR[(long)WITH->KEY]);
      WRITELINE(LINK->F, WITH->LOCQUAL, WITH->LOCQUAL.LEN);
      putc('\n', LINK->F->f);
      TEMP = TEMP->NEXT;
    } while (TEMP != TAIL);   /* FEANODE */
    break;
  }
}  /* WRITELIST */


/****************************************/
/* Write the GenBank entry to the output*/
/****************************************/
Static Void WRITEGB(F_, E)
_TEXT *F_;
ENTRY E;
{
  struct LOC_WRITEGB V;

  V.F = F_;
  /* Write the annotation data */
  /* with E */
  fprintf(V.F->f, "LOCUS       ");
  WRITEWORD(V.F, E.NAME, 18L);
  fprintf(V.F->f, "%11ld bp\n", E.SEQLEN);
  WRITELIST(E.DEFINITION.HEAD, E.DEFINITION.TAIL, "DEFINITION  ", &V);
  WRITELIST(E.ACCESSION.HEAD, E.ACCESSION.TAIL, "ACCESSION   ", &V);
  /* NID field to be removed from GenBank by Dec. 1999. However,
     this statement will write NID fields from old files, if they
     exist. */
  if (E.NID.LEN > 0) {
    fprintf(V.F->f, "NID         ");
    WRITEWORD(V.F, E.NID, E.NID.LEN);
    putc('\n', V.F->f);
  }
  if (E.SEGTOTAL > 0)
    fprintf(V.F->f, "SEGMENT     %ld of %ld\n", E.SEGNUMBER, E.SEGTOTAL);
  WRITELIST(E.FEATURETABLE.HEAD, E.FEATURETABLE.TAIL, "FEATURES    ", &V);
  fprintf(V.F->f, "BASECOUNT   %7ld a%7ld c%7ld g%7ld t\n",
	  E.ACOMP, E.CCOMP, E.GCOMP, E.TCOMP);

  fprintf(V.F->f, "//\n");
}  /* WRITEGB */


/*****************************************************/
/* Reinitialize the entry.                           */
/*****************************************************/
Static Void REINIT(E)
ENTRY *E;
{
  long I, FORLIM;
  OBJECT *WITH1;

  E->NAME.LEN = 0;
  E->SEQLEN = 0;
  E->ENTRYSTATUS.LEN = 0;
  RIDOF(&E->DEFINITION.HEAD, &E->DEFINITION.TAIL);
  RIDOF(&E->ACCESSION.HEAD, &E->ACCESSION.TAIL);
  E->NID.LEN = 0;
  E->SEGTOTAL = 0;
  RIDOF(&E->FEATURETABLE.HEAD, &E->FEATURETABLE.TAIL);
  FORLIM = NUMOBJ;
  for (I = 1; I <= FORLIM; I++) {
    WITH1 = &OBJ[I];
    RIDOF(&WITH1->HEAD, &WITH1->TAIL);
    WITH1->TITLE.LEN = 0;
    WITH1->FEALAB.LEN = 0;
    WITH1->NOLABEL = true;
  }
  NUMOBJ = 0;
}  /* REINIT */


/* -----------------------------------------------------------  */
/* ----------------- MAIN  PROCEDURE  -----------------------  */
main(argc, argv)
int argc;
Char *argv[];
{
  PASCAL_MAIN(argc, argv);
  EXPFILE.f = NULL;
  strcpy(EXPFILE.name, "EXPFILE");
  OUTFILE.f = NULL;
  strcpy(OUTFILE.name, "OUTFILE");
  MESSAGE.f = NULL;
  strcpy(MESSAGE.name, "MESSAGE");
  INDFILE.f = NULL;
  strcpy(INDFILE.name, "INDFILE");
  SEQFILE.f = NULL;
  strcpy(SEQFILE.name, "SEQFILE");
  ANOFILE.f = NULL;
  strcpy(ANOFILE.name, "ANOFILE");
  NAMEFILE.f = NULL;
  strcpy(NAMEFILE.name, "NAMEFILE");
  INFILE.f = NULL;
  /* - - - - Search for specified features in each entry - - - - */


  /*!!!    if not FILES then CLOSE(OUTFILE) */
  strcpy(INFILE.name, "INFILE");
  /* Initialize parameters. */
  INITIALIZE();

  /* Read options from command line and set FILES,RESOLVEXP & ACNO */
  ARGNUM = STARTARGNUM;
  READOPTIONS(&ARGNUM, &FILES, &RESOLVEXP, &ACNO);

  /* Open files.*/
  OPENFILES(&ARGNUM, &INFILE, &NAMEFILE, &ANOFILE, &SEQFILE, &INDFILE,
	    &MESSAGE, &OUTFILE, &EXPFILE);

  /* Read input file containing instructions for processing each
        entry.*/
  READINSTRUCTIONS(&INFILE, TRANSCRIPTS, SITESET, OBJECTTYPE);

  fprintf(MESSAGE.f, "%s\n", VERSION);
  fprintf(MESSAGE.f,
	  "Please cite: Fristensky B. (1993) Feature expressions:\n");
  fprintf(MESSAGE.f, "creating and manipulating sequence datasets.\n");
  fprintf(MESSAGE.f, "Nucl. Acids Res. 21:5997-6003\n\n");
  CA = 1;
  CS = 1;

  if (RESOLVEXP) {   /*ie. -r */
    /*- - - - Resolve objects from external references  - - - - - -*/
    while (!BUFEOF(NAMEFILE.f)) {
      /* Copy lines from NAMEFILE to OUTFILE. These may be objects
         resolved in previous runs of GETOB */
      while ((P_peek(NAMEFILE.f) != '@') & (!BUFEOF(NAMEFILE.f))) {
	while (!P_eoln(NAMEFILE.f)) {
	  CH = getc(NAMEFILE.f);
	  if (CH == '\n')
	    CH = ' ';
	  putc(CH, OUTFILE.f);
	}
	if (!BUFEOF(NAMEFILE.f)) {
	  fscanf(NAMEFILE.f, "%*[^\n]");
	  getc(NAMEFILE.f);
	}
	putc('\n', OUTFILE.f);
      }  /* NAMEFILE^<>'@' */

      /* Process external reference. These are references left un-
          resolved in a previous run of GETOB. */
      if (P_peek(NAMEFILE.f) != '@') {
	continue;
      }  /* NAMEFILE^='@' */
      READLINE(&NAMEFILE, &INDIRREF);
      PARSEREF(&INDIRREF, &ACCESSION, LL, &LLEN, &OKAY);
      FOUND = false;
      if (OKAY)
	FINDDATA(&INDFILE, ACCESSION, &LOCATIONS, &FOUND);
      if (FOUND && OKAY) {
	RGB(&ANOFILE, &SEQFILE, LOCATIONS, &E);
	MODEL(&E, TRANSCRIPTS, SITESET, LL, &LLEN, OBJ);
	WRITEOBJECT(&OUTFILE, OBJ);
	REINIT(&E);
	continue;
      }  /* FOUND */
      fprintf(MESSAGE.f, ">>> Can't resolve indirect reference: ");
      WRITELINE(&MESSAGE, INDIRREF, INDIRREF.LEN);
      putc('\n', MESSAGE.f);
      putc(';', OUTFILE.f);
      WRITELINE(&OUTFILE, INDIRREF, INDIRREF.LEN);
      putc('\n', OUTFILE.f);
    }  /* eof(NAMEFILE) */

  } else {
    /* print message; comment out expression in OUTFILE*/
    /* For each locus in NAMEFILE, find the locations of the annotation
       and sequence parts in ANOFILE and SEQFILE, and then get the data
       and write to OUTFILE format specified by INFILE. */
    while (!BUFEOF(NAMEFILE.f)) {
      READNAME(&NAMEFILE, &SEQNAME);
      FOUND = false;
      if (SEQNAME.LEN > 0)
	FINDDATA(&INDFILE, SEQNAME, &LOCATIONS, &FOUND);
      if (!FOUND) {
	continue;
      }  /* if FOUND */
      if (FILES) {
	MAKEFN(SEQNAME, &FILENAME);
	strcpy(OUTFILE.name, P_trimname(FILENAME.STR, sizeof(CHARARRAY)));
	if (OUTFILE.f != NULL)
	  OUTFILE.f = freopen(OUTFILE.name, "w", OUTFILE.f);
	else
	  OUTFILE.f = fopen(OUTFILE.name, "w");
	if (OUTFILE.f == NULL)
	  _EscIO2(FileNotFound, OUTFILE.name);
	SETUPBUF(OUTFILE.f, Char);
      }  /* FILES */
      RGB(&ANOFILE, &SEQFILE, LOCATIONS, &E);
      if (!strncmp(OBJECTTYPE, "GENBANK        ", sizeof(FK)))
	WRITEGB(&OUTFILE, E);
      else {
	MODEL(&E, TRANSCRIPTS, SITESET, LL, &LLEN, OBJ);
	WRITEOBJECT(&OUTFILE, OBJ);
      }
      REINIT(&E);
      /*!!!        if FILES then CLOSE(OUTFILE) */
    }  /* not eof(NAMEFILE) */
  }
  if (INFILE.f != NULL)
    fclose(INFILE.f);
  if (NAMEFILE.f != NULL)
    fclose(NAMEFILE.f);
  if (ANOFILE.f != NULL)
    fclose(ANOFILE.f);
  if (SEQFILE.f != NULL)
    fclose(SEQFILE.f);
  if (INDFILE.f != NULL)
    fclose(INDFILE.f);
  if (MESSAGE.f != NULL)
    fclose(MESSAGE.f);
  if (OUTFILE.f != NULL)
    fclose(OUTFILE.f);
  if (EXPFILE.f != NULL)
    fclose(EXPFILE.f);
  exit(EXIT_SUCCESS);
}  /* GETOB */



/* End. */
