/* Output from p2c 1.21alpha-07.Dec.93, the Pascal-to-C translator */
/* From input file "intrest.p" */


/**********************************************************/
/*                                                        */
/*     INTREST Version  8/14/2001, Standard Pascal        */
/*             Brian Fristensky                           */
/*             Dept. of Plant Science                     */
/*             University of Manitoba                     */
/*             Winnipeg, MB R3T 2N2  CANADA               */
/*                                                        */
/* Copyright (c) 1986,1987,1988,1990 by Brian Fristensky  */
/* !!! in comment indicates feature which may need change */
/******************************************************** */

#include <p2c.h>


/*,INFILE,OUTFILE*/
/*!!! Some Pascals may require file parameters in program heading */

#define MAXSEQ          750000L

/* BEGIN MODULE INTLIMITS */
/*!!!  MAXINT =  2147483647; */
/*!!!  MININT = -2147483647; */
/* END MODULE INTLIMITS         VERSION= 'SUNMODS     Version  6/26/01'; */

#define MAXPAT          23
#define MAXFRAGS        6000   /*max. # of fragmenst in a single digest */
#define MAXWORD         23
#define MAXLINE         150

#define VERSION         "INTREST   Version  8/14/2001"


typedef enum {
  CANT, A, C, R, D, V, M, K, B, H, Y, G, T, W, S, N, WONT
} NUCLEOTIDE;
typedef long SS;

typedef NUCLEOTIDE SEQUENCE[MAXSEQ + MAXPAT + 1];

/* BEGIN MODULE TYPE.WORD */
/*   <word>::= <non-blank char>[<non-blank char>] */

typedef struct WORD {
  long LEN;
  Char STR[MAXWORD];
} WORD;

/* END MODULE TYPE.WORD         VERSION= 'SUNMODS     Version  6/26/01'; */

/* BEGIN MODULE TYPE.LINE */
typedef Char CHARARRAY[MAXLINE];

typedef struct LINE {
  CHARARRAY STR;
  long LEN;
} LINE;

/* END MODULE TYPE.LINE         VERSION= 'SUNMODS     Version  6/26/01'; */

/* RESTRICTION ENZYME */

typedef struct ENZYME {
  NUCLEOTIDE RSEQ[MAXPAT];
  long PLEN;
  WORD NAME, SITE;
  boolean SYMMETRIC;
  long CUT, ACUT;
} ENZYME;

/* BEGIN MODULE TYPE.FRAG */
/* LIST OF RESTRICTION SITES FOUND */

typedef struct FRAGMENT {
  long START, FINISH, SIZE;
  struct FRAGMENT *PREV, *NEXT;
} FRAGMENT;

typedef struct FRAGSFOUND {
  long LNUM;
  FRAGMENT *HEAD, *TAIL;
} FRAGSFOUND;

/* END MODULE TYPE.FRAG         VERSION= 'SUNMODS     Version  6/26/01'; */


Static _TEXT INFILE, OUTFILE;
Static LINE IFN, OFN;
Static WORD NAME;
Static NUCLEOTIDE COM[17];   /*Complements of NUCLEOTIDE */
Static SS NUCSET[17];
Static SS AMBIGUOUS;
Static NUCLEOTIDE ENDNUC;
Static ENZYME ENZ;   /*Enzyme site to search for */
Static FRAGSFOUND FOUND;   /*List of sites found*/
Static FRAGMENT *FREEFRAG;   /*Points to freelist of sites*/
Static SEQUENCE SEQ;   /*Sequence to be searched*/
Static long SEQLEN;   /*Sequence length*/
Static boolean CIRCULAR;   /*=true if seq. is circular */
Static LINE HLINE;
Static long I, CHOICE;


/* BEGIN MODULE INPLINE */
/*  Read a line from the console.    */
Static Void INPLINE(W_)
LINE *W_;
{
  long I;
  Char CH;
  Char BACKSPACE = '\b';

  /*!!!   if eoln(input) then readln; */
  /*!!!   get(input); */
  W_->LEN = 0;
  while (!P_eoln(stdin)) {
    CH = getchar();
    if (CH == '\n')
      CH = ' ';
    /*!!!       if not eoln then*/
    if (CH == BACKSPACE && W_->LEN > 0) {
      putchar(BACKSPACE);
      W_->LEN--;
    } else if (W_->LEN < MAXLINE) {
      W_->LEN++;
      W_->STR[W_->LEN - 1] = CH;
    }
  }
  scanf("%*[^\n]");
  getchar();
  for (I = W_->LEN; I < MAXLINE; I++)
    W_->STR[I] = ' ';
}  /* INPLINE */


/* END MODULE INPLINE         VERSION= 'SUNMODS     Version  6/26/01'; */

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


/* END MODULE WRITELINE         VERSION= 'SUNMODS     Version  6/26/01'; */

/* BEGIN MODULE GETFILE */
/* The c function access is found in sys/file.h. access must be declared in
   the outermost scope of the Pascal program. See access (2) manual pages */
/*!!!*/
extern long access PP((Char *PATH, long MODE));


#define F_OK            0   /* checks if file exists */


/* - - - - - - - - - - - - - - - - - - - - - - - - - - - */
/* =true if file exists, =false if it doesn't */
/* A comparable procedure is provided with ATT System V Pascal */
/* A special thank you to Mark Foster, CIS Rearch Computing, Univ. of Penn.
   and Glenn Gribble of Synaptics Inc. for showing me how to use access()
   within Pascal programs. B.F. */
Local boolean EXISTS(FILENAME)
LINE FILENAME;
{
  /* with FILENAME */
  FILENAME.STR[FILENAME.LEN] = '\0';   /* c strings terminate with null */
  if (access(FILENAME.STR, (long)F_OK) == -1)
    return false;
  else
    return true;
}  /* EXISTS */

#undef F_OK


/************************************************************/
/* Open  files, checking for I/O errors.                    */
/* Assumes user has already been prompted for filename.     */
/* Syntax: GETFILE(F,'I',<filename>);                       */
/*              opens file F for input                      */
/*         GETFILE(F,'O',<filename>);                       */
/*              opens file F for output                     */
/* !!! NOTE: this non-Standard procedure must be rewritten  */
/* for each different version of Pascal.                    */
/************************************************************/
Static Void GETFILE(F, FTYPE, FILENAME)
_TEXT *F;
Char FTYPE;
LINE *FILENAME;
{  /* GETFILE - - - - - - - - - - - - - - - - - - - - - */
  Char CH;
  boolean OKAY;

  /* Read in filename.*/
  INPLINE(FILENAME);
  switch (FTYPE) {

  /* Input file must exist before being opened.*/
  case 'I':  /* input file */
    while (!EXISTS(*FILENAME)) {
      printf("***FILE NOT FOUND\n");
      printf("Enter another filename:\n");
      INPLINE(FILENAME);
    }
    /*!!!*/
    strcpy(F->name, P_trimname(FILENAME->STR, sizeof(CHARARRAY)));
    if (F->f != NULL)
      F->f = freopen(F->name, "r", F->f);
    else
      F->f = fopen(F->name, "r");
    if (F->f == NULL)
      _EscIO2(FileNotFound, F->name);
    RESETBUF(F->f, Char);
    break;
    /* I */

  /* If output file already exists, ask user if he really wants
     to over write it. */
  case 'O':
    do {
      OKAY = true;
      if (EXISTS(*FILENAME)) {
	do {
	  printf("*** WARNING! File already exists. Overwrite?[Y|N]\n");
	  scanf("%c%*[^\n]", &CH);
	  getchar();
	  if (CH == '\n')
	    CH = ' ';
	} while (CH != 'n' && CH != 'N' && CH != 'y' && CH != 'Y');
	if (CH != 'y' && CH != 'Y') {
	  printf("Enter another filename:\n");
	  INPLINE(FILENAME);
	  OKAY = false;
	}
      }  /* EXISTS(FILENAME) */
    } while (!OKAY);
    strcpy(F->name, P_trimname(FILENAME->STR, sizeof(CHARARRAY)));
    if (F->f != NULL)
      F->f = freopen(F->name, "w", F->f);
    else
      F->f = fopen(F->name, "w");
    if (F->f == NULL)
      _EscIO2(FileNotFound, F->name);
    SETUPBUF(F->f, Char);   /* output file */
    break;
    /* O */
  }
}  /* GETFILE */


/* END MODULE GETFILE         VERSION= 'SUNMODS     Version  6/26/01'; */


/**************************************************/
/*  WORD   I/O  PROCEDURES                        */
/**************************************************/
/* BEGIN MODULE INPWORD */
/* Read a WORD from the terminal. */
Static Void INPWORD(W_)
WORD *W_;
{
  long I;
  Char CH;
  Char BACKSPACE = '\b';

  W_->LEN = 0;
  while (P_peek(stdin) == ' ') {
    CH = getchar();
    if (CH == '\n')
      CH = ' ';
  }
  while ((P_peek(stdin) != ' ') & (!P_eoln(stdin))) {
    CH = getchar();
    if (CH == '\n')
      CH = ' ';
    /*!!!       if not eoln then*/
    if (CH == BACKSPACE && W_->LEN > 0) {
      putchar(BACKSPACE);
      W_->LEN--;
    } else if (W_->LEN < MAXWORD) {
      W_->LEN++;
      W_->STR[W_->LEN - 1] = CH;
    }
  }
  /*readln;*/
  for (I = W_->LEN; I < MAXWORD; I++)
    W_->STR[I] = ' ';
}  /* INPWORD */


/* END MODULE INPWORD         VERSION= 'SUNMODS     Version  6/26/01'; */

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


/* END MODULE READWORD         VERSION= 'SUNMODS     Version  6/26/01'; */

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


/* END MODULE WRITEWORD         VERSION= 'SUNMODS     Version  6/26/01'; */

/************************************************/
/*  Initialization procedures.                  */
/************************************************/
Static Void INITIALIZE()
{
  NUCLEOTIDE B1 = A, B2 = T;

  /* Array NUCSET holds sets of nucleotides  */
  /* using the conventions of IUPAC-IUB      */
  NUCSET[(long)A] = 1L << ((long)A);
  NUCSET[(long)C] = 1L << ((long)C);
  NUCSET[(long)G] = 1L << ((long)G);
  NUCSET[(long)T] = 1L << ((long)T);
  NUCSET[(long)R] = (1L << ((long)A)) | (1L << ((long)G)) | (1L << ((long)R));
  NUCSET[(long)Y] = (1L << ((long)C)) | (1L << ((long)T)) | (1L << ((long)Y));
  NUCSET[(long)S] = (1L << ((long)C)) | (1L << ((long)G)) | (1L << ((long)S));
  NUCSET[(long)M] = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)M));
  NUCSET[(long)K] = (1L << ((long)G)) | (1L << ((long)T)) | (1L << ((long)K));
  NUCSET[(long)W] = (1L << ((long)A)) | (1L << ((long)T)) | (1L << ((long)W));
  NUCSET[(long)B] = (1L << ((long)C)) | (1L << ((long)G)) |
		    (1L << ((long)T)) | (1L << ((long)B)) |
		    (1L << ((long)Y)) | (1L << ((long)K)) | (1L << ((long)S));
  NUCSET[(long)D] = (1L << ((long)A)) | (1L << ((long)G)) |
		    (1L << ((long)T)) | (1L << ((long)W)) |
		    (1L << ((long)R)) | (1L << ((long)D)) | (1L << ((long)K));
  NUCSET[(long)H] = (1L << ((long)A)) | (1L << ((long)C)) |
		    (1L << ((long)T)) | (1L << ((long)M)) |
		    (1L << ((long)Y)) | (1L << ((long)W)) | (1L << ((long)H));
  NUCSET[(long)V] = (1L << ((long)A)) | (1L << ((long)C)) |
		    (1L << ((long)G)) | (1L << ((long)M)) |
		    (1L << ((long)R)) | (1L << ((long)S)) | (1L << ((long)V));
  NUCSET[(long)N] = (1L << ((long)N + 1)) - (1L << ((long)A));
  AMBIGUOUS = NUCSET[(long)N] & (~((1L << ((long)A)) | (1L << ((long)C)) |
				   (1L << ((long)G)) | (1L << ((long)T))));
  /* Array  COM holds the symbol values of nucleotide complements.*/
  /* A complements T, C complements G, R complements Y etc.       */
  while ((long)B1 <= (long)T) {
    COM[(long)B1] = B2;
    B1 = (NUCLEOTIDE)((long)B1 + 1);
    B2 = (NUCLEOTIDE)((long)B2 - 1);
  }
  /* N,S, and W complement themselves */
  for (B1 = W; (long)B1 <= (long)N; B1 = (NUCLEOTIDE)((long)B1 + 1))
    COM[(long)B1] = B1;
  /* Initialize the linked list FOUND */
  FREEFRAG = NULL;
  FOUND.HEAD = (FRAGMENT *)Malloc(sizeof(FRAGMENT));
  FOUND.TAIL = (FRAGMENT *)Malloc(sizeof(FRAGMENT));
  FOUND.HEAD->NEXT = FOUND.TAIL;
  FOUND.TAIL->PREV = FOUND.HEAD;
}  /* INITIALIZE */


/*******************************************************************/
/* If sequence is circular, copy the last MAXPAT bases to the start*/
/*******************************************************************/
Static Void INITENDS()
{
  long Km = 1 - MAXPAT;
  long Sk, FORLIM;

  if (!CIRCULAR)
    return;
  FORLIM = SEQLEN;
  for (Sk = SEQLEN - MAXPAT + 1; Sk <= FORLIM; Sk++) {
    SEQ[Km + MAXPAT] = SEQ[Sk + MAXPAT];
    Km++;
  }
}  /* INITENDS */


/* BEGIN MODULE GETINTEGER */
/* Prompts user for an integer, checks whether all characters are digits,*/
/* and whether number is within desired range; harasses user until valid */
/* integer is received. */
Static Void GETINTEGER(NUM, LBOUND, HBOUND)
long *NUM, LBOUND, HBOUND;
{
  long I, VAL;
  long ORDZERO = '0';
  boolean LEGAL = false, INRANGE = false, NEGATIVE = false;
  Char CH;
  WORD NUMWORD;

  do {
    do {
      INPWORD(&NUMWORD);

      /* Evaluate sign, if any */
      I = 1;
      CH = NUMWORD.STR[I-1];
      if (CH == '-') {
	NEGATIVE = true;
	I++;
	CH = NUMWORD.STR[I-1];
      } else {
	NEGATIVE = false;
	if (CH == '+') {
	  I++;
	  CH = NUMWORD.STR[I-1];
	}
      }

      /* Evaluate unsigned integer */
      *NUM = 0;
      while (isdigit(CH)) {
	VAL = CH - ORDZERO;
	*NUM = *NUM * 10 + VAL;
	I++;
	CH = NUMWORD.STR[I-1];
      }
      if (I > NUMWORD.LEN)
	LEGAL = true;
      else {
	LEGAL = false;
	printf("\nIllegal character encountered.\n");
	printf("Enter new value:  \n");
      }
    } while (!LEGAL);

    /* If the number entered was negative, multiply */
    /* NUM by -1. Check range of number.            */
    if (NEGATIVE)
      *NUM = -*NUM;
    if (*NUM >= LBOUND && *NUM <= HBOUND)
      INRANGE = true;
    else {
      INRANGE = false;
      printf("\nNumber is out of range.\n");
      printf("Please enter new value:\n");
    }
  } while (!INRANGE);
}  /* GETINTEGER */


/* END MODULE GETINTEGER         VERSION= 'SUNMODS     Version  6/26/01'; */

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


/* END MODULE NUC         VERSION= 'SUNMODS     Version  6/26/01'; */

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


/* END MODULE SAMEWORD         VERSION= 'SUNMODS     Version  6/26/01'; */

/* BEGIN MODULE READSEQ */
/*  ******************************************* */
/*  Read a DNA or RNA sequence from SFILE       */
/*   and store it in S.                         */
/*  ******************************************* */
Static Void READSEQ(SFILE, S_, SEQLEN, NAME, CIRCULAR)
_TEXT *SFILE;
NUCLEOTIDE *S_;
long *SEQLEN;
WORD *NAME;
boolean *CIRCULAR;
{
  Char CH, FILETYPE, LASTNUM;
  boolean TOOLONG = false, BADFILE = false;
  WORD ID, LOCUS, ORI, CIRC;

  /* Prompt user for sequence file type */
  printf("The following file formats can be read:\n");
  printf("  F:free format   B:BIONET    G:GENBANK\n");
  do {
    printf("Type letter of format (F|B|G)\n");
    scanf("%c%*[^\n]", &CH);
    getchar();
    if (CH == '\n')
      CH = ' ';
  } while (CH != 'g' && CH != 'G' && CH != 'b' && CH != 'B' && CH != 'f' &&
	   CH != 'F');
  switch (CH) {

  case 'F':
  case 'f':
    FILETYPE = 'F';
    break;

  case 'B':
  case 'b':
    FILETYPE = 'B';
    break;

  case 'G':
  case 'g':
    FILETYPE = 'G';
    break;
  }
  printf("Reading input file...\n");
  NAME->LEN = 0;

  /* For BIONET or GENBANK, read in sequence name and topology */
  /* Advance to beginning of sequence */
  if (FILETYPE == 'B') {
    /* First non-comment line is name. Name may be blank.*/
    while (P_peek(SFILE->f) == ';') {
      fscanf(SFILE->f, "%*[^\n]");
      getc(SFILE->f);
    }
    while ((P_peek(SFILE->f) == ' ') & (!P_eoln(SFILE->f)))
      getc(SFILE->f);
    if (!P_eoln(SFILE->f))
      READWORD(SFILE, NAME);
    else
      NAME->LEN = 0;
    if (!BUFEOF(SFILE->f)) {
      fscanf(SFILE->f, "%*[^\n]");
      getc(SFILE->f);
    }
  }  /* BIONET */

  else if (FILETYPE == 'G') {
    /* Initialize identifiers */
    LOCUS.STR[0] = 'L';
    LOCUS.STR[1] = 'O';
    LOCUS.STR[2] = 'C';
    LOCUS.STR[3] = 'U';
    LOCUS.STR[4] = 'S';
    LOCUS.LEN = 5;
    ORI.STR[0] = 'O';
    ORI.STR[1] = 'R';
    ORI.STR[2] = 'I';
    ORI.STR[3] = 'G';
    ORI.STR[4] = 'I';
    ORI.STR[5] = 'N';
    ORI.LEN = 6;
    CIRC.STR[0] = 'c';
    CIRC.STR[1] = 'i';
    CIRC.STR[2] = 'r';
    CIRC.STR[3] = 'c';
    CIRC.STR[4] = 'u';
    CIRC.STR[5] = 'l';
    CIRC.STR[6] = 'a';
    CIRC.STR[7] = 'r';
    CIRC.LEN = 8;
    /* Advance to LOCUS line. Read in NAME and topology. */
    while (!((P_peek(SFILE->f) == 'L') | BUFEOF(SFILE->f))) {
      fscanf(SFILE->f, "%*[^\n]");
      getc(SFILE->f);
    }
    if (!BUFEOF(SFILE->f)) {
      READWORD(SFILE, &ID);
      if (SAMEWORD(&ID, &LOCUS)) {
	if (!BUFEOF(SFILE->f))
	  READWORD(SFILE, NAME);
	if (!BUFEOF(SFILE->f))   /* skip seq. length */
	  READWORD(SFILE, &ID);
	if (!BUFEOF(SFILE->f))   /* skip. 'bp' */
	  READWORD(SFILE, &ID);
	/* After 'bp', there is an optional field telling
	   the type of molecule (ss-RNA, ds-DNA etc.)
	   Since this field is optional, we must test the
	   next two tokens to see if they are 'circular' */
	*CIRCULAR = false;
	if (!BUFEOF(SFILE->f))
	  READWORD(SFILE, &ID);
	if (SAMEWORD(&ID, &CIRC))
	  *CIRCULAR = true;
	else {
	  if (!BUFEOF(SFILE->f))
	    READWORD(SFILE, &ID);
	  if (SAMEWORD(&ID, &CIRC))
	    *CIRCULAR = true;
	}
      }  /* SAMEWORD(ID,LOCUS) */
      else
	BADFILE = true;
    }

    /* Advance to ORIGIN line. Sequence begins on next line */
    if (!BUFEOF(SFILE->f)) {
      do {
	fscanf(SFILE->f, "%*[^\n]");
	getc(SFILE->f);
	if (P_peek(SFILE->f) == 'O')
	  READWORD(SFILE, &ID);
      } while (!(SAMEWORD(&ID, &ORI) | BUFEOF(SFILE->f)));
      if (SAMEWORD(&ID, &ORI)) {
	fscanf(SFILE->f, "%*[^\n]");
	getc(SFILE->f);
      }
    }
    if (BUFEOF(SFILE->f))
      BADFILE = true;
  }  /* GENBANK */

  /* Read in sequence */
  *SEQLEN = 0;
  if (BADFILE) {
    printf(">>> ERROR: Not a GENBANK file. No sequence read.\n");
    /*!!!   CLOSE(SFILE);*/
    return;
  }
  while (!BUFEOF(SFILE->f)) {
    while (!P_eoln(SFILE->f)) {
      CH = getc(SFILE->f);
      if (CH == '\n')
	CH = ' ';
      if (CH == 'b' || CH == 'B' || CH == 'v' || CH == 'V' || CH == 'h' ||
	  CH == 'H' || CH == 'd' || CH == 'D' || CH == 'k' || CH == 'K' ||
	  CH == 's' || CH == 'S' || CH == 'w' || CH == 'W' || CH == 'm' ||
	  CH == 'M' || CH == 'y' || CH == 'Y' || CH == 'r' || CH == 'R' ||
	  CH == 'n' || CH == 'N' || CH == 'u' || CH == 'U' || CH == 't' ||
	  CH == 'T' || CH == 'g' || CH == 'G' || CH == 'c' || CH == 'C' ||
	  CH == 'a' || CH == 'A') {
/* p2c: intrest.p, line 528: 
 * Note: Line breaker spent 0.0 seconds, 5000 tries on line 784 [251] */
	if (*SEQLEN < MAXSEQ - 2) {
	  /*write(CH);*/
	  (*SEQLEN)++;
	  S_[*SEQLEN + MAXPAT] = NUC(CH);
	} else
	  TOOLONG = true;
	continue;
      }
      if (CH == ';') {   /*begin    comment in input file */
	fscanf(SFILE->f, "%*[^\n]");
	getc(SFILE->f);
      }
      /*;writeln end*/
      else if (CH == '2' || CH == '1')
	LASTNUM = CH;
    }
    fscanf(SFILE->f, "%*[^\n]");
    getc(SFILE->f);   /*writeln*/
  }
  if (TOOLONG)
    printf(">>> WARNING! Sequence length exceeds MAXSEQ-2. Seq. truncated.\n");

  if (FILETYPE == 'F') {
    do {
      printf("Is sequence circular or linear? (Type C or L)\n");
      scanf("%c%*[^\n]", &CH);
      getchar();
      if (CH == '\n')
	CH = ' ';
    } while (CH != 'l' && CH != 'L' && CH != 'c' && CH != 'C');
    switch (CH) {

    case 'C':
    case 'c':
      *CIRCULAR = true;
      break;

    case 'L':
    case 'l':
      *CIRCULAR = false;
      break;
    }
    return;
  }
  if (FILETYPE != 'B')
    return;
  if (LASTNUM == '1')
    *CIRCULAR = false;
  else
    *CIRCULAR = true;
}  /* READSEQ */


/* END MODULE READSEQ         VERSION= 'SUNMODS     Version  6/26/01'; */


/*********************************************************/
/*  Read in a restriction enzyme site from the terminal. */
/*********************************************************/
Static Void READSITE(ENZ)
ENZYME *ENZ;
{
  long LOC;
  boolean LEGAL;
  WORD *WITH1;
  long FORLIM;

  /* with ENZ */
  do {
    LEGAL = true;
    printf("Enter name of enzyme:\n");
    INPWORD(&ENZ->NAME);
    scanf("%*[^\n]");
    getchar();
    printf("Enter recognition site using the symbols below:\n");
    printf("Nucleotides:  A,C,G,T\n");
    printf("Ambiguous nucleotides:\n");
    printf("R = A,G    M = A,C    B = C,G,T   N = A,C,G,T\n");
    printf("Y = C,T    K = G,T    D = A,G,T   V = A,C,G\n");
    printf("S = C,G    W = A,T    H = A,C,T\n\n");
    printf("Site must be <=%12ld characters:\n", MAXPAT - 1L);
    INPWORD(&ENZ->SITE);
    scanf("%*[^\n]");
    getchar();
    /* Check for illegal symbols */
    WITH1 = &ENZ->SITE;
    FORLIM = WITH1->LEN;
    for (LOC = 0; LOC < FORLIM; LOC++) {
      if (WITH1->STR[LOC] == 'w' || WITH1->STR[LOC] == 'W' ||
	  WITH1->STR[LOC] == 's' || WITH1->STR[LOC] == 'S' ||
	  WITH1->STR[LOC] == 'n' || WITH1->STR[LOC] == 'N' ||
	  WITH1->STR[LOC] == 't' || WITH1->STR[LOC] == 'T' ||
	  WITH1->STR[LOC] == 'g' || WITH1->STR[LOC] == 'G' ||
	  WITH1->STR[LOC] == 'y' || WITH1->STR[LOC] == 'Y' ||
	  WITH1->STR[LOC] == 'h' || WITH1->STR[LOC] == 'H' ||
	  WITH1->STR[LOC] == 'b' || WITH1->STR[LOC] == 'B' ||
	  WITH1->STR[LOC] == 'k' || WITH1->STR[LOC] == 'K' ||
	  WITH1->STR[LOC] == 'm' || WITH1->STR[LOC] == 'M' ||
	  WITH1->STR[LOC] == 'v' || WITH1->STR[LOC] == 'V' ||
	  WITH1->STR[LOC] == 'd' || WITH1->STR[LOC] == 'D' ||
	  WITH1->STR[LOC] == 'r' || WITH1->STR[LOC] == 'R' ||
	  WITH1->STR[LOC] == 'c' || WITH1->STR[LOC] == 'C' ||
	  WITH1->STR[LOC] == 'a' || WITH1->STR[LOC] == 'A')
/* p2c: intrest.p, line 569: 
 * Note: Line breaker spent 0.0 seconds, 5000 tries on line 889 [251] */
	ENZ->RSEQ[LOC] = NUC(WITH1->STR[LOC]);
      else {
	printf("Illegal character  %c\n", WITH1->STR[LOC]);
	LEGAL = false;
      }
    }
  } while (!LEGAL);
  ENZ->PLEN = ENZ->SITE.LEN;
  printf("Position of cut?\n");
  GETINTEGER(&ENZ->CUT, -LONG_MAX, LONG_MAX);
  scanf("%*[^\n]");
  getchar();
}  /* READSITE */


/* BEGIN MODULE LINKED */
/*********************************************************/
/*  Linked-list operations for restriction fragment list.*/
/*********************************************************/

/*Get a new fragment from freelist.*/
Static Void GETFRAG(NEWFRAG)
FRAGMENT **NEWFRAG;
{
  if (FREEFRAG == NULL)
    *NEWFRAG = (FRAGMENT *)Malloc(sizeof(FRAGMENT));
  else {
    *NEWFRAG = FREEFRAG;
    FREEFRAG = FREEFRAG->NEXT;
  }
}


/*Add a fragment after DEST*/
Static Void ADDFRAG(AFRAG, DEST)
FRAGMENT **AFRAG, **DEST;
{
  FRAGMENT *TEMP;

  TEMP = (*DEST)->NEXT;
  (*DEST)->NEXT = *AFRAG;
  (*AFRAG)->PREV = *DEST;
  (*AFRAG)->NEXT = TEMP;
  TEMP->PREV = *AFRAG;
}


/*Return a list to the top of freelist*/
Static Void RIDOF(HEAD, TAIL)
FRAGMENT **HEAD, **TAIL;
{
  FRAGMENT *TEMPHEAD, *TEMPTAIL;

  if ((*HEAD)->NEXT == *TAIL)
    return;
  TEMPHEAD = (*HEAD)->NEXT;
  TEMPTAIL = (*TAIL)->PREV;
  (*HEAD)->NEXT = *TAIL;
  (*TAIL)->PREV = *HEAD;
  TEMPHEAD->PREV = NULL;
  TEMPTAIL->NEXT = FREEFRAG;
  FREEFRAG = TEMPHEAD;
}


/* Local variables for KMP: */
struct LOC_KMP {
  ENZYME *ENZ;
  long NEXT[MAXPAT], GOBACK[MAXPAT];
  SS PATTERN[MAXPAT];
  long PLEN1;
  FRAGMENT *RIGHTMOST;
} ;

/* The NEXT array is used when  PATTERN[Px] does not match  STR[Sk].*/
/* After a mismatch is found, instead of backing up to PATTERN[1],  */
/* the search resumes at PATTERN[NEXT[Pj]].                         */
Local Void MAKENEXT(LINK)
struct LOC_KMP *LINK;
{
  long Pj;
  long FPj = 0;
  long LEFTAM;
  ENZYME *WITH;
  long FORLIM;

  WITH = LINK->ENZ;
  LINK->PLEN1 = WITH->PLEN + 1;
  /* Ambiguous nucleotide in first position is a special case */
  if ((AMBIGUOUS & LINK->PATTERN[0]) != 0) {
    LINK->NEXT[0] = 0;
    LINK->GOBACK[0] = 0;
    FORLIM = LINK->PLEN1;
    for (Pj = 2; Pj <= FORLIM; Pj++) {
      LINK->NEXT[Pj-1] = 1;
      LINK->GOBACK[Pj-1] = Pj - 2;
    }
    return;
  }
  /*Set LEFTAM equal to position of leftmost ambiguous symbol*/
  /* after the first position in pattern.                    */
  LEFTAM = LINK->PLEN1;
  Pj = 2;
  while (Pj < LINK->PLEN1) {
    if ((AMBIGUOUS & LINK->PATTERN[Pj-1]) != 0) {
      LEFTAM = Pj;
      goto _L1;
    }
    Pj++;
  }
  /*Calculate NEXT & GOBACK tables*/
_L1:
  Pj = 1;
  LINK->NEXT[0] = 0;
  LINK->GOBACK[0] = 0;
  while (Pj < LINK->PLEN1) {  /*FPj = F[Pj] */
_L2:
    if (FPj > 0) {
      if ((LINK->PATTERN[Pj-1] & (~LINK->PATTERN[FPj-1])) != 0) {
	FPj = LINK->NEXT[FPj-1];
	goto _L2;
      }
    }
    FPj++;
    Pj++;
    if (Pj <= LEFTAM) {
      if ((LINK->PATTERN[Pj-1] & (~LINK->PATTERN[FPj-1])) == 0)
	LINK->NEXT[Pj-1] = LINK->NEXT[FPj-1];
      else
	LINK->NEXT[Pj-1] = FPj;
      LINK->GOBACK[Pj-1] = 0;
    } else {
      LINK->NEXT[Pj-1] = LINK->NEXT[LEFTAM-1];
      LINK->GOBACK[Pj-1] = Pj - LEFTAM;
    }
  }
  LINK->GOBACK[LINK->PLEN1-1] = Pj - LEFTAM;
}  /*MAKENEXT*/

/* Local variables for SEARCH: */
struct LOC_SEARCH {
  struct LOC_KMP *LINK;
  long CUTSITE, Sk;
} ;

/*  Add a node to the linked list FOUND, telling the cutting site.*/
Local Void PATTERNMATCHED(LINK)
struct LOC_SEARCH *LINK;
{
  FRAGMENT *NEWFRAG;
  long POSN;
  ENZYME *WITH;

  WITH = LINK->LINK->ENZ;
  POSN = LINK->Sk - WITH->PLEN + LINK->CUTSITE +
	 LINK->LINK->GOBACK[LINK->LINK->PLEN1-1];
  if (CIRCULAR) {
    if (POSN < 1)
      POSN += SEQLEN;
    else if (POSN > SEQLEN)
      POSN -= SEQLEN;
  }
  if (POSN < 1 || POSN > SEQLEN)
    return;
  GETFRAG(&NEWFRAG);
  NEWFRAG->START = POSN;
  if (LINK->LINK->RIGHTMOST->START > POSN)
    LINK->LINK->RIGHTMOST = FOUND.HEAD;
  while (LINK->LINK->RIGHTMOST->NEXT->START < POSN)
    LINK->LINK->RIGHTMOST = LINK->LINK->RIGHTMOST->NEXT;
  ADDFRAG(&NEWFRAG, &LINK->LINK->RIGHTMOST);
  LINK->LINK->RIGHTMOST = NEWFRAG;
  FOUND.LNUM++;
}  /* PATTERNMATCHED */

/* Search the sequence.*/
Local Void SEARCH(CUTSITE_, LINK)
long CUTSITE_;
struct LOC_KMP *LINK;
{
  struct LOC_SEARCH V;
  long Pj = 1;
  long RESUME;
  SS ENDSET;
  ENZYME *WITH;

  V.LINK = LINK;
  V.CUTSITE = CUTSITE_;
  WITH = LINK->ENZ;
  /* Initialize search parameters */
  ENDSET = LINK->PATTERN[0];
  SEQ[SEQLEN + MAXPAT + 1] = WONT;   /*WONT will never be found in CANT set*/
  RESUME = LINK->NEXT[WITH->PLEN];
      /* Resume search at this position after match*/
  LINK->NEXT[WITH->PLEN] = -1;
  SEQ[SEQLEN + MAXPAT + 2] = ENDNUC;   /*BASE value of PATTERN[1]*/
  /*Begin at leftmost chars of PATTERN and STR*/
  if (CIRCULAR)
    V.Sk = 2 - WITH->PLEN;
  else
    V.Sk = 1;

_L10:   /* GETSTARTED, Pj = 1 */
  /*Advance until first match */
  while (((1L << ((long)SEQ[V.Sk + MAXPAT])) & ENDSET) == 0)
    V.Sk++;
  if (V.Sk > SEQLEN)   /*INPUTEXHAUSTED*/
    goto _L99;

_L20:   /* CHARMATCHED */
  Pj++;
  V.Sk++;

_L30:   /* LOOP, Pj > 0 */
  if (((1L << ((long)SEQ[V.Sk + MAXPAT])) & LINK->PATTERN[Pj-1]) != 0)
	/*CHARMATCHED*/
	  goto _L20;

  V.Sk -= LINK->GOBACK[Pj-1];
  Pj = LINK->NEXT[Pj-1];   /* char not matched */
  if (Pj == 1)   /*GETSTARTED*/
    goto _L10;
  if (Pj == 0) {
    Pj = 1;
    V.Sk++;
    goto _L10;   /*GETSTARTED*/
  }
  if (Pj > 0)   /*LOOP*/
    goto _L30;

  /* Pj = -1, PATTERN matched */
  PATTERNMATCHED(&V);
  Pj = RESUME;
  goto _L30;   /*LOOP*/

_L99: ;   /* INPUTEXHAUSTED */
}  /*SEARCH*/


/* END MODULE LINKED         VERSION= 'SUNMODS     Version  6/26/01'; */

/****************************************************************/
/*  Search for site ENZ in SEQ, using the algorithm found in    */
/*    "Fast Pattern Matching in Strings"                        */
/*     Knuth, Morris, and Pratt, SIAM Journal of Computing,     */
/*     Vol.6, No.2, June 1977.  pp323-350.                      */
/****************************************************************/
Static Void KMP(ENZ_, OUTFILE)
ENZYME *ENZ_;
_TEXT *OUTFILE;
{
  struct LOC_KMP V;
  long Pj, Cj;
  NUCLEOTIDE CRP;
  ENZYME *WITH;
  long FORLIM;


  V.ENZ = ENZ_;
  FOUND.LNUM = 0;
  FOUND.TAIL->START = MAXSEQ + 1;   /*Always > POSN of site found*/
  WITH = V.ENZ;

  /* with ENZ */
  /* Search SEQ using PATTERN derived from RSEQ  */
  ENDNUC = WITH->RSEQ[0];
  FORLIM = WITH->PLEN;
  for (Pj = 0; Pj < FORLIM; Pj++)
    V.PATTERN[Pj] = NUCSET[(long)WITH->RSEQ[Pj]];
  V.PATTERN[WITH->PLEN] = 1L << ((long)CANT);   /* Can not be matched */
  V.RIGHTMOST = FOUND.HEAD;
  MAKENEXT(&V);
  SEARCH(WITH->CUT, &V);

  /* Let PATTERN = the inverse complement of PATTERN*/
  /* Check for symmetry.                            */
  Cj = WITH->PLEN;
  WITH->SYMMETRIC = true;
  FORLIM = WITH->PLEN;
  for (Pj = 0; Pj < FORLIM; Pj++) {
    CRP = COM[(long)WITH->RSEQ[Pj]];
    V.PATTERN[Cj-1] = NUCSET[(long)CRP];
    if (WITH->RSEQ[Cj-1] != CRP)
      WITH->SYMMETRIC = false;
    Cj--;
  }
  if (WITH->SYMMETRIC) {
    return;
  }  /* not SYMMETRIC */
  printf("Asymmetric Site!\n");
  printf("Type position of cut on opposite strand\n");
  GETINTEGER(&WITH->ACUT, -LONG_MAX, LONG_MAX);
  scanf("%*[^\n]");
  getchar();
  ENDNUC = CRP;
  MAKENEXT(&V);
  V.RIGHTMOST = FOUND.HEAD;
  SEARCH(WITH->PLEN - WITH->ACUT, &V);
}  /*KMP*/


/* Local variables for REPORT: */
struct LOC_REPORT {
  _TEXT *OUTFILE;
  FRAGMENT *ORDER[MAXFRAGS];
  boolean TOOMANY_FRAGS;   /* =true if SITES > MAXFRAGS */
  long NUMSITES;
} ;

/*Calculate sizes and ends of fragments*/
Local Void CALCULATE(LINK)
struct LOC_REPORT *LINK;
{
  FRAGMENT *CURRENTFRAG;
  long SITE = 1;

  LINK->NUMSITES = FOUND.LNUM;
  /* if LINEAR add a fragment to head of list*/
  if (!CIRCULAR && FOUND.HEAD->NEXT->START > 1) {
    GETFRAG(&CURRENTFRAG);
    ADDFRAG(&CURRENTFRAG, &FOUND.HEAD);
    CURRENTFRAG->START = 1;
    FOUND.LNUM++;
  } else
    CURRENTFRAG = FOUND.HEAD->NEXT;

  if (FOUND.LNUM <= 0)
    return;
  /*Calculate ends and size of each fragment and assign it */
  /* a place in the ORDER array, to be sorted later.    */
  while (CURRENTFRAG->NEXT != FOUND.TAIL) {
    CURRENTFRAG->FINISH = CURRENTFRAG->NEXT->START - 1;
    CURRENTFRAG->SIZE = CURRENTFRAG->FINISH - CURRENTFRAG->START + 1;
    if (SITE <= MAXFRAGS)
      LINK->ORDER[SITE-1] = CURRENTFRAG;
    else
      LINK->TOOMANY_FRAGS = true;
    SITE++;
    CURRENTFRAG = CURRENTFRAG->NEXT;
  }
  /*Last fragment in list is a special case*/
  if (CIRCULAR) {
    CURRENTFRAG->FINISH = FOUND.HEAD->NEXT->START - 1;
    CURRENTFRAG->SIZE = SEQLEN - CURRENTFRAG->START + CURRENTFRAG->FINISH + 1;
    if (CURRENTFRAG->FINISH == 0)
      CURRENTFRAG->FINISH = SEQLEN;
  } else {
    CURRENTFRAG->FINISH = SEQLEN;
    CURRENTFRAG->SIZE = CURRENTFRAG->FINISH - CURRENTFRAG->START + 1;
  }
  if (SITE <= MAXFRAGS)
    LINK->ORDER[SITE-1] = CURRENTFRAG;
  else
    LINK->TOOMANY_FRAGS = true;
}  /* CALCULATE */

Local Void SWAP(FIRST, SECOND)
FRAGMENT **FIRST, **SECOND;
{
  FRAGMENT *TEMP;

  TEMP = *FIRST;
  *FIRST = *SECOND;
  *SECOND = TEMP;
}

/* BEGIN MODULE SORT */
/*  Invariant:  The array elements > TOP are sorted.*/
/*    TOP >= unsorted elements.                     */
Local Void BUBBLESORT(TOP, BOTTOM, LINK)
long TOP, BOTTOM;
struct LOC_REPORT *LINK;
{
  long SMALLEST, NEXT;

  while (TOP >= BOTTOM) {
    /*bubble smallest unsorted number to the top of sorted list*/
    SMALLEST = BOTTOM;
    NEXT = BOTTOM + 1;
    while (NEXT <= TOP) {
      if (LINK->ORDER[SMALLEST-1]->SIZE < LINK->ORDER[NEXT-1]->SIZE)
	SWAP(&LINK->ORDER[SMALLEST-1], &LINK->ORDER[NEXT-1]);
      SMALLEST = NEXT;
      NEXT = SMALLEST + 1;
    }
    TOP--;
  }
}  /* BUBBLESORT */

/* END MODULE SORT         VERSION= 'SUNMODS     Version  6/26/01'; */

/*  Print the list of sites found.        */
Local Void PRINTLIST(LINK)
struct LOC_REPORT *LINK;
{
  FRAGMENT *THISFRAG;
  long SITE = 1;
  FRAGMENT *WITH;

  /* Write a heading for the site */
  WRITEWORD(LINK->OUTFILE, ENZ.NAME, 10L);
  WRITEWORD(LINK->OUTFILE, ENZ.SITE, (long)MAXPAT);
  fprintf(LINK->OUTFILE->f, "%3ld", ENZ.CUT);
  if (ENZ.SYMMETRIC)
    fprintf(LINK->OUTFILE->f, "%6c", ' ');
  else
    fprintf(LINK->OUTFILE->f, " (%3ld)", ENZ.ACUT);
  fprintf(LINK->OUTFILE->f, "%5ld\n", LINK->NUMSITES);
  if (LINK->NUMSITES == FOUND.LNUM)
    THISFRAG = FOUND.HEAD;
  else
    THISFRAG = FOUND.HEAD->NEXT;
  while (THISFRAG->NEXT != FOUND.TAIL) {
    THISFRAG = THISFRAG->NEXT;
    fprintf(LINK->OUTFILE->f, "%*ld", (int)(MAXPAT + 30), THISFRAG->START);
    WITH = LINK->ORDER[SITE-1];
    fprintf(LINK->OUTFILE->f, "%8ld%8ld%8ld\n",
	    WITH->SIZE, WITH->START, WITH->FINISH);
    SITE++;
  }
  /* Print last fragment. If no sites are found, a linear */
  /* sequence will still yeild one fragment.              */
  if (LINK->NUMSITES < FOUND.LNUM) {
    WITH = LINK->ORDER[SITE-1];
    fprintf(LINK->OUTFILE->f, "%*ld%8ld%8ld\n",
	    (int)(MAXPAT + 38), WITH->SIZE, WITH->START, WITH->FINISH);
  }
  putc('\n', LINK->OUTFILE->f);
}  /* PRINTLIST */


/*********************************************/
/* Compile a report for output.              */
/*********************************************/
Static Void REPORT(OUTFILE_)
_TEXT *OUTFILE_;
{
  struct LOC_REPORT V;


  V.OUTFILE = OUTFILE_;
  V.TOOMANY_FRAGS = false;
  CALCULATE(&V);
  if (V.TOOMANY_FRAGS) {
    printf(">>> Number of fragments exceeds %12ld\n", (long)MAXFRAGS);
    printf(">>> Increase MAXFRAGS and recompile\n");
    fprintf(V.OUTFILE->f, "\n>>> Number of fragments exceeds %12ld\n",
	    (long)MAXFRAGS);
    fprintf(V.OUTFILE->f, ">>> Increase MAXFRAGS and recompile\n");
  } else {
    BUBBLESORT(FOUND.LNUM, 1L, &V);
    PRINTLIST(&V);
  }
  RIDOF(&FOUND.HEAD, &FOUND.TAIL);
}  /*REPORT*/


/***************************************************************/
/* Prompt user for sites and call SEARCH.                      */
/***************************************************************/
Static Void PROMPT(OUTFILE)
_TEXT *OUTFILE;
{
  Char ANSWER;

  /* Write heading */
  fprintf(OUTFILE->f, "%s\n", VERSION);
  WRITEWORD(OUTFILE, NAME, NAME.LEN);
  fprintf(OUTFILE->f, "  Configuration: ");
  switch (CIRCULAR) {

  case true:
    fprintf(OUTFILE->f, " CIRCULAR");
    break;

  case false:
    fprintf(OUTFILE->f, " LINEAR");
    break;
  }
  fprintf(OUTFILE->f, "  Length: %12ld bp\n", SEQLEN);
  fprintf(OUTFILE->f, "%*s\n", (int)(MAXPAT + 24), "# of");
  fprintf(OUTFILE->f, "%*s\n",
	  (int)(MAXPAT + 54), "Cut     Sites Sites   Frags   Begin     End");
  do {
    READSITE(&ENZ);
    KMP(&ENZ, OUTFILE);
    REPORT(OUTFILE);
    do {
      printf("Type  S to search for a site, Q to return to main menu:\n");
      scanf("%c%*[^\n]", &ANSWER);
      getchar();
      if (ANSWER == '\n')
	ANSWER = ' ';
    } while (ANSWER != 'q' && ANSWER != 'Q' && ANSWER != 's' && ANSWER != 'S');
  } while (ANSWER != 'q' && ANSWER != 'Q');
}  /* PROMPT */


/* ----------------------------------------------------------  */
/* ----------------- MAIN  PROCEDURE  -----------------------  */
main(argc, argv)
int argc;
Char *argv[];
{
  _TEXT TEMP;

  /* BEGIN MODULE STARTUP */
  /* Peform operations which must be done at beginning of main
     procedure. */
  /*!!!   TERMIN(input);    Open input for interactive use */
  /*!!!   TERMOUT(output);   "   output "      "        "  */
  PASCAL_MAIN(argc, argv);
  OUTFILE.f = NULL;
  *OUTFILE.name = '\0';
  INFILE.f = NULL;
  *INFILE.name = '\0';
  printf("%50s\n\n", VERSION);
  /* END MODULE STARTUP         VERSION= 'SUNMODS     Version  6/26/01'; */
  INITIALIZE();
  /* Initialize horizontal output line */

  /* Read in initial sequence and set up arrays. */
  for (I = 1; I <= MAXLINE; I++)
    HLINE.STR[I-1] = '_';
  HLINE.LEN = MAXLINE;
  printf("Enter sequence filename:\n");
  GETFILE(&INFILE, 'I', &IFN);
  NAME.LEN = 0;
  READSEQ(&INFILE, SEQ, &SEQLEN, &NAME, &CIRCULAR);
  if (NAME.LEN == 0) {
    printf("Type name to appear on output:\n");
    INPWORD(&NAME);
    scanf("%*[^\n]");
    getchar();
  }
  INITENDS();

  OFN.LEN = 0;   /* indicates that OUTFILE is not yet open */

  /* MAIN MENU */
  do {
    putchar('\n');
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, HLINE, 80L);
    printf("\nINTREST%30s\n", "MAIN MENU");
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, HLINE, 80L);
    printf("\nInput file:        ");
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, IFN, IFN.LEN);
    printf("\nOutput file:       ");
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, OFN, OFN.LEN);
    putchar('\n');
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, HLINE, 80L);
    printf("\n%20c1) Read in a new sequence\n", ' ');
    printf("%20c2) Open a new output file\n", ' ');
    printf("%20c3) Search for sites (output to screen)\n", ' ');
    printf("%20c4) Search for sites (output to file)\n", ' ');
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, HLINE, 80L);
    printf("\nType the number of your choice  (0 to quit program)\n");
    GETINTEGER(&CHOICE, 0L, 4L);
    scanf("%*[^\n]");
    getchar();

    switch (CHOICE) {

    case 0:
      /* blank case */
      break;

    case 1:
      printf("Enter sequence filename:\n");
      GETFILE(&INFILE, 'I', &IFN);
      NAME.LEN = 0;
      READSEQ(&INFILE, SEQ, &SEQLEN, &NAME, &CIRCULAR);
      if (NAME.LEN == 0) {
	printf("Type name to appear on output:\n");
	INPWORD(&NAME);
	scanf("%*[^\n]");
	getchar();
      }
      INITENDS();
      break;

    case 2:
      printf("Type output filename:\n");
      GETFILE(&OUTFILE, 'O', &OFN);
      break;

    /*!!!         CLOSE(OUTFILE); */
    case 3:
      TEMP.f = stdout;
      *TEMP.name = '\0';
      PROMPT(&TEMP);
      break;

    case 4:
      if (OFN.LEN > 0)
	PROMPT(&OUTFILE);
      else {
	printf(">>> NO OUTPUT FILE CURRENTLY OPEN\n");
	scanf("%*[^\n]");
	getchar();
      }
      break;
    }
  } while (CHOICE != 0);
  /*!!!*/
  if (OUTFILE.f != NULL)
    fclose(OUTFILE.f);
  OUTFILE.f = NULL;
  if (INFILE.f != NULL)
    fclose(INFILE.f);
  if (OUTFILE.f != NULL)
    fclose(OUTFILE.f);
  exit(EXIT_SUCCESS);
}  /* INTREST */




/* End. */
