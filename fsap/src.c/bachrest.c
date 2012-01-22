/* Output from p2c 1.21alpha-07.Dec.93, the Pascal-to-C translator */
/* From input file "bachrest.p" */


/**********************************************************/
/*                                                        */
/*    BACHREST Version  3/29/2006,  Standard Pascal       */
/*             Brian Fristensky                           */
/*             Dept. of Plant Science                     */
/*             University of Manitoba                     */
/*             Winnipeg, MB R3T 2N2  CANADA               */
/*                                                        */
/* Copyright (c) 1986, 1987, 1990  by Brian Fristensky    */
/* !!! in comment indicates feature which may need change */
/******************************************************** */
/* REVISION HISTORY
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
*/

#include <p2c.h>


/*,INFILE,RESTFILE,OUTFILE*/
/*!!! Other Pascals may require file parameters in heading */

#define MAXSEQ          20000000L

/* BEGIN MODULE INTLIMITS */
/*!!!  MAXINT =  2147483647; */
/*!!!  MININT = -2147483647; */
/* END MODULE INTLIMITS         VERSION= 'SUNMODS     Version  6/26/01'; */

#define MAXPAT          23
#define MAXFRAGS        6000
#define MAXWORD         23
#define MAXLINE         150

#define VERSION         "BACHREST   Version  3/29/2006"


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
  WORD NAME, PROTONAME, RECSTR;   /* string value of recog. seq.*/
  long PLEN;
  NUCLEOTIDE RSEQ[MAXPAT];
  boolean SYMMETRIC;
  long CUT, CUTOPP;
  Char FRAGENDS;
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


Static _TEXT INFILE, RESTFILE, OUTFILE;
Static LINE IFN, RFN, OFN;   /* File names */
Static boolean REBASE;
    /*=true if RESTFILE is REBASE, false if FSAP format */
Static WORD NAME;
Static NUCLEOTIDE COM[17];   /*Complements of NUCLEOTIDE */
Static SS NUCSET[17];
Static long LEGALNUCS[9];
Static SS AMBIGUOUS;
Static NUCLEOTIDE ENDNUC;
Static ENZYME ENZ;   /*Enzyme site to search for */
Static FRAGSFOUND FOUND;   /*List of sites found*/
Static FRAGMENT *FREEFRAG;   /*Points to freelist of sites*/
Static SEQUENCE SEQ;   /*Sequence to be searched*/
Static long SEQLEN;   /*Sequence length */
Static boolean CIRCULAR;   /*=true if seq. is circular */
Static LINE TEMPLINE;   /* dummy input line */
Static LINE HLINE;   /* horizontal line for menus */
Static long I, J, FIRSTITEM;   /* line number of first REBASE entry */
Static long CHOICE;   /* # of menu choice */

/* Global search parameters */
Static Char SOURCE;   /* Search for commercial or all */
Static Char PROTOTYPE;   /* Search for prototype, or all isoscizomers */
Static Char PROT3;   /* Search for 3' protruding end sites */
Static Char BLUNT;   /* Search for blunt end sites */
Static Char PROT5;   /* Search for 5' protruding end sites */
Static Char SYMM;   /* Search for symmetric, assym. or both */
Static long MINSITE;   /* Minimum length of RE seq. to search for */
Static long MAXSITE;   /* Maximum length of RE seq. to search for */
Static long FRAGLEAST;   /* Minimum number of fragments in a digest */
Static long FRAGMOST;   /* Maximum number of fragments in a digest */
Static long FRAGPRINT;   /* Maximum number of fragments to print */
Static boolean FITSPARAMS;   /* Last enzyme read fits parameter settings */


/* for type of RE to search for*/


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


/* BEGIN MODULE SKIPBLANKS */
/* Skip to next non-blank character on the current line. */
Static Void SKIPBLANKS(F)
_TEXT *F;
{
  boolean DONE = false;
  Char CH;

  do {
    if (P_eoln(F->f) | BUFEOF(F->f))
      DONE = true;
    else if (P_peek(F->f) == ' ') {
      CH = getc(F->f);
      if (CH == '\n')
	CH = ' ';
    } else
      DONE = true;
  } while (!DONE);   /* SKIPBLANKS */
}


/* END MODULE SKIPBLANKS */


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


/* END MODULE READLINE         VERSION= 'SUNMODS     Version  6/26/01'; */

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
/* p2c: bachrest.p, line 563: 
 * Note: Line breaker spent 0.0 seconds, 5000 tries on line 813 [251] */
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

/************************************************/
/*  Initialize arrays.                          */
/************************************************/
Static Void INITIALIZE()
{
  NUCLEOTIDE B1 = A, B2 = T;

  P_addset(P_expset(LEGALNUCS, 0L), 'A');
  P_addset(LEGALNUCS, 'a');
  P_addset(LEGALNUCS, 'C');
  P_addset(LEGALNUCS, 'c');
  P_addset(LEGALNUCS, 'R');
  P_addset(LEGALNUCS, 'r');
  P_addset(LEGALNUCS, 'D');
  P_addset(LEGALNUCS, 'd');
  P_addset(LEGALNUCS, 'V');
  P_addset(LEGALNUCS, 'v');
  P_addset(LEGALNUCS, 'M');
  P_addset(LEGALNUCS, 'm');
  P_addset(LEGALNUCS, 'K');
  P_addset(LEGALNUCS, 'k');
  P_addset(LEGALNUCS, 'B');
  P_addset(LEGALNUCS, 'b');
  P_addset(LEGALNUCS, 'H');
  P_addset(LEGALNUCS, 'h');
  P_addset(LEGALNUCS, 'Y');
  P_addset(LEGALNUCS, 'y');
  P_addset(LEGALNUCS, 'G');
  P_addset(LEGALNUCS, 'g');
  P_addset(LEGALNUCS, 'T');
  P_addset(LEGALNUCS, 't');
  P_addset(LEGALNUCS, 'N');
  P_addset(LEGALNUCS, 'n');
  P_addset(LEGALNUCS, 'S');
  P_addset(LEGALNUCS, 's');
  P_addset(LEGALNUCS, 'W');
  P_addset(LEGALNUCS, 'w');
  /* Legal nucleotide characters */
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

  /* Initialize search parameters */
  SOURCE = 'C';
  PROTOTYPE = 'A';
  PROT3 = 'Y';
  BLUNT = 'Y';
  PROT5 = 'Y';
  SYMM = 'B';
  MINSITE = 4;
  MAXSITE = MAXPAT;
  FRAGLEAST = 0;
  FRAGMOST = MAXFRAGS;
  FRAGPRINT = 30;

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


#define HLINE_          "-----------------------------------------------------------"


/*******************************************************************/
/* Print a header to OUTFILE.                                      */
/*******************************************************************/
Static Void HEADER(OUTFILE)
_TEXT *OUTFILE;
{
  fprintf(OUTFILE->f, "%s\n", HLINE_);
  fprintf(OUTFILE->f, "%50s\n\n", VERSION);
  WRITEWORD(OUTFILE, NAME, 18L);
  fprintf(OUTFILE->f, "  Topology: ");
  switch (CIRCULAR) {

  case true:
    fprintf(OUTFILE->f, "CIRCULAR");
    break;

  case false:
    fprintf(OUTFILE->f, "LINEAR");
    break;
  }
  fprintf(OUTFILE->f, "  Length: %8ld bp\n", SEQLEN);
  fprintf(OUTFILE->f, "%s\n", HLINE_);
  fprintf(OUTFILE->f, "Search parameters:\n");
  fprintf(OUTFILE->f, "   Recognition sequences between %4ld and %4ld bp\n",
	  MINSITE, MAXSITE);

  fprintf(OUTFILE->f, "   Ends: ");
  if (PROT5 == 'Y')
    fprintf(OUTFILE->f, "5' protruding");
  if (PROT5 == 'Y' && (BLUNT == 'Y' || PROT3 == 'Y'))
    fprintf(OUTFILE->f, ", ");
  if (BLUNT == 'Y')
    fprintf(OUTFILE->f, "Blunt");
  if (BLUNT == 'Y' && PROT3 == 'Y')
    fprintf(OUTFILE->f, ", ");
  if (PROT3 == 'Y')
    fprintf(OUTFILE->f, "3' protruding");
  fprintf(OUTFILE->f, "\n   Type: ");

  if (SYMM == 'S')
    fprintf(OUTFILE->f, "Symmetric\n");
  else if (SYMM == 'B')
    fprintf(OUTFILE->f, "Symmetric, Asymmetric\n");
  else
    fprintf(OUTFILE->f, "Asymmetric\n");

  fprintf(OUTFILE->f, "   Minimum fragments: %5ld     ", FRAGLEAST);
  fprintf(OUTFILE->f, "Maximum fragments: %5ld\n", FRAGMOST);
  fprintf(OUTFILE->f, "   Maximum fragments to print: %5ld\n", FRAGPRINT);

  fprintf(OUTFILE->f, "%s\n\n", HLINE_);
  fprintf(OUTFILE->f, "%*s\n", (int)(MAXPAT + 22), "# of");
  fprintf(OUTFILE->f,
    "Enzyme          Recognition Sequence     Sites     Sites   Frags   Begin     End\n");
  fprintf(OUTFILE->f, "%s---------------------\n\n", HLINE_);
}  /* HEADER */

#undef HLINE_


/**************************************************************/
/* Determine whether RESTFILE is in REBASE or FSAP format.    */
/* if REBASE = true, then RESTFILE^ points to the first entry.*/
/* if REBASE = false, the end of file will be reached.        */
/**************************************************************/
Static Void CHECKFORMAT(RESTFILE, REBASE)
_TEXT *RESTFILE;
boolean *REBASE;
{
  long I = 1;

  while ((P_peek(RESTFILE->f) != '<') & (!BUFEOF(RESTFILE->f))) {
    fscanf(RESTFILE->f, "%*[^\n]");
    getc(RESTFILE->f);
    I++;
  }
  if (BUFEOF(RESTFILE->f))
    *REBASE = false;   /* FSAP format */
  else {
    *REBASE = true;
    /* REBASE */

  }
}  /* CHECKFORMAT */


/******************************************************/
/*  Read in a restriction enzyme site from RESTFILE.  */
/******************************************************/
Static Void READSITE(RESTFILE, ENZ)
_TEXT *RESTFILE;
ENZYME *ENZ;
{
  long LOC;
  boolean LEGAL;
  long FORLIM;

  do {
    LEGAL = true;
    /* Read enzyme name, recognition site, and cutting site from  file */
    READWORD(RESTFILE, &ENZ->NAME);
    READWORD(RESTFILE, &ENZ->RECSTR);
    fscanf(RESTFILE->f, "%ld", &ENZ->CUT);
    SKIPBLANKS(RESTFILE);
    if (P_eoln(RESTFILE->f)) {
      fscanf(RESTFILE->f, "%*[^\n]");
      getc(RESTFILE->f);
    } else if (!BUFEOF(RESTFILE->f)) {
      fscanf(RESTFILE->f, "%ld", &ENZ->CUTOPP);
      if (P_eoln(RESTFILE->f)) {
	fscanf(RESTFILE->f, "%*[^\n]");
	getc(RESTFILE->f);
      }
    }
    FORLIM = ENZ->RECSTR.LEN;

    /* Check for illegal symbols */
    for (LOC = 0; LOC < FORLIM; LOC++) {
      if (P_inset(ENZ->RECSTR.STR[LOC], LEGALNUCS))
	ENZ->RSEQ[LOC] = NUC(ENZ->RECSTR.STR[LOC]);
      else
	LEGAL = false;
    }
    if (!LEGAL)
      fprintf(OUTFILE.f, ">>>> Illegal input ignored\n");
  } while (!LEGAL);
  ENZ->PLEN = ENZ->RECSTR.LEN;


}  /* READSITE */


/* Local variables for READREBASE: */
struct LOC_READREBASE {
  long ORDZERO;
} ;

/* Extract an integer from a WORD */
Local Void GETNUM(W_, POSN, NUM, LINK)
WORD W_;
long *POSN, *NUM;
struct LOC_READREBASE *LINK;
{
  long SIGN = 1;

  /* with W */
  *NUM = 0;
  if (W_.STR[*POSN - 1] == '-') {
    SIGN = -1;
    (*POSN)++;
  }
  while (isdigit(W_.STR[*POSN - 1]) && *POSN <= W_.LEN) {
    *NUM = *NUM * 10 + W_.STR[*POSN - 1] - LINK->ORDZERO;
    (*POSN)++;
  }
  if (SIGN == -1)
    *NUM = -*NUM;
}  /* GETNUM */


/******************************************************/
/*  Read in a restriction enzyme site from RESTFILE.  */
/******************************************************/
Static Void READREBASE(RESTFILE, ENZ, FITSPARAMS)
_TEXT *RESTFILE;
ENZYME *ENZ;
boolean *FITSPARAMS;
{
  struct LOC_READREBASE V;
  WORD COMSOURCE;
  boolean CARET = true;   /* G^AATTC format, rather than GAAGA(8/7) */
  long POSN = 1;
  long I;
  Char CH;

  /* with ENZ */
  V.ORDZERO = '0';

  /* Name begins at 4th char. or 1st line of a record */
  for (I = 1; I <= 3; I++) {
    CH = getc(RESTFILE->f);
    if (CH == '\n')
      CH = ' ';
  }
  READWORD(RESTFILE, &ENZ->NAME);
  fscanf(RESTFILE->f, "%*[^\n]");
  getc(RESTFILE->f);

  /* Prototype begins at 4th char. on 2nd line of a record */
  for (I = 1; I <= 3; I++) {
    CH = getc(RESTFILE->f);
    if (CH == '\n')
      CH = ' ';
  }
  if (!P_eoln(RESTFILE->f))
    READWORD(RESTFILE, &ENZ->PROTONAME);
  else
    ENZ->PROTONAME.LEN = 0;
  fscanf(RESTFILE->f, "%*[^\n]");
  getc(RESTFILE->f);
  if (PROTOTYPE == 'P' && ENZ->PROTONAME.LEN > 0)
    *FITSPARAMS = false;

  /* Site begins at 4th char. of line 3 */
  for (I = 1; I <= 3; I++) {
    CH = getc(RESTFILE->f);
    if (CH == '\n')
      CH = ' ';
  }
  READWORD(RESTFILE, &ENZ->RECSTR);
  fscanf(RESTFILE->f, "%*[^\n]");
  getc(RESTFILE->f);

  /* Extract cut sites from RECSTR, and write the site to RSEQ */
  ENZ->PLEN = 0;

  ENZ->CUT = 0;   /* if no cut site specified, assume zero */
  while (POSN <= ENZ->RECSTR.LEN) {
    if (P_inset(ENZ->RECSTR.STR[POSN-1], LEGALNUCS)) {  /* recognition seq. */
      ENZ->PLEN++;
      ENZ->RSEQ[ENZ->PLEN - 1] = NUC(ENZ->RECSTR.STR[POSN-1]);
      POSN++;
      continue;
    }
    if (ENZ->RECSTR.STR[POSN-1] == '^') {  /* symmetric cut */
      ENZ->CUT = POSN - 1;
      ENZ->CUTOPP = ENZ->CUT;
      POSN++;
      continue;
    }
    if (ENZ->RECSTR.STR[POSN-1] != '(') {  /* asymmetric cut */
      POSN++;
      continue;
    }
    CARET = false;
    POSN++;
    GETNUM(ENZ->RECSTR, &POSN, &ENZ->CUT, &V);
    POSN++;
    GETNUM(ENZ->RECSTR, &POSN, &ENZ->CUTOPP, &V);
  }

  /* To keep the specification of cut sites consistent, cut sites
  specified using caret (eg. G^AATTC) are converted to coordinates
  consistent with parenthesis notation (eg. GAAGA(8/7)). In caret
  notation, CUT is with reference to the position 5' to the start
  of the recognition sequence. In parenthesis notation, CUT is
  with reference to the 3' end of the sequence. */
  if (CARET) {
    ENZ->CUT -= ENZ->PLEN;
    ENZ->CUTOPP = ENZ->CUT;
  }


  fscanf(RESTFILE->f, "%*[^\n]");
  getc(RESTFILE->f);   /* skip methylation data */

  /* Read list of commercial suppliers. */
  for (I = 1; I <= 3; I++) {
    CH = getc(RESTFILE->f);
    if (CH == '\n')
      CH = ' ';
  }
  if (!P_eoln(RESTFILE->f))
    READWORD(RESTFILE, &COMSOURCE);
  else
    COMSOURCE.LEN = 0;
  fscanf(RESTFILE->f, "%*[^\n]");
  getc(RESTFILE->f);
  if (SOURCE == 'C' && COMSOURCE.LEN == 0)
    *FITSPARAMS = false;
  fscanf(RESTFILE->f, "%*[^\n]");
  getc(RESTFILE->f);   /* skip reference data */
  while ((P_peek(RESTFILE->f) != '<') & (!BUFEOF(RESTFILE->f))) {
    fscanf(RESTFILE->f, "%*[^\n]");
    getc(RESTFILE->f);
  }

}  /* READREBASE */


typedef Char LETTERS[10];
typedef long CHSET[4];


/*  Read an integer parameter from the console and check */
/*    that it is in range.                               */
Local Void GETNUMBER(P, PNAME, LOW, HIGH)
long *P;
Char *PNAME;
long LOW, HIGH;
{
  printf("Type new value for %.10s  (CURRENT VALUE: %12ld)\n", PNAME, *P);
  GETINTEGER(P, LOW, HIGH);
  scanf("%*[^\n]");
  getchar();
}  /* GETNUMBER */

Local Char TOUPPER(CH)
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

/* BEGIN MODULE GETCHAR */
/* Read a character from the console and check */
/*  for correct response.                      */
Local Void GETCHAR(CH, PNAME, ALLOWED)
Char *CH;
Char *PNAME;
long *ALLOWED;
{
  /* BEGIN MODULE TOUPPER */
  /* Change a character from lower to uppercase */
  /* END MODULE TOUPPER         VERSION= 'SUNMODS     Version  6/26/01'; */
  putchar('\n');
  do {
    printf("Type new value for %.10s  (CURRENT VALUE: %c)\n", PNAME, *CH);
    scanf("%c%*[^\n]", CH);
    getchar();
    if (*CH == '\n')
      *CH = ' ';
    *CH = TOUPPER(*CH);
    if (!P_inset(*CH, ALLOWED))
      printf("Inappropriate response: %c\n", *CH);
  } while (!P_inset(*CH, ALLOWED));   /* GETCHAR */
}

/* END MODULE GETCHAR         VERSION= 'SUNMODS     Version  6/26/01'; */

/* Display  parameters on screen */
Local Void DISPLAY()
{
  _TEXT TEMP;

  TEMP.f = stdout;
  *TEMP.name = '\0';
  WRITELINE(&TEMP, HLINE, 80L);
  printf("\n%12cParameter   Description/Response                     Value\n",
	 ' ');
  TEMP.f = stdout;
  *TEMP.name = '\0';
  WRITELINE(&TEMP, HLINE, 80L);
  printf("\n%12c 1)SOURCE    C:Commercial only  A:all%20c\n", ' ', SOURCE);
  printf("%12c 2)PROTOTYPE P:prototypes only  A:all%20c\n", ' ', PROTOTYPE);
  printf("%12c 3)PROT3     3' protruding end cutters (Y/N)%13c\n", ' ', PROT3);
  printf("%12c 4)BLUNT     Blunt end cutters         (Y/N)%13c\n", ' ', BLUNT);
  printf("%12c 5)PROT5     5' protruding end cutters (Y/N)%13c\n", ' ', PROT5);
  printf("%12c 6)SYMM      S:symmetric A:asymmetric B:both%13c\n", ' ', SYMM);
  printf("%12c 7)MINSITE   Minimum RE site length           %11ld\n",
	 ' ', MINSITE);
  printf("%12c 8)MAXSITE   Maximum RE site length           %11ld\n",
	 ' ', MAXSITE);
  printf("%12c 9)FRAGLEAST Min. # of fragments to print a digest%7ld\n",
	 ' ', FRAGLEAST);
  printf("%12c10)FRAGMOST  Max. # of fragments to print a digest%7ld\n",
	 ' ', FRAGMOST);
  printf("%12c11)FRAGPRINT Print a number if > FRAGPRINT frags%9ld\n",
	 ' ', FRAGPRINT);
  TEMP.f = stdout;
  *TEMP.name = '\0';
  WRITELINE(&TEMP, HLINE, 80L);
  putchar('\n');
}  /* DISPLAY */

Local Void REBASEMSG()
{
  printf(">>> This parameter has no effect unless REBASE is used\n");
  printf(">>> as the restriction site file. Press ENTER to continue\n");
  scanf("%*[^\n]");
  getchar();
}  /* REBASEMSG */


/* **************************************************** */
/* Prompt user for parameters used by program.          */
/* **************************************************** */
Static Void PARAMETERS()
{
  long RESPONSE;
  long SET[4];
  long SET1[4];
  long SET2[4];
  long SET3[4];

  /* Prompt user for new parameter values */
  do {
    fprintf(stdout, "\f");
    DISPLAY();
    printf("Type number of parameter you wish to change (0 to continue)\n");
    GETINTEGER(&RESPONSE, 0L, 11L);
    scanf("%*[^\n]");
    getchar();
    if ((unsigned long)RESPONSE < 32 && ((1L << RESPONSE) & 0xffe) != 0) {
      switch (RESPONSE) {

      case 1:
	if (REBASE) {
	  P_addset(P_expset(SET, 0L), 'C');
	  GETCHAR(&SOURCE, "SOURCE    ", P_addset(SET, 'A'));
	} else
	  REBASEMSG();
	break;

      case 2:
	if (REBASE) {
	  P_addset(P_expset(SET1, 0L), 'P');
	  GETCHAR(&PROTOTYPE, "PROTOTYPE ", P_addset(SET1, 'A'));
	} else
	  REBASEMSG();
	break;

      case 3:
	P_addset(P_expset(SET2, 0L), 'Y');
	GETCHAR(&PROT3, "PROT3     ", P_addset(SET2, 'N'));
	break;

      case 4:
	P_addset(P_expset(SET2, 0L), 'Y');
	GETCHAR(&BLUNT, "BLUNT     ", P_addset(SET2, 'N'));
	break;

      case 5:
	P_addset(P_expset(SET2, 0L), 'Y');
	GETCHAR(&PROT5, "PROT5     ", P_addset(SET2, 'N'));
	break;

      case 6:
	P_addset(P_expset(SET3, 0L), 'S');
	P_addset(SET3, 'A');
	GETCHAR(&SYMM, "SYMM      ", P_addset(SET3, 'B'));
	break;

      case 7:
	GETNUMBER(&MINSITE, "MINSITE   ", 4L, MAXSITE);
	break;

      case 8:
	GETNUMBER(&MAXSITE, "MAXSITE   ", MINSITE, (long)MAXPAT);
	break;

      case 9:
/* p2c: bachrest.p, line 920:
 * Warning: Too many characters for packed array of char [162] */
	GETNUMBER(&FRAGLEAST, "FRAGLEAST   ", 0L, FRAGMOST);
	break;

      case 10:
/* p2c: bachrest.p, line 921:
 * Warning: Too many characters for packed array of char [162] */
	GETNUMBER(&FRAGMOST, "FRAGMOST   ", FRAGLEAST, (long)MAXFRAGS);
	break;

      case 11:
/* p2c: bachrest.p, line 922:
 * Warning: Too many characters for packed array of char [162] */
	GETNUMBER(&FRAGPRINT, "FRAGPRINT  ", 1L, FRAGMOST);
	break;
      }
    }
  } while (RESPONSE != 0);   /* PARAMETERS */
}


/**************************************************************/
/* Determine whether the recognition sequence is symmetric .  */
/**************************************************************/
Static Void TESTSYMMETRY(ENZ)
ENZYME *ENZ;
{
  long Cj, Pj;
  NUCLEOTIDE CRP;
  long FORLIM;

  ENZ->SYMMETRIC = true;
  Cj = ENZ->PLEN;
  FORLIM = ENZ->PLEN;
  for (Pj = 0; Pj < FORLIM; Pj++) {
    CRP = COM[(long)ENZ->RSEQ[Pj]];
    if (ENZ->RSEQ[Cj-1] != CRP)
      ENZ->SYMMETRIC = false;
    Cj--;
  }
}  /* TESTSYMMETRY */


/**************************************************************/
/* Determine the type of cutting site.                        */
/**************************************************************/
Static Void TESTCUT(ENZ)
ENZYME *ENZ;
{
  long CUTSAFTER;

  if (ENZ->SYMMETRIC) {
    if (REBASE)
      CUTSAFTER = ENZ->PLEN + ENZ->CUT;
    else
      CUTSAFTER = ENZ->CUT;
    if (CUTSAFTER < ENZ->PLEN / 2.0) {
      ENZ->FRAGENDS = '5';
      return;
    }
    if (CUTSAFTER > ENZ->PLEN / 2.0)
      ENZ->FRAGENDS = '3';
    else
      ENZ->FRAGENDS = 'B';
    return;
  }
  if (ENZ->CUT < ENZ->CUTOPP) {
    ENZ->FRAGENDS = '5';
    return;
  }
  if (ENZ->CUT > ENZ->CUTOPP)
    ENZ->FRAGENDS = '3';
  else
    ENZ->FRAGENDS = 'B';
}  /* TESTCUT */


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
}  /* MAKENEXT */

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
  long POSN;   /* # of nucleotide 5' to the cut on input strand */
  ENZYME *WITH;

  WITH = LINK->LINK->ENZ;
  POSN = LINK->Sk - WITH->PLEN + LINK->CUTSITE +
	 LINK->LINK->GOBACK[LINK->LINK->PLEN1-1];
  if (CIRCULAR) {   /* Test for enz. which cuts beyond recog. site*/
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
  _TEXT TEMP;
  ENZYME *WITH;

  V.LINK = LINK;
  V.CUTSITE = CUTSITE_;
  TEMP.f = stdout;
  *TEMP.name = '\0';
  WRITEWORD(&TEMP, LINK->ENZ->NAME, 10L);
  putchar('\n');
  WITH = LINK->ENZ;
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

_L99:   /* INPUTEXHAUSTED */
  putchar('\n');
}  /* SEARCH */


/* END MODULE LINKED         VERSION= 'SUNMODS     Version  6/26/01'; */


/****************************************************************/
/*  Search for site ENZ in SEQ, using the algorithm found in    */
/*    "Fast Pattern Matching in Strings"                        */
/*     Knuth, Morris, and Pratt, SIAM Journal of Computing,     */
/*     Vol.6, No.2, June 1977.  pp323-350.                      */
/*                                                              */
/*     The algorithm works in O(m+n) units of time, where n is  */
/*     the length of the sequence searched, and m is the length */
/*     of the pattern to search for.                            */
/****************************************************************/
Static Void KMP(ENZ_)
ENZYME *ENZ_;
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
  V.PATTERN[WITH->PLEN] = 1L << ((long)CANT);   /* Can not match anything */
  V.RIGHTMOST = FOUND.HEAD;
  MAKENEXT(&V);
  if (REBASE)
    SEARCH(WITH->CUT + WITH->PLEN, &V);
  else
    SEARCH(WITH->CUT, &V);

  /* Let PATTERN = the inverse complement of PATTERN*/
  /* Check for symmetry.                            */
  Cj = WITH->PLEN;
  FORLIM = WITH->PLEN;
  for (Pj = 0; Pj < FORLIM; Pj++) {
    CRP = COM[(long)WITH->RSEQ[Pj]];
    V.PATTERN[Cj-1] = NUCSET[(long)CRP];
    if (WITH->RSEQ[Cj-1] != CRP)
      WITH->SYMMETRIC = false;
    Cj--;
  }

  if (WITH->SYMMETRIC)
    return;
  ENDNUC = CRP;
  MAKENEXT(&V);
  V.RIGHTMOST = FOUND.HEAD;
  if (REBASE)
    SEARCH(-WITH->CUTOPP, &V);
  else
    SEARCH(WITH->PLEN - WITH->CUTOPP, &V);
}  /* KMP */


typedef FRAGMENT *SORTEDLIST[MAXFRAGS];


/* Local variables for REPORT: */
struct LOC_REPORT {
  SORTEDLIST ORDER;
  boolean TOOMANY_FRAGS;   /* =true if SITES > MAXFRAGS */
  long NUMSITES;
} ;

/*Calculate sizes and ends of fragments*/
Local Void CALCULATE(ORDER, LINK)
FRAGMENT **ORDER;
struct LOC_REPORT *LINK;
{
  FRAGMENT *CURRENTFRAG;
  long SITE = 1;
  long I;

  /* Initialize order array, to prevent optimizer errors */
  for (I = 0; I < MAXFRAGS; I++)
    ORDER[I] = NULL;

  /* with FOUND */
  LINK->NUMSITES = FOUND.LNUM;
  /* if LINEAR add a fragment to head of list*/
  if (!CIRCULAR && FOUND.HEAD->NEXT->START > 1) {
    GETFRAG(&CURRENTFRAG);
    ADDFRAG(&CURRENTFRAG, &FOUND.HEAD);
    CURRENTFRAG->START = 1;
    FOUND.LNUM++;
  } else
    CURRENTFRAG = FOUND.HEAD->NEXT;

  if (FOUND.LNUM <= 0) {
    return;
  }  /* if LNUM > 0 */
  /*Calculate ends and size of each fragment and assign it */
  /* a place in the ORDER array, to be sorted later.    */
  while (CURRENTFRAG->NEXT != FOUND.TAIL) {
    CURRENTFRAG->FINISH = CURRENTFRAG->NEXT->START - 1;
    CURRENTFRAG->SIZE = CURRENTFRAG->FINISH - CURRENTFRAG->START + 1;
    if (SITE <= MAXFRAGS)
      ORDER[SITE-1] = CURRENTFRAG;
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
    ORDER[SITE-1] = CURRENTFRAG;
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
Local Void PRINTLIST(ORDER, LINK)
FRAGMENT **ORDER;
struct LOC_REPORT *LINK;
{
  FRAGMENT *THISFRAG;
  long SITE;
  FRAGMENT *WITH;

  fprintf(OUTFILE.f, "%6ld\n", LINK->NUMSITES);
  if (LINK->NUMSITES <= FRAGPRINT) {
    if (LINK->NUMSITES == FOUND.LNUM)
      THISFRAG = FOUND.HEAD;
    else
      THISFRAG = FOUND.HEAD->NEXT;
    SITE = 1;
    while (THISFRAG->NEXT != FOUND.TAIL) {
      THISFRAG = THISFRAG->NEXT;
      fprintf(OUTFILE.f, "%*ld", (int)(MAXPAT + 33), THISFRAG->START);
      WITH = ORDER[SITE-1];
      fprintf(OUTFILE.f, "%8ld%8ld%8ld\n",
	      WITH->SIZE, WITH->START, WITH->FINISH);
      SITE++;
    }
    /* Print last fragment. If no sites are found, a linear */
    /* sequence will still yeild one fragment.              */
    if (LINK->NUMSITES < FOUND.LNUM) {
      WITH = ORDER[SITE-1];
      fprintf(OUTFILE.f, "%*ld%8ld%8ld\n",
	      (int)(MAXPAT + 41), WITH->SIZE, WITH->START, WITH->FINISH);
    }
  } else {
    fprintf(OUTFILE.f, "%*c(Fragments not shown)\n", (int)(MAXPAT + 33), ' ');

    /* This now looks unnecessary, but we'll keep it as a comment.
    if TOOMANY_FRAGS then begin
       writeln('>>> Number of fragments exceeds ',MAXFRAGS);
       writeln('>>> Increase MAXFRAGS and recompile');
       writeln(OUTFILE);
       writeln(OUTFILE, '>>> Number of fragments exceeds ',MAXFRAGS);
       writeln(OUTFILE, '>>> Increase MAXFRAGS and recompile');
       end
       */
  }
  putc('\n', OUTFILE.f);
}  /* PRINTLIST */


/*********************************************/
/* Compile a report for output.              */
/*********************************************/
/* Note (3/4/93). The original BACHREST had ORDER as a global array
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
   Regardless of the cause, it works. */

Static Void REPORT()
{
  struct LOC_REPORT V;

  V.TOOMANY_FRAGS = false;
  CALCULATE(V.ORDER, &V);
  if (FOUND.LNUM >= FRAGLEAST &&
      (FOUND.LNUM <= FRAGMOST || FRAGMOST == MAXFRAGS)) {
    /* Print a header for each the enzyme */
    WRITEWORD(&OUTFILE, ENZ.NAME, 15L);
    putc(' ', OUTFILE.f);
    WRITEWORD(&OUTFILE, ENZ.RECSTR, (long)MAXPAT);
    if (REBASE)
      putc(' ', OUTFILE.f);
    else {
      fprintf(OUTFILE.f, "%3ld", ENZ.CUT);
      if (ENZ.SYMMETRIC)
	fprintf(OUTFILE.f, "%6c", ' ');
      else
	fprintf(OUTFILE.f, " (%3ld)", ENZ.CUTOPP);
    }
    if (!V.TOOMANY_FRAGS)
      BUBBLESORT(FOUND.LNUM, 1L, &V);
    PRINTLIST(V.ORDER, &V);
  }
  RIDOF(&FOUND.HEAD, &FOUND.TAIL);
}  /* REPORT */



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
  RESTFILE.f = NULL;
  *RESTFILE.name = '\0';
  INFILE.f = NULL;
  *INFILE.name = '\0';
  printf("%50s\n\n", VERSION);
  /* END MODULE STARTUP         VERSION= 'SUNMODS     Version  6/26/01'; */

  /* Initialize horizontal output line */

  for (I = 1; I <= MAXLINE; I++)
    HLINE.STR[I-1] = '_';
  HLINE.LEN = MAXLINE;
  printf("Enter filename of sequence to search:\n");
  GETFILE(&INFILE, 'I', &IFN);
  NAME.LEN = 0;
  READSEQ(&INFILE, SEQ, &SEQLEN, &NAME, &CIRCULAR);
  if (NAME.LEN == 0) {
    printf("Type name to appear on output:\n");
    INPWORD(&NAME);
    scanf("%*[^\n]");
    getchar();
  }
  INITIALIZE();
  INITENDS();

  /* Open the restriction site file, skip the first two title lines, */
  /*  and advance to first non-blank character.                      */
  printf("Enter restriction site filename:\n");
  GETFILE(&RESTFILE, 'I', &RFN);
  CHECKFORMAT(&RESTFILE, &REBASE);

  /* Open an output file. */
  printf("Enter output filename:\n");
  GETFILE(&OUTFILE, 'O', &OFN);

  /* MAIN MENU */
  do {
    putchar('\n');
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, HLINE, 80L);
    printf("\nBACHREST%30s\n", "MAIN MENU");
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, HLINE, 80L);
    printf("\nInput file:        ");
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, IFN, IFN.LEN);
    printf("\nRest. enz. file:   ");
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, RFN, RFN.LEN);
    printf("\nOutput file:       ");
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, OFN, OFN.LEN);
    putchar('\n');
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, HLINE, 80L);
    printf("\n%20c1) Read in a new sequence\n", ' ');
    printf("%20c2) Open a new rest. enz. file\n", ' ');
    printf("%20c3) Open a new output file\n", ' ');
    printf("%20c4) Set search parameters\n", ' ');
    printf("%20c5) Search for sites and write output to file\n", ' ');
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, HLINE, 80L);
    printf("\nType the number of your choice  (0 to quit program)\n");
    GETINTEGER(&CHOICE, 0L, 5L);
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
      printf("Type Rest. Enz. filename:\n");
      GETFILE(&RESTFILE, 'I', &RFN);
      CHECKFORMAT(&RESTFILE, &REBASE);
      break;

    /*!!!         CLOSE(RESTFILE); */
    case 3:
      printf("Type output filename:\n");
      GETFILE(&OUTFILE, 'O', &OFN);
      break;

    /*!!!         CLOSE(OUTFILE); */
    case 4:
      PARAMETERS();
      break;

    case 5:
      if (OFN.LEN > 0) {
	HEADER(&OUTFILE);

	/* Reset file pointer to beginning of RESTFILE */
	if (*RESTFILE.name != '\0') {
	  if (RESTFILE.f != NULL)
	    RESTFILE.f = freopen(RESTFILE.name, "r", RESTFILE.f);
	  else
	    RESTFILE.f = fopen(RESTFILE.name, "r");
	} else
	  rewind(RESTFILE.f);
	if (RESTFILE.f == NULL)
	  _EscIO2(FileNotFound, RESTFILE.name);
	RESETBUF(RESTFILE.f, Char);
	if (REBASE) {
	  J = 1;
	  while ((P_peek(RESTFILE.f) != '<') & (!BUFEOF(RESTFILE.f))) {
	    fscanf(RESTFILE.f, "%*[^\n]");
	    getc(RESTFILE.f);
	    J++;
	  }
	  if (!BUFEOF(RESTFILE.f)) {
	    /* The .lst files distributed with REBASE begin with header information
	       and an entry formula for each of six fields. If the header has
	       not been deleted, we need to read past it to the first true entry.*/
	    /* Set FIRSTITEM to # of line on which first data item begins */
	    do {
	      READLINE(&RESTFILE, &TEMPLINE);
	      J++;
	    } while (TEMPLINE.STR[1] != '1' || TEMPLINE.STR[3] == '<');
	    FIRSTITEM = J - 1;
	    /* Reset the file and advance to the first data item. */
	    if (*RESTFILE.name != '\0') {
	      if (RESTFILE.f != NULL)
		RESTFILE.f = freopen(RESTFILE.name, "r", RESTFILE.f);
	      else
		RESTFILE.f = fopen(RESTFILE.name, "r");
	    } else
	      rewind(RESTFILE.f);
	    if (RESTFILE.f == NULL)
	      _EscIO2(FileNotFound, RESTFILE.name);
	    RESETBUF(RESTFILE.f, Char);
	    J = 1;
	    while (J < FIRSTITEM) {
	      fscanf(RESTFILE.f, "%*[^\n]");
	      getc(RESTFILE.f);
	      J++;
	    }
	  }  /* if not eof */

	}  /* REBASE */
	else {
	  if (!BUFEOF(RESTFILE.f)) {
	    fscanf(RESTFILE.f, "%*[^\n]");
	    getc(RESTFILE.f);
	  }
	  if (!BUFEOF(RESTFILE.f)) {
	    fscanf(RESTFILE.f, "%*[^\n]");
	    getc(RESTFILE.f);
	  }
	}

	/* Search loop */
	while (!BUFEOF(RESTFILE.f)) {
	  FITSPARAMS = true;
	  if (REBASE)
	    READREBASE(&RESTFILE, &ENZ, &FITSPARAMS);
	  else
	    READSITE(&RESTFILE, &ENZ);

	  /* Check minimim and maximum lengths of cut sites
	  if ENZ.PLEN < MINSITE then FITSPARAMS:=false;
	  if ENZ.PLEN > MAXSITE then FITSPARAMS:=false;

	  (* Check for symmetry */
	  TESTSYMMETRY(&ENZ);
	  if (SYMM == 'S' && !ENZ.SYMMETRIC)
	    FITSPARAMS = false;
	  if (SYMM == 'A' && ENZ.SYMMETRIC)
	    FITSPARAMS = false;

	  /* Check for cutting position */
	  TESTCUT(&ENZ);
	  switch (ENZ.FRAGENDS) {

	  case '5':
	    if (PROT5 == 'N')
	      FITSPARAMS = false;
	    break;

	  case '3':
	    if (PROT3 == 'N')
	      FITSPARAMS = false;
	    break;

	  case 'B':
	    if (BLUNT == 'N')
	      FITSPARAMS = false;
	    break;
	  }

	  /* Run the search. */
	  if (FITSPARAMS) {
	    KMP(&ENZ);
	    REPORT();
	  }
	}  /* while not eof */
      }

      else {
	printf(">>> NO OUTPUT FILE CURRENTLY OPEN\n");
	scanf("%*[^\n]");
	getchar();
      }
      break;
    }
  } while (CHOICE != 0);

  /*!!!*/
  if (RESTFILE.f != NULL)
    fclose(RESTFILE.f);
  RESTFILE.f = NULL;
  /*!!!*/
  if (OUTFILE.f != NULL)
    fclose(OUTFILE.f);
  OUTFILE.f = NULL;
  if (INFILE.f != NULL)
    fclose(INFILE.f);
  if (RESTFILE.f != NULL)
    fclose(RESTFILE.f);
  if (OUTFILE.f != NULL)
    fclose(OUTFILE.f);
  exit(EXIT_SUCCESS);
}  /* BACHREST */




/* End. */
