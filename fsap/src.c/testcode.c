/* Output from p2c 1.21alpha-07.Dec.93, the Pascal-to-C translator */
/* From input file "testcode.p" */


/* ********************************************************  */
/*                                                           */
/*   TESTCODE  Version   8/13/2001 Standard Pascal           */
/*             Brian Fristensky                              */
/*             Dept. of Plant Science                        */
/*             University of Manitoba                        */
/*             Winnipeg, MB  R3T 2N2  CANADA                 */
/*                                                           */
/*  Copyright (c) 1984-1990  by Brian Fristensky             */
/*  !!! in comment indicates feature which may need change.  */
/*  *******************************************************  */
/* REVISION HISTORY
Aug. 13, 2001 Rebuilt using new Readseq procedure which reads
new GenBank format.

*/

/*!!!*/

#include <p2c.h>


/*,INFILE,OUTFILE*/
/*!!! Some Pascals require file parameters in program heading */

#define MAXSEQ          500000L

#define MINCODONS       30
#define MAXLINE         70
#define MAXWORD         25
#define LINESIZE        70

#define VERSION         "TESTCODE          Version   8/13/2001"


/* NUCLEOTIDE types */
typedef enum {
  T, C, A, G, N, R, Y, M, W, S, K, D, H, V, B, Z
} NUCLEOTIDE;
typedef NUCLEOTIDE NA[(long)N - (long)T + 1];
typedef double RA[(long)G - (long)T + 1];
typedef NUCLEOTIDE SEQUENCE[MAXSEQ];

/* table for probability values */

typedef struct PARAM {
  double LIMIT;
  double PROB[(long)G - (long)T + 1];
} PARAM;

typedef PARAM TABLE[10];

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


Static _TEXT INFILE, OUTFILE;
Static LINE IFN, OFN, HLINE;
Static SEQUENCE SEQ;
Static long SEQLEN;   /* sequence length */
Static WORD SEQNAME;
Static TABLE POSTABLE, CONTABLE;
Static double INDVAL[10], CODPROB[10];
Static NA INPSTRAND, COMPSTRAND;
Static LINE TITLE;
/* GLOBAL PARAMETERS */
Static long START, FINISH, WINDOW, SKIP, I, CHOICE;
Static Char WHICH, FORMAT;
Static boolean CIRCULAR;


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
/* p2c: testcode.p, line 377: 
 * Note: Line breaker spent 1.0 seconds, 5000 tries on line 592 [251] */
	if (*SEQLEN < MAXSEQ - 2) {
	  /*write(CH);*/
	  (*SEQLEN)++;
	  S_[*SEQLEN - 1] = NUC(CH);
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

/******************************************************/
/* Convert ambiguous nucleotides to N's.              */
/******************************************************/
Static Void AMBIGTON(S_, SEQLEN)
NUCLEOTIDE *S_;
long SEQLEN;
{
  long I;

  for (I = 0; I < SEQLEN; I++) {
    if ((long)S_[I] > (long)N)
      S_[I] = N;
  }
}  /* AMBIGTON */


Local Void INITPOST()
{
  PARAM *WITH;

  WITH = POSTABLE;
  WITH->LIMIT = 0.0;
  WITH->PROB[(long)A - (long)T] = 0.22;
  WITH->PROB[(long)C - (long)T] = 0.23;
  WITH->PROB[(long)G - (long)T] = 0.08;
  WITH->PROB[0] = 0.09;
  WITH = &POSTABLE[1];
  WITH->LIMIT = 1.1;
  WITH->PROB[(long)A - (long)T] = 0.20;
  WITH->PROB[(long)C - (long)T] = 0.30;
  WITH->PROB[(long)G - (long)T] = 0.08;
  WITH->PROB[0] = 0.09;
  WITH = &POSTABLE[2];
  WITH->LIMIT = 1.2;
  WITH->PROB[(long)A - (long)T] = 0.34;
  WITH->PROB[(long)C - (long)T] = 0.33;
  WITH->PROB[(long)G - (long)T] = 0.16;
  WITH->PROB[0] = 0.20;
  WITH = &POSTABLE[3];
  WITH->LIMIT = 1.3;
  WITH->PROB[(long)A - (long)T] = 0.45;
  WITH->PROB[(long)C - (long)T] = 0.51;
  WITH->PROB[(long)G - (long)T] = 0.27;
  WITH->PROB[0] = 0.54;
  WITH = &POSTABLE[4];
  WITH->LIMIT = 1.4;
  WITH->PROB[(long)A - (long)T] = 0.68;
  WITH->PROB[(long)C - (long)T] = 0.48;
  WITH->PROB[(long)G - (long)T] = 0.48;
  WITH->PROB[0] = 0.44;
  WITH = &POSTABLE[5];
  WITH->LIMIT = 1.5;
  WITH->PROB[(long)A - (long)T] = 0.58;
  WITH->PROB[(long)C - (long)T] = 0.66;
  WITH->PROB[(long)G - (long)T] = 0.53;
  WITH->PROB[0] = 0.69;
  WITH = &POSTABLE[6];
  WITH->LIMIT = 1.6;
  WITH->PROB[(long)A - (long)T] = 0.93;
  WITH->PROB[(long)C - (long)T] = 0.81;
  WITH->PROB[(long)G - (long)T] = 0.64;
  WITH->PROB[0] = 0.68;
  WITH = &POSTABLE[7];
  WITH->LIMIT = 1.7;
  WITH->PROB[(long)A - (long)T] = 0.84;
  WITH->PROB[(long)C - (long)T] = 0.70;
  WITH->PROB[(long)G - (long)T] = 0.74;
  WITH->PROB[0] = 0.91;
  WITH = &POSTABLE[8];
  WITH->LIMIT = 1.8;
  WITH->PROB[(long)A - (long)T] = 0.68;
  WITH->PROB[(long)C - (long)T] = 0.70;
  WITH->PROB[(long)G - (long)T] = 0.88;
  WITH->PROB[0] = 0.97;
  WITH = &POSTABLE[9];
  WITH->LIMIT = 1.9;
  WITH->PROB[(long)A - (long)T] = 0.94;
  WITH->PROB[(long)C - (long)T] = 0.80;
  WITH->PROB[(long)G - (long)T] = 0.90;
  WITH->PROB[0] = 0.97;
}  /* INITPOST */

/* Initialize the Content table */
Local Void INITCONT()
{
  PARAM *WITH;

  WITH = CONTABLE;
  WITH->LIMIT = 0.00;
  WITH->PROB[(long)A - (long)T] = 0.21;
  WITH->PROB[(long)C - (long)T] = 0.31;
  WITH->PROB[(long)G - (long)T] = 0.29;
  WITH->PROB[0] = 0.58;
  WITH = &CONTABLE[1];
  WITH->LIMIT = 0.17;
  WITH->PROB[(long)A - (long)T] = 0.81;
  WITH->PROB[(long)C - (long)T] = 0.39;
  WITH->PROB[(long)G - (long)T] = 0.33;
  WITH->PROB[0] = 0.51;
  WITH = &CONTABLE[2];
  WITH->LIMIT = 0.19;
  WITH->PROB[(long)A - (long)T] = 0.65;
  WITH->PROB[(long)C - (long)T] = 0.44;
  WITH->PROB[(long)G - (long)T] = 0.41;
  WITH->PROB[0] = 0.69;
  WITH = &CONTABLE[3];
  WITH->LIMIT = 0.21;
  WITH->PROB[(long)A - (long)T] = 0.67;
  WITH->PROB[(long)C - (long)T] = 0.43;
  WITH->PROB[(long)G - (long)T] = 0.41;
  WITH->PROB[0] = 0.56;
  WITH = &CONTABLE[4];
  WITH->LIMIT = 0.23;
  WITH->PROB[(long)A - (long)T] = 0.49;
  WITH->PROB[(long)C - (long)T] = 0.59;
  WITH->PROB[(long)G - (long)T] = 0.73;
  WITH->PROB[0] = 0.75;
  WITH = &CONTABLE[5];
  WITH->LIMIT = 0.25;
  WITH->PROB[(long)A - (long)T] = 0.62;
  WITH->PROB[(long)C - (long)T] = 0.59;
  WITH->PROB[(long)G - (long)T] = 0.64;
  WITH->PROB[0] = 0.55;
  WITH = &CONTABLE[6];
  WITH->LIMIT = 0.27;
  WITH->PROB[(long)A - (long)T] = 0.55;
  WITH->PROB[(long)C - (long)T] = 0.64;
  WITH->PROB[(long)G - (long)T] = 0.64;
  WITH->PROB[0] = 0.40;
  WITH = &CONTABLE[7];
  WITH->LIMIT = 0.29;
  WITH->PROB[(long)A - (long)T] = 0.44;
  WITH->PROB[(long)C - (long)T] = 0.51;
  WITH->PROB[(long)G - (long)T] = 0.47;
  WITH->PROB[0] = 0.39;
  WITH = &CONTABLE[8];
  WITH->LIMIT = 0.31;
  WITH->PROB[(long)A - (long)T] = 0.49;
  WITH->PROB[(long)C - (long)T] = 0.64;
  WITH->PROB[(long)G - (long)T] = 0.54;
  WITH->PROB[0] = 0.24;
  WITH = &CONTABLE[9];
  WITH->LIMIT = 0.33;
  WITH->PROB[(long)A - (long)T] = 0.28;
  WITH->PROB[(long)C - (long)T] = 0.82;
  WITH->PROB[(long)G - (long)T] = 0.40;
  WITH->PROB[0] = 0.28;
}  /* INITCONT */

/* Initialize the indicator arrays */
Local Void INITIND()
{
  INDVAL[0] = 0.00;
  CODPROB[0] = 0.00;
  INDVAL[1] = 0.43;
  CODPROB[1] = 0.04;
  INDVAL[2] = 0.53;
  CODPROB[2] = 0.07;
  INDVAL[3] = 0.64;
  CODPROB[3] = 0.29;
  INDVAL[4] = 0.74;
  CODPROB[4] = 0.40;
  INDVAL[5] = 0.84;
  CODPROB[5] = 0.77;
  INDVAL[6] = 0.95;
  CODPROB[6] = 0.92;
  INDVAL[7] = 1.05;
  CODPROB[7] = 0.98;
  INDVAL[8] = 1.16;
  CODPROB[8] = 1.00;
  INDVAL[9] = 1.26;
  CODPROB[9] = 1.00;
}  /* INITIND */


/******************************************************/
/* Initialize the probability tables.                 */
/******************************************************/
Static Void INITTABLES()
{
  /* Initialize the Position table */
  INITPOST();
  INITCONT();
  INITIND();
}  /* INITTABLES */


/* **************************************************** */
/* Set default values for global parameters             */
/* **************************************************** */
Static Void INITPARAM()
{
  START = 1;
  FINISH = SEQLEN;
  WHICH = 'I';
  FORMAT = 'G';
  WINDOW = 67;
  SKIP = 10;
}  /* INITPARAM */


typedef Char LETTERS[10];
typedef long CHSET[4];


/*  Read an integer parameter from the console and check */
/*    that it is in range.                               */
Local Void GETNUMBER(P, PNAME, LOW, HIGH)
long *P;
Char *PNAME;
long LOW, HIGH;
{
  printf("\nType new value for %.10s  (CURRENT VALUE: %12ld)\n", PNAME, *P);
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

  printf("Name: ");
  TEMP.f = stdout;
  *TEMP.name = '\0';
  WRITEWORD(&TEMP, SEQNAME, 20L);
  switch (CIRCULAR) {

  case true:
    printf("Topology: CIRCULAR");
    break;

  case false:
    printf("Topology:   LINEAR");
    break;
  }
  printf("%13s%10ld nt\n", "Length: ", SEQLEN);
  TEMP.f = stdout;
  *TEMP.name = '\0';
  WRITELINE(&TEMP, HLINE, 80L);
  printf("\n%12cParameter   Description/Response                 Value\n",
	 ' ');
  TEMP.f = stdout;
  *TEMP.name = '\0';
  WRITELINE(&TEMP, HLINE, 80L);
  printf("\n%12c 1)START    first nucleotide evaluated%16ld\n", ' ', START);
  printf("%12c 2)FINISH   last  nucleotide evaluated%16ld\n", ' ', FINISH);
  printf("%12c 3)WHICH    I: input strand  O: opposite strand%7c\n",
	 ' ', WHICH);
  printf("%12c 4)FORMAT   T:tabular output G:graphic output %8c\n",
	 ' ', FORMAT);
  printf("%12c 5)WINDOW   #codons in search window%18ld\n", ' ', WINDOW);
  printf("%12c 6)SKIP     #codons to skip for each window%11ld\n", ' ', SKIP);
  TEMP.f = stdout;
  *TEMP.name = '\0';
  WRITELINE(&TEMP, HLINE, 80L);
  putchar('\n');
}  /* DISPLAY */


/* **************************************************** */
/* Prompt user for parameters used by program.          */
/* **************************************************** */
Static Void PARAMETERS()
{
  long RESPONSE;
  long SET[4];
  long SET1[4];

  /* Prompt user for new parameter values */
  do {
    fprintf(stdout, "\f");
    DISPLAY();
    if (WHICH == 'O') {
      printf("Be sure START and FINISH values are appropriate for\n");
      printf("WHICH=O\n");
    }
    printf("Type number of parameter you wish to change (0 to continue)\n");
    GETINTEGER(&RESPONSE, 0L, 6L);
    scanf("%*[^\n]");
    getchar();
    if ((unsigned long)RESPONSE < 32 && ((1L << RESPONSE) & 0x7e) != 0) {
      switch (RESPONSE) {

      case 1:
	GETNUMBER(&START, "START     ", 1L, SEQLEN);
	break;

      case 2:
	GETNUMBER(&FINISH, "FINISH    ", 1L, SEQLEN);
	break;

      case 3:
	P_addset(P_expset(SET, 0L), 'I');
	GETCHAR(&WHICH, "WHICH     ", P_addset(SET, 'O'));
	break;

      case 4:
	P_addset(P_expset(SET1, 0L), 'T');
	GETCHAR(&FORMAT, "FORMAT    ", P_addset(SET1, 'G'));
	break;

      case 5:
	GETNUMBER(&WINDOW, "WINDOW    ", 10L, SEQLEN / 3);
	break;

      case 6:
	GETNUMBER(&SKIP, "SKIP      ", 1L, SEQLEN / 3);
	break;
      }
    }
  } while (RESPONSE != 0);   /* PARAMETERS */
}


/* **************************************************** */
/*   Initialize INPSTRAND and COMPSTRAND, the arrays    */
/*   which hold the valaues of the NUCLEOTIDES.         */
/* **************************************************** */
Static Void INITNUCLEOTIDES()
{
  INPSTRAND[(long)A - (long)T] = A;
  COMPSTRAND[(long)A - (long)T] = T;
  INPSTRAND[(long)C - (long)T] = C;
  COMPSTRAND[(long)C - (long)T] = G;
  INPSTRAND[(long)G - (long)T] = G;
  COMPSTRAND[(long)G - (long)T] = C;
  INPSTRAND[0] = T;
  COMPSTRAND[0] = A;
  INPSTRAND[(long)N - (long)T] = N;
  COMPSTRAND[(long)N - (long)T] = N;
}  /* INITNUCLEOTIDES */


/* Local variables for TEST: */
struct LOC_TEST {
  _TEXT *OUTFILE;
  long START, I, DISTANCE, CENTER, SENSE;
  NUCLEOTIDE NUC;
  long FREQ[(long)N - (long)T + 1][3];
  RA POS, CONT;
  double INDICATOR;
  Char OUTLINE[71], TEMPLATE[71];
} ;

Local Void SETPARAM(LINK)
struct LOC_TEST *LINK;
{
  /* Set parameters dependant on the strand */
  switch (WHICH) {

  case 'I':
    LINK->SENSE = 1;
    if (LINK->START < FINISH)
      LINK->DISTANCE = (FINISH - LINK->START + 1) / 3;
    else
      LINK->DISTANCE = (SEQLEN - LINK->START + FINISH + 1) / 3;
    break;

  case 'O':
    LINK->SENSE = -1;
    if (LINK->START > FINISH)
      LINK->DISTANCE = (LINK->START - FINISH + 1) / 3;
    else
      LINK->DISTANCE = (LINK->START + SEQLEN - FINISH + 1) / 3;
    break;
  }
}  /* SETPARAM */

/* Set up variables and print header for graph */
Local Void SETUPGRAPH(LINK)
struct LOC_TEST *LINK;
{
  long I;
  double INDEX = 0.0;

  fprintf(LINK->OUTFILE->f, "%10cWINDOW=%5ld%15s%5ld\n\n",
	  ' ', WINDOW, "SKIP=", SKIP);
  fprintf(LINK->OUTFILE->f, "%35s%22s%14s\n\n",
	  "NON-CODING", "NO OPINION", "CODING");
  for (I = 0; I <= 7; I++) {
    fprintf(LINK->OUTFILE->f, "%10.1f", INDEX);
    INDEX += 0.2;
  }
  fprintf(LINK->OUTFILE->f, "\n%10c", ' ');
  for (I = 1; I <= 7; I++)
    fprintf(LINK->OUTFILE->f, "---------+");
  putc('\n', LINK->OUTFILE->f);

  /* Initialize line templates */
  for (I = 1; I <= 70; I++)
    LINK->TEMPLATE[I] = ' ';
  LINK->TEMPLATE[0] = '|';
  LINK->TEMPLATE[37] = '|';
  LINK->TEMPLATE[47] = '|';
  memcpy(LINK->OUTLINE, LINK->TEMPLATE, 71L);
}  /* SETUPGRAPH */

/* Compute the next nucleotide position */
Local long NEXT(POS, LINK)
long POS;
struct LOC_TEST *LINK;
{
  POS += LINK->SENSE;
  if (CIRCULAR) {
    if (POS > SEQLEN)
      return 1;
    else if (POS < 1)
      return SEQLEN;
    else
      return POS;
  } else if (POS > SEQLEN || POS < 1)
    return 0;
  else
    return POS;
}  /* NEXT */

/* Local variables for INTERMEDIATE: */
struct LOC_INTERMEDIATE {
  struct LOC_TEST *LINK;
} ;

Local Void PRINTARRAY(AR, LINK)
double *AR;
struct LOC_INTERMEDIATE *LINK;
{
  NUCLEOTIDE NUC;

  for (NUC = T; (long)NUC <= (long)G; NUC = (NUCLEOTIDE)((long)NUC + 1))
    fprintf(LINK->LINK->OUTFILE->f, "%10.2f", AR[(long)NUC - (long)T]);
  putc('\n', LINK->LINK->OUTFILE->f);
}  /* PRINTARRAY */

/*Print the intermediate results */
Local Void INTERMEDIATE(LINK)
struct LOC_TEST *LINK;
{
  struct LOC_INTERMEDIATE V;
  long I;
  NUCLEOTIDE NUC;

  V.LINK = LINK;
  fprintf(LINK->OUTFILE->f, "\nOpen reading frame from %12ld to %12ld\n",
	  LINK->START, FINISH);
  fprintf(LINK->OUTFILE->f, "%20c%10c%10c%10c\n", 'T', 'C', 'A', 'G');
  fprintf(LINK->OUTFILE->f, "Pos.Freq.\n");
  for (I = 1; I <= 3; I++) {
    fprintf(LINK->OUTFILE->f, "%10ld", I);
    for (NUC = T; (long)NUC <= (long)G; NUC = (NUCLEOTIDE)((long)NUC + 1))
      fprintf(LINK->OUTFILE->f, "%10ld", LINK->FREQ[(long)NUC - (long)T][I-1]);
    putc('\n', LINK->OUTFILE->f);
  }
  for (I = 1; I <= 60; I++)
    putc('-', LINK->OUTFILE->f);
  fprintf(LINK->OUTFILE->f, "\nCont.Param");
  PRINTARRAY(LINK->CONT, &V);
  fprintf(LINK->OUTFILE->f, "Posn.Param");
  PRINTARRAY(LINK->POS, &V);
}  /* INTERMEDIATE */

/* Calculate the maximum of three integers */
Local long MAX(I1, I2, I3, LINK)
long I1, I2, I3;
struct LOC_TEST *LINK;
{
  long TEMP;

  if (I1 > I2)
    TEMP = I1;
  else
    TEMP = I2;
  if (I3 > TEMP)
    return I3;
  else
    return TEMP;
}  /* MAX */

/* Calculate the minimum of three integers */
Local long MIN(I1, I2, I3, LINK)
long I1, I2, I3;
struct LOC_TEST *LINK;
{
  long TEMP;

  if (I1 < I2)
    TEMP = I1;
  else
    TEMP = I2;
  if (I3 < TEMP)
    return I3;
  else
    return TEMP;
}  /* MIN */

/* Look up the probability of a given parameter in the table specified */
Local Void FINDPROB(P, T_, LINK)
double *P;
PARAM *T_;
struct LOC_TEST *LINK;
{
  long I = 10;

  while (*P < T_[I-1].LIMIT)
    I--;
  *P = T_[I-1].PROB[(long)LINK->NUC - (long)T];
}  /* FINDPROB */

/* Print the results */
Local Void RESULTS(LINK)
struct LOC_TEST *LINK;
{
  long NUMCHAR, J;

  switch (FORMAT) {

  case 'T':
    fprintf(LINK->OUTFILE->f, "TESTCODE indicator: % .5E\n", LINK->INDICATOR);
    fprintf(LINK->OUTFILE->f, "Probability of coding: % .5E\n",
	    CODPROB[LINK->I-1]);
    fprintf(LINK->OUTFILE->f, "Prediction: \n");
    switch (LINK->I) {

    case 1:
    case 2:
    case 3:
    case 4:
      fprintf(LINK->OUTFILE->f, "NONCODING\n");
      break;

    case 5:
    case 6:
      fprintf(LINK->OUTFILE->f, "NO OPINION\n");
      break;

    case 7:
    case 8:
    case 9:
    case 10:
      fprintf(LINK->OUTFILE->f, "CODING\n");
      break;
    }
    break;

  case 'G':
    NUMCHAR = (long)floor(LINK->INDICATOR / 1.4 * LINESIZE + 0.5);
    if (NUMCHAR > LINESIZE)
      NUMCHAR = LINESIZE;
    fprintf(LINK->OUTFILE->f, "%9ld", LINK->CENTER);
    for (J = 1; J <= NUMCHAR; J++)
      LINK->OUTLINE[J] = '=';
    for (J = 0; J <= LINESIZE; J++)
      putc(LINK->OUTLINE[J], LINK->OUTFILE->f);
    putc('\n', LINK->OUTFILE->f);
    break;
  }
}  /* RESULTS */


/* *********************************************************** */
/*  Test an open reading frame specified by the user to        */
/*  to determine whether or not it is a protein coding region  */
/*  using the algorithm found in:                              */
/*   Fickett,James, "Recognition of protein coding regions in  */
/*   DNA sequences", Nucleic Acids Research 10 No.17,5303-5318.*/
/* *********************************************************** */
Static Void TEST(OUTFILE_, STRAND, START_)
_TEXT *OUTFILE_;
NUCLEOTIDE *STRAND;
long START_;
{
  struct LOC_TEST V;
  long POSITION, CODONSDONE, NUMCYCLES;
  long CYCLESDONE = 0;
  long S2, S3, FORLIM;

  /* Print the header */
  V.OUTFILE = OUTFILE_;
  V.START = START_;
  SETPARAM(&V);

  /* Only test regions >= MINCODONS */
  if (V.DISTANCE < MINCODONS) {
    printf("Reading frame must be >= %12ld codons (%12ld bp)\n",
	   (long)MINCODONS, MINCODONS * 3L);
    return;
  }
  fprintf(V.OUTFILE->f, "\n\n%45s\n\n", VERSION);
  fprintf(V.OUTFILE->f, "%10c", ' ');
  WRITELINE(V.OUTFILE, TITLE, TITLE.LEN);
  putc('\n', V.OUTFILE->f);

  if (FORMAT == 'G') {
    SETUPGRAPH(&V);
    NUMCYCLES = (V.DISTANCE - WINDOW) / SKIP + 1;
  } else {
    WINDOW = V.DISTANCE;
    NUMCYCLES = 1;
  }

  while (CYCLESDONE < NUMCYCLES) {
    /* Initialize the frequency table */
    for (V.NUC = T; (long)V.NUC <= (long)N; V.NUC = (NUCLEOTIDE)((long)V.NUC + 1)) {
      for (V.I = 1; V.I <= 3; V.I++)
	V.FREQ[(long)V.NUC - (long)T][V.I-1] = 0;
    }

    /* For the reading frame specified, calculate the sums of */
    /* Ai, Gi, Ci, Ti for positions i= 1 to 3 in the codons   */
    POSITION = V.START;
    V.CENTER = WINDOW / 2;
    CODONSDONE = 0;
    while (CODONSDONE < WINDOW) {
      S2 = NEXT(POSITION, &V);
      S3 = NEXT(S2, &V);
      V.FREQ[(long)STRAND[(long)SEQ[POSITION-1] - (long)T] - (long)T][0]++;
      V.FREQ[(long)STRAND[(long)SEQ[S2-1] - (long)T] - (long)T][1]++;
      V.FREQ[(long)STRAND[(long)SEQ[S3-1] - (long)T] - (long)T][2]++;
      /* CENTER determines coordinate for center of winwow */
      CODONSDONE++;
      if (CODONSDONE == V.CENTER)
	V.CENTER = POSITION;
      POSITION = NEXT(S3, &V);
    }  /* CODONSDONE < WINDOW */

    /* Calculate the position and content parameters */
    /* Content =  % A,C,G,T                          */
    /*                   MAX(A1,A2,A3)               */
    /* A-Posn. =        -----------------     etc.   */
    /*                   MIN(A1,A2,A3) + 1           */
    for (V.NUC = T; (long)V.NUC <= (long)G; V.NUC = (NUCLEOTIDE)((long)V.NUC + 1)) {
      V.CONT[(long)V.NUC - (long)T] = (V.FREQ[(long)V.NUC - (long)T]
				       [0] + V.FREQ[(long)V.NUC - (long)T]
				       [1] + V.FREQ[(long)V.NUC - (long)T]
				       [2]) / (V.DISTANCE * 3.0);
      V.POS[(long)V.NUC - (long)T] = MAX(V.FREQ[(long)V.NUC - (long)T]
	    [0], V.FREQ[(long)V.NUC - (long)T]
	    [1], V.FREQ[(long)V.NUC - (long)T][2],
	    &V) / (MIN(V.FREQ[(long)V.NUC - (long)T]
		       [0], V.FREQ[(long)V.NUC - (long)T]
		       [1], V.FREQ[(long)V.NUC - (long)T][2], &V) + 1.0);
    }

    if (FORMAT == 'T')   /*Print intermediate results */
      INTERMEDIATE(&V);

    /* Look up the probabilities in Table 1.  Substitute the eight  */
    /* parameters with their corresponding probabilities            */
    for (V.NUC = T; (long)V.NUC <= (long)G; V.NUC = (NUCLEOTIDE)((long)V.NUC + 1)) {
      FINDPROB(&V.POS[(long)V.NUC - (long)T], POSTABLE, &V);
      FINDPROB(&V.CONT[(long)V.NUC - (long)T], CONTABLE, &V);
    }

    /* Multiply by weighting factors in Table 2 */
    V.POS[(long)A - (long)T] *= 0.26;
    V.CONT[(long)A - (long)T] *= 0.11;
    V.POS[(long)C - (long)T] *= 0.18;
    V.CONT[(long)C - (long)T] *= 0.12;
    V.POS[(long)G - (long)T] *= 0.31;
    V.CONT[(long)G - (long)T] *= 0.15;
    V.POS[0] *= 0.33;
    V.CONT[0] *= 0.14;

    /* Calculate the value of the TESTCODE indicator as */
    /* the sum of the weighted products.                */
    V.INDICATOR = 0.0;
    for (V.NUC = T; (long)V.NUC <= (long)G; V.NUC = (NUCLEOTIDE)((long)V.NUC + 1))
      V.INDICATOR += V.POS[(long)V.NUC - (long)T] + V.CONT[(long)V.NUC - (long)T];
    V.I = 10;
    while (V.INDICATOR < INDVAL[V.I-1])
      V.I--;

    /* Print the results */
    RESULTS(&V);

    CYCLESDONE++;
    memcpy(V.OUTLINE, V.TEMPLATE, 71L);
    FORLIM = SKIP * 3;
    for (V.I = 1; V.I <= FORLIM; V.I++)
      V.START = NEXT(V.START, &V);
  }  /* while CYCLESDONE < NUMCYCLES */
}  /* TEST */


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
  printf("Type sequence filename:\n");
  GETFILE(&INFILE, 'I', &IFN);
  READSEQ(&INFILE, SEQ, &SEQLEN, &SEQNAME, &CIRCULAR);
  printf("Type output filename:\n");
  GETFILE(&OUTFILE, 'O', &OFN);
  putc('\n', OUTFILE.f);
  putchar('\n');
  AMBIGTON(SEQ, SEQLEN);
  INITPARAM();
  INITTABLES();
  INITNUCLEOTIDES();

  /* Initialize horizontal output line */

  for (I = 1; I <= MAXLINE; I++)
    HLINE.STR[I-1] = '_';
  HLINE.LEN = MAXLINE;
  /* MAIN LOOP */
  do {
    putchar('\n');
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, HLINE, 80L);
    printf("\nTESTCODE%30s\n", "MAIN MENU");
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
    printf("\nTitle:             ");
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, TITLE, TITLE.LEN);
    putchar('\n');
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, HLINE, 80L);
    printf("\n%20c1) Read in a new sequence\n", ' ');
    printf("%20c2) Open a new output file\n", ' ');
    printf("%20c3) Type in a title line for output\n", ' ');
    printf("%20c4) Change parameters\n", ' ');
    printf("%20c5) Search sequence (output to screen)\n", ' ');
    printf("%20c6) Search sequence (output to file)\n", ' ');
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, HLINE, 80L);
    printf("\nType the number of your choice  (0 to quit program)\n");
    GETINTEGER(&CHOICE, 0L, 6L);
    scanf("%*[^\n]");
    getchar();

    switch (CHOICE) {

    case 0:
      /* blank case */
      break;

    case 1:
      printf("Enter sequence filename:\n");
      GETFILE(&INFILE, 'I', &IFN);
      READSEQ(&INFILE, SEQ, &SEQLEN, &SEQNAME, &CIRCULAR);
      AMBIGTON(SEQ, SEQLEN);
      INITPARAM();
      break;

    case 2:
      if (OUTFILE.f != NULL)
	fclose(OUTFILE.f);
      OUTFILE.f = NULL;
      printf("Type output filename:\n");
      GETFILE(&OUTFILE, 'O', &OFN);
      break;

    /*!!!*/
    case 3:
      printf("Type a title to appear in output (<RETURN> for blank)\n");
      INPLINE(&TITLE);
      break;

    case 4:
      PARAMETERS();
      break;

    case 5:
      switch (WHICH) {

      case 'I':
	TEMP.f = stdout;
	*TEMP.name = '\0';
	TEST(&TEMP, INPSTRAND, START);
	break;

      case 'O':
	TEMP.f = stdout;
	*TEMP.name = '\0';
	TEST(&TEMP, COMPSTRAND, START);
	break;
      }
      printf("\nPress RETURN to continue");
      scanf("%*[^\n]");
      getchar();
      break;

    case 6:
      switch (WHICH) {

      case 'I':
	TEST(&OUTFILE, INPSTRAND, START);
	break;

      case 'O':
	TEST(&OUTFILE, COMPSTRAND, START);
	break;
      }
      fprintf(OUTFILE.f, "\n\n");
      printf("Press RETURN to continue");
      scanf("%*[^\n]");
      getchar();
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
}  /* TESTCODE */




/* End. */
