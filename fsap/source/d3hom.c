/* Output from p2c 1.21alpha-07.Dec.93, the Pascal-to-C translator */
/* From input file "d3hom.p" */


/* ********************************************************  */
/*                                                           */
/*   D3HOM     VERSION  8/13/2001, Standard Pascal           */
/*             Brian Fristensky                              */
/*             Dept. of Plant Science                        */
/*             University of Manitoba                        */
/*             Winnipeg, MB R3T 2N2   CANADA                 */
/*                                                           */
/* Copyright (c) 1984,1986,1987, 1990 by Brian Fristensky.   */
/* !!! in comment indicates feature which may need change.   */
/*  *******************************************************  */

/*!!!*/

#include <p2c.h>


/*, SFILEX,SFILEY,OUTFILE*/
/*!!! Some Pascals require file parameters in heading */

#define MAXSEQ          750000L
/* BEGIN MODULE REALLIMITS */
/*!!!*/

#define MAXREAL         1.7e38
/* END MODULE REALLIMITS         VERSION= 'SUNMODS     Version  6/26/01'; */

#define MAXWORD         25
#define MAXRANGE        30
#define MAXLINE         150

#define VERSION         "D3HOM         Version  8/13/2001"


typedef enum {
  T, C, A, G, N, R, Y, M, W, S, K, D, H, V, B, Z
} NUCLEOTIDE;
typedef NUCLEOTIDE SEQUENCE[MAXSEQ + MAXRANGE + 1];

typedef struct NODE {
  long POS;
  struct NODE *NEXT;
} NODE;

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


Static _TEXT SFILEX, SFILEY, OUTFILE;
Static LINE XFN, YFN, OFN, HLINE;
Static SEQUENCE SEQX, SEQY;
Static long LENX, LENY;
Static WORD NAMEX, NAMEY;
Static NODE *TRIPLET[(long)G - (long)T + 1][(long)G - (long)T + 1]
	    [(long)G - (long)T + 1];
Static NODE *FREELIST;
Static Char NUCHAR[16];
Static long STARTX, FINISHX, STARTY, FINISHY, HOMRANGE, MINPER, COMPFACT,
	    LINESIZE, MAXSCORE, THREEMATCH, MINSCORE;
Static double SCALEFAC, MAXFACTOR, TEMP;
Static long SUBSCORE[MAXRANGE + 1];
Static Char KIND;
Static Char PC[102];   /*Printing characters for graph*/
Static boolean INVERSEX, CIRCX, CIRCY;
Static long CHOICE, J;


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


/*!!!*/

#define HIGHEXP         38
/*!!!*/
#define LOWEXP          (-38)


typedef char EXPONENT;




Local double TEN(E)
EXPONENT E;
{
  /* = 10**E */
  long I = 0;
  double T_ = 1.0;

  do {
    if (E & 1) {
      switch (I) {

      case 0:
	T_ *= 1.0e1;
	break;

      case 1:
	T_ *= 1.0e2;
	break;

      case 2:
	T_ *= 1.0e4;
	break;

      case 3:
	T_ *= 1.0e8;
	break;

      case 4:
	T_ *= 1.0e16;
	break;

      case 5:
	T_ *= 1.0e32;
	break;
	/*!!!   6: T:= T * 1.0E64; Max. exponent is 38
	        7: T:= T * 1.0E128;
	        8: T:= T * 1.0E256; */
      }
    }
    E /= 2;
    I++;
  } while (E != 0);
  return T_;
}  /* TEN */


/* END MODULE WRITEWORD         VERSION= 'SUNMODS     Version  6/26/01'; */

/* BEGIN MODULE GETREAL */
/* The procedure GETREAL has been adapted from the procedure 'rdr', from
   Jensen,K., and Wirth,N. (1974) Pascal User Manual and Report, 2nd
   Ed., Springer-Verlag, pp122-123.  The scope of real numbers considered
   legal by GETREAL includes all possible real numbers as defined on p111
   with two additional allowances: (i) numbers may begin with a decimal
   point (ii) 'E' or 'e' may be used to indicate exponentiation. */

Static Void GETREAL(VAL, LOW, HIGH)
double *VAL, LOW, HIGH;
{
  Char CH;
  double LIMIT = MAXREAL / 10;
  double Y_;
  long ORDZERO = '0';
  long I, E;
  boolean LEGAL, INRANGE;
  WORD NUMWORD;
  boolean NEG1, NEG2;
  double A_, S_;
  long LCOUNT, RCOUNT;
  long DIGITS[3];

  P_addsetr(P_expset(DIGITS, 0L), '0', '9');
  do {
    do {
      LEGAL = true;
      INRANGE = true;
      INPWORD(&NUMWORD);
      /* Evaluate sign, if any */
      I = 1;
      CH = NUMWORD.STR[I-1];
      if (CH == '-') {
	NEG1 = true;
	I++;
	CH = NUMWORD.STR[I-1];
      } else {
	NEG1 = false;
	if (CH == '+') {
	  I++;
	  CH = NUMWORD.STR[I-1];
	}
      }
      if (!(CH == '.' || isdigit(CH)))
	goto _L777;

      /* Evaluate whole number part (optional) */
      A_ = 0.0;
      E = 0;
      LCOUNT = 1;
      while (P_inset(CH, DIGITS)) {
	if (A_ < LIMIT)
	  A_ = 10 * A_ + CH - ORDZERO;
	else
	  E++;
	I++;
	CH = NUMWORD.STR[I-1];
	LCOUNT++;
      }

      /* Evaluate fractional part. */
      RCOUNT = 0;
      if (CH == '.') {
	I++;
	CH = NUMWORD.STR[I-1];
	while (P_inset(CH, DIGITS)) {
	  if (A_ < LIMIT) {
	    A_ = 10 * A_ + CH - ORDZERO;
	    E--;
	  }
	  I++;
	  CH = NUMWORD.STR[I-1];
	  RCOUNT++;
	}
      }

      /* Evaluate exponent */
      if (CH == 'E' || CH == 'e') {
	I++;
	CH = NUMWORD.STR[I-1];
	S_ = 0.0;
	if (CH == '-') {
	  NEG2 = true;
	  I++;
	  CH = NUMWORD.STR[I-1];
	} else {
	  NEG2 = false;
	  if (CH == '+') {
	    I++;
	    CH = NUMWORD.STR[I-1];
	  }
	}
	if (!P_inset(CH, DIGITS))
	  goto _L777;
	S_ = CH - ORDZERO;
	I++;
	CH = NUMWORD.STR[I-1];
	while (P_inset(CH, DIGITS)) {
	  if (S_ < LIMIT)
	    S_ = 10 * S_ + CH - ORDZERO;
	  I++;
	  CH = NUMWORD.STR[I-1];
	}
	if (NEG2)
	  E -= (long)floor(S_ + 0.5);
	else
	  E += (long)floor(S_ + 0.5);
      }

      /* Check for errors  */
      if (I <= NUMWORD.LEN)   /*illegal char*/
	goto _L777;
      if (E - RCOUNT < LOWEXP || E + LCOUNT > HIGHEXP)
	goto _L776;

      /* Calculate final value */
      Y_ = A_;
      if (NEG1)
	Y_ = -Y_;
      if (E < 0)
	*VAL = Y_ / TEN((int)(-E));
      else if (E > 0)
	*VAL = Y_ * TEN((int)E);
      else
	*VAL = Y_;
      goto _L778;

      /* Error handling statements */
_L776:
      LEGAL = false;
      printf("\n Exponent Error.\n");
      printf(" Enter new value:\n");
      goto _L778;
_L777:
      LEGAL = false;
      printf("\n Illegal character encountered.\n");
      printf("Enter new value:");
_L778: ;   /* no error, do nothing */
    } while (!LEGAL);


    /* Make sure number is in range */
    if (*VAL < LOW || *VAL > HIGH) {
      INRANGE = false;
      printf("\n Number is out of range\n");
      printf(" Enter new value:\n");
    } else
      INRANGE = true;
  } while (!INRANGE);   /* GETREAL */
}

#undef HIGHEXP
#undef LOWEXP


/* END MODULE GETREAL         VERSION= 'SUNMODS     Version  6/26/01'; */

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
/* p2c: d3hom.p, line 544: 
 * Note: Line breaker spent 0.0 seconds, 5000 tries on line 842 [251] */
	if (*SEQLEN < MAXSEQ - 2) {
	  /*write(CH);*/
	  (*SEQLEN)++;
	  S_[*SEQLEN + MAXRANGE] = NUC(CH);
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

/*****************************************************/
/* Initialization procedures.                        */
/*****************************************************/
Static Void INITIALIZE()
{
  NUCLEOTIDE N1, N2, N3;

  /*Set default values for parameters */
  HOMRANGE = 10;
  SCALEFAC = 0.95;
  MINPER = 60;
  COMPFACT = 10;
  KIND = 'D';
  LINESIZE = 70;

  /*   Initialize NUCHAR array which holds the            */
  /*   character values of the NUCLEOTIDES.               */
  NUCHAR[(long)A] = 'A';
  NUCHAR[(long)C] = 'C';
  NUCHAR[(long)G] = 'G';
  NUCHAR[(long)T] = 'T';
  NUCHAR[(long)R] = 'R';
  NUCHAR[(long)Y] = 'Y';
  NUCHAR[(long)M] = 'M';
  NUCHAR[(long)W] = 'W';
  NUCHAR[(long)S] = 'S';
  NUCHAR[(long)K] = 'K';
  NUCHAR[(long)D] = 'D';
  NUCHAR[(long)H] = 'H';
  NUCHAR[(long)V] = 'V';
  NUCHAR[(long)B] = 'B';
  NUCHAR[(long)N] = 'N';
  NUCHAR[(long)Z] = 'N';

  /* Initialize the PC array which holds the       */
  /* characers which symbolize percent identity.   */
  PC[0] = ' ';
  PC[1] = '.';
  PC[26] = 'm';
  PC[27] = 'l';
  PC[28] = 'l';
  PC[29] = 'k';
  PC[30] = 'k';
  PC[31] = 'j';
  PC[32] = 'j';
  PC[33] = 'i';
  PC[34] = 'i';
  PC[35] = 'h';
  PC[36] = 'h';
  PC[37] = 'g';
  PC[38] = 'g';
  PC[39] = 'f';
  PC[40] = 'f';
  PC[41] = 'e';
  PC[42] = 'e';
  PC[43] = 'd';
  PC[44] = 'd';
  PC[45] = 'c';
  PC[46] = 'c';
  PC[47] = 'b';
  PC[48] = 'b';
  PC[49] = 'a';
  PC[50] = 'a';
  PC[51] = 'Z';
  PC[52] = 'Z';
  PC[53] = 'Y';
  PC[54] = 'Y';
  PC[55] = 'X';
  PC[56] = 'X';
  PC[57] = 'W';
  PC[58] = 'W';
  PC[59] = 'V';
  PC[60] = 'V';
  PC[61] = 'U';
  PC[62] = 'U';
  PC[63] = 'T';
  PC[64] = 'T';
  PC[65] = 'S';
  PC[66] = 'S';
  PC[67] = 'R';
  PC[68] = 'R';
  PC[69] = 'Q';
  PC[70] = 'Q';
  PC[71] = 'P';
  PC[72] = 'P';
  PC[73] = 'O';
  PC[74] = 'O';
  PC[75] = 'N';
  PC[76] = 'N';
  PC[77] = 'M';
  PC[78] = 'M';
  PC[79] = 'L';
  PC[80] = 'L';
  PC[81] = 'K';
  PC[82] = 'K';
  PC[83] = 'J';
  PC[84] = 'J';
  PC[85] = 'I';
  PC[86] = 'I';
  PC[87] = 'H';
  PC[88] = 'H';
  PC[89] = 'G';
  PC[90] = 'G';
  PC[91] = 'F';
  PC[92] = 'F';
  PC[93] = 'E';
  PC[94] = 'E';
  PC[95] = 'D';
  PC[96] = 'D';
  PC[97] = 'C';
  PC[98] = 'C';
  PC[99] = 'B';
  PC[100] = 'B';
  PC[101] = 'A';

  /* Initialize TRIPLET table to nil */
  for (N1 = T; (long)N1 <= (long)G; N1 = (NUCLEOTIDE)((long)N1 + 1)) {
    for (N2 = T; (long)N2 <= (long)G; N2 = (NUCLEOTIDE)((long)N2 + 1)) {
      for (N3 = T; (long)N3 <= (long)G; N3 = (NUCLEOTIDE)((long)N3 + 1))
	TRIPLET[(long)N1 - (long)T][(long)N2 - (long)T]
	  [(long)N3 - (long)T] = NULL;
    }
  }
}  /* INITIALIZE */


/* Initialize negative and last parts of the sequence*/
Local Void SEQENDS(S_, LEN, CIRCULAR, UNKNOWN)
NUCLEOTIDE *S_;
long *LEN;
boolean *CIRCULAR;
NUCLEOTIDE UNKNOWN;
{
  long I, J, FORLIM;

  if (*CIRCULAR) {
    J = *LEN - MAXRANGE - 1;
    for (I = -MAXRANGE; I <= 0; I++) {
      S_[I + MAXRANGE] = S_[J + MAXRANGE];
      J++;
    }
    J = 1;
    FORLIM = *LEN + MAXRANGE;
    for (I = *LEN + 1; I <= FORLIM; I++) {
      S_[I + MAXRANGE] = S_[J + MAXRANGE];
      J++;
    }
    return;
  }
  for (I = -MAXRANGE; I <= 0; I++)
    S_[I + MAXRANGE] = UNKNOWN;
  FORLIM = *LEN + MAXRANGE;
  for (I = *LEN + 1; I <= FORLIM; I++)
    S_[I + MAXRANGE] = UNKNOWN;
}  /* SEQENDS */


/* **************************************************** */
/* Open a sequence file and read in sequence.           */
/* **************************************************** */
Static Void NEWSEQ(SFILE, FN, S_, LEN, START, FINISH, NAME, CIRC, AXIS)
_TEXT *SFILE;
LINE *FN;
NUCLEOTIDE *S_;
long *LEN, *START, *FINISH;
WORD *NAME;
boolean *CIRC;
Char AXIS;
{
  Char ANSWER;
  long J, FORLIM;

  printf("Enter %c-axis sequence filename:\n", AXIS);
  GETFILE(SFILE, 'I', FN);
  NAME->LEN = 0;
  READSEQ(SFILE, S_, LEN, NAME, CIRC);
  if (NAME->LEN == 0) {
    printf("Type name for %c-axis sequence to appear on output:\n", AXIS);
    INPWORD(NAME);
    scanf("%*[^\n]");
    getchar();
  }
  *START = 1;
  *FINISH = *LEN;
  switch (AXIS) {

  case 'X':
    SEQENDS(S_, LEN, CIRC, Z);
    FORLIM = *LEN;
    for (J = 1; J <= FORLIM; J++) {
      if (S_[J + MAXRANGE] == N)
	S_[J + MAXRANGE] = Z;
    }
    do {
      printf("Is seq. on X-axis the inverse strand? [Y/N]\n");
      scanf("%c%*[^\n]", &ANSWER);
      getchar();
      if (ANSWER == '\n')
	ANSWER = ' ';
    } while (ANSWER != 'n' && ANSWER != 'N' && ANSWER != 'y' && ANSWER != 'Y');
    if (ANSWER == 'y' || ANSWER == 'Y') {
      INVERSEX = true;
      *START = *LEN;
      *FINISH = 1;
    } else
      INVERSEX = false;
    break;

  case 'Y':
    SEQENDS(S_, LEN, CIRC, N);
    break;
  }
}  /* NEWSEQ */


typedef Char LETTERS[10];
typedef long CHSET[4];



/*  Read an integer parameter from the console and check */
/*    that it is in range.                               */
Local Void GETNUMBER(P, PNAME, LOW, HIGH)
long *P;
Char *PNAME;
long LOW, HIGH;
{
  double TEMP;

  printf("\nType new value for %.10s  (CURRENT VALUE: %12ld)\n", PNAME, *P);
  GETREAL(&TEMP, (double)LOW, (double)HIGH);
  scanf("%*[^\n]");
  getchar();
  *P = (long)floor(TEMP + 0.5);
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

Local Void HEADLINE(NAME, CIRC, INVERSE, LEN, AXIS)
WORD NAME;
boolean CIRC, INVERSE;
long LEN;
Char AXIS;
{
  _TEXT TEMP1;

  printf("%c-axis: ", AXIS);
  TEMP1.f = stdout;
  *TEMP1.name = '\0';
  WRITEWORD(&TEMP1, NAME, 20L);
  switch (CIRC) {

  case true:
    printf("Topology: CIRCULAR");
    break;

  case false:
    printf("Topology:   LINEAR");
    break;
  }
  switch (INVERSE) {

  case true:
    printf(" (INVERSE)");
    break;

  case false:
    printf("          ");
    break;
  }
  printf("%10s%10ld nt\n", "Length: ", LEN);
}  /* HEADLINE */

/* END MODULE GETCHAR         VERSION= 'SUNMODS     Version  6/26/01'; */

/* Display  parameters on screen */
Local Void DISPLAY()
{
  _TEXT TEMP1;

  fprintf(stdout, "\f");
  HEADLINE(NAMEX, CIRCX, INVERSEX, LENX, 'X');
  HEADLINE(NAMEY, CIRCY, false, LENY, 'Y');
  TEMP1.f = stdout;
  *TEMP1.name = '\0';
  WRITELINE(&TEMP1, HLINE, 80L);
  printf("\n%66s\n", "Parameter   Description/Response                 Value");
  TEMP1.f = stdout;
  *TEMP1.name = '\0';
  WRITELINE(&TEMP1, HLINE, 80L);
  printf("\n%60s%6ld\n",
	 " 1)STARTX   first nucleotide position in SEQX   ", STARTX);
  printf("%60s%6ld\n",
	 " 2)FINISHX  last  nucleotide position in SEQX   ", FINISHX);
  printf("%60s%6ld\n",
	 " 3)STARTY   first nucleotide position in SEQY   ", STARTY);
  printf("%60s%6ld\n",
	 " 4)FINISHY  last  nucleotide position in SEQY   ", FINISHY);
  printf("%60s%6ld\n",
	 " 5)HOMRANGE dist.from central triplet in a match", HOMRANGE);
  printf("%60s%6.2f\n",
	 " 6)SCALEFAC scale factor for exponential curve  ", SCALEFAC);
  printf("%60s%6ld\n",
	 " 7)MINPER   minimum percent similarity printed  ", MINPER);
  printf("%60s%6ld\n",
	 " 8)COMPFACT graph compression factor            ", COMPFACT);
  printf("%60s%6c\n",
	 " 9)KIND     D:DNA            R:RNA              ", KIND);
  printf("%60s%6ld\n",
	 "10)LINESIZE width of output line (eg.70,120)    ", LINESIZE);
  TEMP1.f = stdout;
  *TEMP1.name = '\0';
  WRITELINE(&TEMP1, HLINE, 80L);
  printf("\nType number of parameter you wish to change (0 to continue)\n");
}  /* DISPLAY */


/* **************************************************** */
/* Prompt user for parameters used by program.          */
/* **************************************************** */
Static Void PARAMETERS()
{
  long RESPONSE;
  double TEMP;
  long SET[4];

  /* Prompt user for new parameter values */
  do {
    DISPLAY();
    GETREAL(&TEMP, 0.0, 10.0);
    scanf("%*[^\n]");
    getchar();
    RESPONSE = (long)floor(TEMP + 0.5);
    if ((unsigned long)RESPONSE < 32 && ((1L << RESPONSE) & 0x7fe) != 0) {
      switch (RESPONSE) {

      case 1:
	GETNUMBER(&STARTX, "STARTX    ", 1L, LENX);
	break;

      case 2:
	GETNUMBER(&FINISHX, "FINISHX   ", 1L, LENX);
	break;

      case 3:
	GETNUMBER(&STARTY, "STARTY    ", 1L, LENY);
	break;

      case 4:
	GETNUMBER(&FINISHY, "FINISHY   ", 1L, LENY);
	break;

      case 5:
	GETNUMBER(&HOMRANGE, "HOMRANGE  ", 1L, (long)MAXRANGE);
	break;

      case 6:
	printf("Type new value for SCALEFAC:\n");
	GETREAL(&SCALEFAC, 0.0, 1.0);
	scanf("%*[^\n]");
	getchar();
	break;

      case 7:
	GETNUMBER(&MINPER, "MINPER    ", 40L, 100L);
	break;

      case 8:
	GETNUMBER(&COMPFACT, "COMPFACT  ", 1L, 500L);
	break;

      case 9:
	P_addset(P_expset(SET, 0L), 'D');
	GETCHAR(&KIND, "KIND      ", P_addset(SET, 'R'));
	break;

      case 10:
	GETNUMBER(&LINESIZE, "LINESIZE  ", 40L, (long)MAXLINE);
	break;
      }
    }
  } while (RESPONSE != 0);
  switch (KIND) {

  case 'D':
    NUCHAR[(long)T] = 'T';
    break;

  case 'R':
    NUCHAR[(long)T] = 'U';
    break;
  }
}  /* PARAMETERS */


/****************************************************************/
/* Calculate a table of scores for matches at each distance from*/
/* the central triplet in a local homology.                     */
/****************************************************************/
Static Void CALCSCORES()
{
  long I;
  long V_ = 100;
  double S_;
  long FORLIM;

  SUBSCORE[0] = V_;
  MAXSCORE = V_;
  S_ = SCALEFAC;
  FORLIM = HOMRANGE;
  for (I = 1; I <= FORLIM; I++) {
    SUBSCORE[I] = (long)floor(V_ * S_ + 0.5);
    MAXSCORE += SUBSCORE[I] * 2;
    S_ *= SCALEFAC;
  }
  THREEMATCH = SUBSCORE[0] + SUBSCORE[1] * 2;
  MINSCORE = (long)floor(MINPER / 100.0 * MAXSCORE + 0.5);
  MAXFACTOR = 100.0 / MAXSCORE;
}  /*CALCSCORES*/


/* Local variables for QUICKSEARCH: */
struct LOC_QUICKSEARCH {
  _TEXT *OUTFILE;
  long POSX, POSY, LEFT;
  short THISLINE[MAXLINE];
} ;

/* Print the header listing graph parameters. */
Local Void HEADER(LINK)
struct LOC_QUICKSEARCH *LINK;
{
  fprintf(LINK->OUTFILE->f, "%s\n", VERSION);
  fprintf(LINK->OUTFILE->f, "X-axis: ");
  WRITEWORD(LINK->OUTFILE, NAMEX, NAMEX.LEN);
  fprintf(LINK->OUTFILE->f, "\nY-axis: ");
  WRITEWORD(LINK->OUTFILE, NAMEY, NAMEY.LEN);
  fprintf(LINK->OUTFILE->f, "\nSIMILARITY RANGE:%4ld%29s%4ld\n",
	  HOMRANGE, "MIN.PERCENT SIMILARITY:", MINPER);
  fprintf(LINK->OUTFILE->f, "SCALE FACTOR:%8.2f%18s%15ld\n",
	  SCALEFAC, "COMPRESSION:", COMPFACT);
}  /*HEADER*/

/* Print a horizontal axis */
Local Void HORAXIS(LEFT, RIGHT, LINK)
long LEFT, RIGHT;
struct LOC_QUICKSEARCH *LINK;
{
  long NUMBER, DELTA, I, NUMPRINTED;

  fprintf(LINK->OUTFILE->f, "\n%7c", ' ');
  /* Write numbers */
  DELTA = COMPFACT * 10;
  NUMPRINTED = (RIGHT - LEFT + 1) / DELTA;
  if (INVERSEX) {
    DELTA = -DELTA;
    NUMBER = LENX - LEFT + 2;
  } else
    NUMBER = LEFT - 1;
  for (I = 1; I <= NUMPRINTED; I++) {
    NUMBER += DELTA;
    fprintf(LINK->OUTFILE->f, "%10ld", NUMBER);
  }
  putc('\n', LINK->OUTFILE->f);

  /* Write sequence if COMPFACT = 1*/
  if (COMPFACT == 1) {
    fprintf(LINK->OUTFILE->f, "%7c", ' ');
    for (I = LEFT; I <= RIGHT; I++)
      putc(NUCHAR[(long)SEQX[I + MAXRANGE]], LINK->OUTFILE->f);
  }
  putc('\n', LINK->OUTFILE->f);
}  /*HORAXIS*/

/* Local variables for MAKETABLE: */
struct LOC_MAKETABLE {
  struct LOC_QUICKSEARCH *LINK;
  long I;
} ;

Local Void ADDNODE(N_, LINK)
NODE **N_;
struct LOC_MAKETABLE *LINK;
{
  NODE *TEMP;

  TEMP = *N_;
  if (FREELIST == NULL)
    *N_ = (NODE *)Malloc(sizeof(NODE));
  else {
    *N_ = FREELIST;
    FREELIST = FREELIST->NEXT;
  }
  (*N_)->NEXT = TEMP;
  (*N_)->POS = LINK->I;
}

/*  Make a table of locations of trinucleotides in SEQX.  Each triplet*/
/*  has a stack of nodes, each of which holds a location in SEQX      */
/*  at which the trinucleotide occurs.                                */
Local Void MAKETABLE(LEFT, RIGHT, LINK)
long LEFT, RIGHT;
struct LOC_QUICKSEARCH *LINK;
{
  struct LOC_MAKETABLE V;
  NUCLEOTIDE N1, N2, N3;

  V.LINK = LINK;
  for (V.I = RIGHT; V.I >= LEFT; V.I--) {
    N1 = SEQX[V.I + MAXRANGE - 1];
    N2 = SEQX[V.I + MAXRANGE];
    N3 = SEQX[V.I + MAXRANGE + 1];
    if ((long)N1 < (long)N && (long)N2 < (long)N && (long)N3 < (long)N)
      ADDNODE(&TRIPLET[(long)N1 - (long)T][(long)N2 - (long)T]
	      [(long)N3 - (long)T], &V);
  }
}  /*MAKETABLE*/

/* At each occurrence of the triplet in SEQX, compare a region  */
/* HOMRANGE bases on either side.  If the match is good enough, */
/* print a character at the corresponding point in the matrix.  */
Local Void SEARCHFOR(N_, LINK)
NODE *N_;
struct LOC_QUICKSEARCH *LINK;
{
  long LX, LY, RX, RY, DISTANCE, SCORE, PERCENT, X;

  while (N_ != NULL) {
    LINK->POSX = N_->POS;
    SCORE = THREEMATCH;
    LX = LINK->POSX - 2;
    RX = LINK->POSX + 2;
    LY = LINK->POSY - 2;
    RY = LINK->POSY + 2;
    DISTANCE = 2;
    while (DISTANCE <= HOMRANGE) {
      if (SEQX[LX + MAXRANGE] == SEQY[LY + MAXRANGE])
	SCORE += SUBSCORE[DISTANCE];
      if (SEQX[RX + MAXRANGE] == SEQY[RY + MAXRANGE])
	SCORE += SUBSCORE[DISTANCE];
      LX--;
      RX++;
      LY--;
      RY++;
      DISTANCE++;
    }
    if (SCORE >= MINSCORE) {
      PERCENT = (long)floor(SCORE * MAXFACTOR + 0.5);
      X = (LINK->POSX - LINK->LEFT) / COMPFACT + 1;
      if (PERCENT > LINK->THISLINE[X-1])
	LINK->THISLINE[X-1] = PERCENT;
    }
    N_ = N_->NEXT;
  }
}  /* SEARCHFOR */

/* Push the linked-list (if there is one) for each TRIPLET to */
/* the top of the FREELIST and reset TRIPLET to nil.          */
Local Void RIDOF(LINK)
struct LOC_QUICKSEARCH *LINK;
{
  NUCLEOTIDE N1, N2, N3;
  NODE *HEAD, *TAIL;

  for (N1 = T; (long)N1 <= (long)G; N1 = (NUCLEOTIDE)((long)N1 + 1)) {
    for (N2 = T; (long)N2 <= (long)G; N2 = (NUCLEOTIDE)((long)N2 + 1)) {
      for (N3 = T; (long)N3 <= (long)G; N3 = (NUCLEOTIDE)((long)N3 + 1)) {
	if (TRIPLET[(long)N1 - (long)T][(long)N2 - (long)T]
	    [(long)N3 - (long)T] != NULL) {
	  HEAD = TRIPLET[(long)N1 - (long)T][(long)N2 - (long)T]
	    [(long)N3 - (long)T];
	  TAIL = HEAD;
	  while (TAIL->NEXT != NULL)
	    TAIL = TAIL->NEXT;
	  TAIL->NEXT = FREELIST;
	  FREELIST = HEAD;
	  TRIPLET[(long)N1 - (long)T][(long)N2 - (long)T]
	    [(long)N3 - (long)T] = NULL;
	}
      }
    }
  }
}  /*RIDOF*/


/* ******************************************************************* */
/*  Find similarities (2*HOMRANGE)+1 long with MINPER or better match. */
/* ******************************************************************* */
Static Void QUICKSEARCH(OUTFILE_)
_TEXT *OUTFILE_;
{
  struct LOC_QUICKSEARCH V;
  long RIGHT, RIGHTLIM, I;
  long J = 1;
  long K_, INDEX, NEXTLINE;
  short ONELINE[MAXLINE], TENLINE[MAXLINE];
  NUCLEOTIDE N1, N2, N3;
  long FORLIM, FORLIM1;

  V.OUTFILE = OUTFILE_;
  /* Initialize LEFT & RIGHT, which define the part of SEQX */
  /* to be compared with SEQY.                              */
  if (INVERSEX) {
    V.LEFT = LENX - STARTX + 1;
    RIGHTLIM = LENX - FINISHX + 1;
  } else {
    V.LEFT = STARTX;
    RIGHTLIM = FINISHX;
  }
  RIGHT = V.LEFT + LINESIZE * COMPFACT - 1;
  if (RIGHT > RIGHTLIM)
    RIGHT = RIGHTLIM;

  FORLIM = LINESIZE;
  /* Initialize line templates */
  for (I = 0; I < FORLIM; I++) {
    TENLINE[I] = 0;   /* period */
    if (J == 10) {
      ONELINE[I] = 0;
      J = 0;
    } else
      ONELINE[I] = -1;
    /* blank */
    J++;
  }

  /* SIMILARITY SEARCH */
  printf("Search begins...\n");
  HEADER(&V);
  while (V.LEFT <= RIGHTLIM) {
    HORAXIS(V.LEFT, RIGHT, &V);
    MAKETABLE(V.LEFT, RIGHT, &V);
    INDEX = 0;
    K_ = STARTY % 10;
/* p2c: d3hom.p, line 935:
 * Note: Using % for possibly-negative arguments [317] */
    if (K_ == 0)
      memcpy(V.THISLINE, TENLINE, MAXLINE * sizeof(short));
    else
      memcpy(V.THISLINE, ONELINE, MAXLINE * sizeof(short));
    FORLIM = FINISHY;
    for (V.POSY = STARTY; V.POSY <= FORLIM; V.POSY++) {
      N1 = SEQY[V.POSY + MAXRANGE - 1];
      N2 = SEQY[V.POSY + MAXRANGE];
      N3 = SEQY[V.POSY + MAXRANGE + 1];
      if ((long)N1 < (long)N && (long)N2 < (long)N && (long)N3 < (long)N)
	SEARCHFOR(TRIPLET[(long)N1 - (long)T][(long)N2 - (long)T]
		  [(long)N3 - (long)T], &V);
      INDEX++;
      if (INDEX == COMPFACT) {
	if (K_ < 9) {
	  fprintf(V.OUTFILE->f, "%6c", ' ');
	  NEXTLINE = 1;
	} else if (K_ == 9) {
	  fprintf(V.OUTFILE->f, "%6c", ' ');
	  NEXTLINE = 10;
	} else {
	  fprintf(V.OUTFILE->f, "%6ld", V.POSY);
	  NEXTLINE = 1;
	  K_ = 0;
	}
	if (COMPFACT == 1)
	  putc(NUCHAR[(long)SEQY[V.POSY + MAXRANGE]], V.OUTFILE->f);
	else
	  putc(' ', V.OUTFILE->f);
	FORLIM1 = LINESIZE;
	for (I = 0; I < FORLIM1; I++)
	  putc(PC[V.THISLINE[I] + 1], V.OUTFILE->f);
	putc('\n', V.OUTFILE->f);
	if (NEXTLINE == 1)
	  memcpy(V.THISLINE, ONELINE, MAXLINE * sizeof(short));
	else
	  memcpy(V.THISLINE, TENLINE, MAXLINE * sizeof(short));
	INDEX = 0;
	K_++;
      }
    }
    RIDOF(&V);
    V.LEFT = RIGHT + 1;
    RIGHT += LINESIZE * COMPFACT;
    if (RIGHT > RIGHTLIM)
      RIGHT = RIGHTLIM;
  }
  fprintf(V.OUTFILE->f, "\n\n");
}  /* QUICKSEARCH */


/* ----------------------------------------------------------  */
/* ----------------- MAIN  PROCEDURE  -----------------------  */
main(argc, argv)
int argc;
Char *argv[];
{
  _TEXT TEMP1;

  /* BEGIN MODULE STARTUP */
  /* Peform operations which must be done at beginning of main
     procedure. */
  /*!!!   TERMIN(input);    Open input for interactive use */
  /*!!!   TERMOUT(output);   "   output "      "        "  */
  PASCAL_MAIN(argc, argv);
  OUTFILE.f = NULL;
  *OUTFILE.name = '\0';
  SFILEY.f = NULL;
  *SFILEY.name = '\0';
  SFILEX.f = NULL;
  *SFILEX.name = '\0';
  printf("%50s\n\n\n", VERSION);
  /* END MODULE STARTUP         VERSION= 'SUNMODS     Version  6/26/01'; */

  /* Read in X-axis sequence */
  NEWSEQ(&SFILEX, &XFN, SEQX, &LENX, &STARTX, &FINISHX, &NAMEX, &CIRCX, 'X');

  /* Read in Y-axis sequence */
  NEWSEQ(&SFILEY, &YFN, SEQY, &LENY, &STARTY, &FINISHY, &NAMEY, &CIRCY, 'Y');

  /* Open output file. */
  printf("Type output filename:\n");
  GETFILE(&OUTFILE, 'O', &OFN);

  INITIALIZE();

  /* Initialize horizontal output line */

  for (J = 1; J <= MAXLINE; J++)
    HLINE.STR[J-1] = '_';
  HLINE.LEN = MAXLINE;
  /* MAIN LOOP */
  do {
    putchar('\n');
    TEMP1.f = stdout;
    *TEMP1.name = '\0';
    WRITELINE(&TEMP1, HLINE, 80L);
    printf("\nD3HOM   %30s\n", "MAIN MENU");
    TEMP1.f = stdout;
    *TEMP1.name = '\0';
    WRITELINE(&TEMP1, HLINE, 80L);
    printf("\nX-axis file:       ");
    TEMP1.f = stdout;
    *TEMP1.name = '\0';
    WRITELINE(&TEMP1, XFN, XFN.LEN);
    printf("\nY-axis file:       ");
    TEMP1.f = stdout;
    *TEMP1.name = '\0';
    WRITELINE(&TEMP1, YFN, YFN.LEN);
    printf("\nOutput file:       ");
    TEMP1.f = stdout;
    *TEMP1.name = '\0';
    WRITELINE(&TEMP1, OFN, OFN.LEN);
    putchar('\n');
    TEMP1.f = stdout;
    *TEMP1.name = '\0';
    WRITELINE(&TEMP1, HLINE, 80L);
    printf("\n%20c1) Read in a new X-axis sequence\n", ' ');
    printf("%20c2) Read in a new Y-axis sequence\n", ' ');
    printf("%20c3) Open a new output file\n", ' ');
    printf("%20c4) Change parameters\n", ' ');
    printf("%20c5) Compare sequences and write output screen\n", ' ');
    printf("%20c6) Compare sequences and write output to file\n", ' ');
    TEMP1.f = stdout;
    *TEMP1.name = '\0';
    WRITELINE(&TEMP1, HLINE, 80L);
    printf("\nType the number of your choice  (0 to quit program)\n");
    GETREAL(&TEMP, 0.0, 6.0);
    scanf("%*[^\n]");
    getchar();
    CHOICE = (long)floor(TEMP + 0.5);

    switch (CHOICE) {

    case 0:
      /* blank case */
      break;

    case 1:
      NEWSEQ(&SFILEX, &XFN, SEQX, &LENX, &STARTX, &FINISHX, &NAMEX, &CIRCX,
	     'X');
      break;

    case 2:
      NEWSEQ(&SFILEY, &YFN, SEQY, &LENY, &STARTY, &FINISHY, &NAMEY, &CIRCY,
	     'Y');
      break;

    case 3:
      if (OUTFILE.f != NULL)
	fclose(OUTFILE.f);
      OUTFILE.f = NULL;
      printf("Type output filename:\n");
      GETFILE(&OUTFILE, 'O', &OFN);
      break;

    /*!!!*/
    case 4:
      PARAMETERS();
      break;

    case 5:
      CALCSCORES();
      TEMP1.f = stdout;
      *TEMP1.name = '\0';
      QUICKSEARCH(&TEMP1);
      printf("Press RETURN to continue");
      scanf("%*[^\n]");
      getchar();
      break;

    case 6:
      CALCSCORES();
      QUICKSEARCH(&OUTFILE);
      break;
    }
  } while (CHOICE != 0);
  /*!!!*/
  if (OUTFILE.f != NULL)
    fclose(OUTFILE.f);
  OUTFILE.f = NULL;
  if (SFILEX.f != NULL)
    fclose(SFILEX.f);
  if (SFILEY.f != NULL)
    fclose(SFILEY.f);
  if (OUTFILE.f != NULL)
    fclose(OUTFILE.f);
  exit(EXIT_SUCCESS);
}  /* D3HOM  */



/* End. */
