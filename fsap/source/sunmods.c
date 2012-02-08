/* Output from p2c 1.21alpha-07.Dec.93, the Pascal-to-C translator */
/* From input file "sunmods.p" */


/* ********************************************************  */
/*                                                           */
/* FRISTENSKY SEQUENCE ANALYSIS PACKAGE SUBROUTINE MODULES   */
/*             VERSION  6/26/01  SUN    Pascal               */
/*             Brian Fristensky                              */
/*             Dept. of Plant Science                        */
/*             University of Manitoba                        */
/*             Winnipeg, MB R3T 2N2 CANADA                   */
/*                                                           */
/* Copyright (c) 1986-2001               by Brian Fristensky */
/* !!! in comment indicates feature which may need change.   */
/*  *******************************************************  */

/* --------  REVISION HISTORY --------------------------------
26 Jun 01 - Readseq updated to accommodate changes to the
GenBank LOCUS line format, which should take effect with
GenBank 126.0. The new format allows names up to 18 char. long,
and sequences up to 9,999,999 bases. One side effect is that
the topology has moved from columns 43-52 to columns 56-63.
------------------------------------------------------------*/


#include <p2c.h>


/*,INFILE*/
/*!!! Some Pascals require file parameters in program heading */

#define MAXSEQ          750000L
/* BEGIN MODULE VERSION */

#define VERSION         "SUNMODS     Version  6/26/01"
/* END MODULE VERSION */

/* BEGIN MODULE REALLIMITS */
/*!!!*/

#define MAXREAL         1.7e38
/* END MODULE REALLIMITS */

/* BEGIN MODULE INTLIMITS */
/*!!!  MAXINT =  2147483647; */
/*!!!  MININT = -2147483647; */
/* END MODULE INTLIMITS */

/* BEGIN MODULE STARTARGNUM */

#define STARTARGNUM     1
    /* SUN Pascal: ARG(1) is 1st command line argument*/
/*STARTARGNUM=2;       ATT Pascal: ARG(2) is 1st command line argument*/
/* END MODULE STARTARGNUM */

#define MAXWORD         25
#define MAXLINE         70
#define MAXRANGE        0


typedef enum {
  CANT, A, C, R, D, V, M, K, B, H, Y, G, T, W, S, N, Z, WONT
} NUCLEOTIDE;
typedef enum {
  GLY, ALA, VALINE, LEU, ILE, MET, PHE, PRO, SER, THR, CYS, ASN, GLN, TYR,
  TRP, ASP, GLU, HIS, LYS, ARG, ASX, GLX, TERM, UNKX, UNKY
} AMINOACID;
typedef NUCLEOTIDE SEQUENCE[MAXSEQ];
typedef AMINOACID PROTEIN[3000];
typedef long CHSET[4];

typedef Char LETTERS[10];
typedef long VECTOR[MAXLINE];

/* BEGIN MODULE TYPE.WORD */
/*   <word>::= <non-blank char>[<non-blank char>] */

typedef struct WORD {
  long LEN;
  Char STR[MAXWORD];
} WORD;

/* END MODULE TYPE.WORD         VERSION= 'CSAPMODS     Version  4/13/88'; */

/* BEGIN MODULE TYPE.LINE */
typedef Char CHARARRAY[MAXLINE];

typedef struct LINE {
  CHARARRAY STR;
  long LEN;
} LINE;

/* END MODULE TYPE.LINE         VERSION= 'SUNMODS     Version  3/21/90'; */

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

/* END MODULE TYPE.FRAG */


Static _TEXT INFILE;
Static Char CH;
Static WORD WO;
Static LINE LI;
Static double RE;
Static SEQUENCE SEQ;
Static long SEQLEN, INT;
Static boolean CIRCULAR;
Static FRAGMENT *FREEFRAG;
Static FRAGMENT *ORDER[10];


/************ DUMMY LINES FOR MODULE STRIPPING **********/
/* BEGIN MODULE DUMMY */
/* BEGIN MODULE REALLIMITS */
/* END MODULE REALLIMITS */

/* BEGIN MODULE INTLIMITS */
/* END MODULE INTLIMITS */

/* BEGIN MODULE TYPE.WORD */
/* END MODULE TYPE.WORD */

/* BEGIN MODULE TYPE.LINE */
/* END MODULE TYPE.LINE */

/* BEGIN MODULE TYPE.FRAG */
/* END MODULE TYPE.FRAG */

/* BEGIN MODULE INPLINE */
/* END MODULE INPLINE */

/* BEGIN MODULE READLINE */
/* END MODULE READLINE */

/* BEGIN MODULE WRITELINE */
/* END MODULE WRITELINE */

/* BEGIN MODULE FILEARGS */
/* END MODULE FILEARGS */

/* BEGIN MODULE GETFILE */
/* END MODULE GETFILE */

/* BEGIN MODULE INPWORD */
/* END MODULE INPWORD */

/* BEGIN MODULE READWORD */
/* END MODULE READWORD */

/* BEGIN MODULE WRITEWORD */
/* END MODULE WRITEWORD */

/* BEGIN MODULE TOUPPER */
/* END MODULE TOUPPER */

/* BEGIN MODULE TOLOWER */
/* END MODULE TOLOWER */

/* BEGIN MODULE GETCHAR */
/* END MODULE GETCHAR */

/* BEGIN MODULE GETREAL */
/* END MODULE GETREAL */

/* BEGIN MODULE GETINTEGER */
/* END MODULE GETINTEGER */

/* BEGIN MODULE NUC */
/* END MODULE NUC */

/* BEGIN MODULE AA */
/* END MODULE AA */

/* BEGIN MODULE SAMEWORD */
/* END MODULE SAMEWORD */

/* BEGIN MODULE READSEQ */
/* END MODULE READSEQ */

/* BEGIN MODULE READPRO */
/* END MODULE READPRO */

/* BEGIN MODULE LINKED */
/* END MODULE LINKED */

/* BEGIN MODULE SORT */
/* END MODULE SORT */

/* BEGIN MODULE MV */
/* END MODULE MV */

/* BEGIN MODULE STARTUP */
/* END MODULE STARTUP */
/* END MODULE DUMMY */

/********************* ACTUAL PROCEDURES ***************/
/**********************************************************/
/*  LINE I/O Procedures.                                  */
/**********************************************************/
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


/* END MODULE INPLINE */

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


/* END MODULE READLINE */

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


/* END MODULE WRITELINE         VERSION= 'CSAPMODS     Version  4/13/88'; */

/* BEGIN MODULE FILEARGS */
/* This procedure overcomes one of the stupidest aspects of UNIX Pascal,
   namely the fact that filenames in the program statement are supposed to
   be actual UNIX filenames!  To overcome this, the 2-argument version of
   reset and rewrite must be used with string variables.  This module
   need only contain the reset and rewrite statements in any normal
   implementation of Pascal. */
Static Void FILEARGS(F, FTYPE, ARGNUM)
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


/* END MODULE FILEARGS */


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


/* END MODULE GETFILE         VERSION= 'SUNMODS     Version  3/21/90'; */

/**********************************************************/
/*  WORD I/O Procedures.                                  */
/**********************************************************/
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


/* END MODULE INPWORD */

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


/* END MODULE READWORD         VERSION= 'CSAPMODS     Version  4/13/88'; */

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


/* END MODULE WRITEWORD         VERSION= 'CSAPMODS     Version  4/13/88'; */

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


/* END MODULE TOUPPER */

/* BEGIN MODULE TOLOWER */
/* Change a character from upper to lowercase */
Static Char TOLOWER(CH)
Char CH;
{
  Char Result;

  if (!isupper(CH))
    return CH;
  switch (CH) {

  case 'A':
    Result = 'a';
    break;

  case 'B':
    Result = 'b';
    break;

  case 'C':
    Result = 'c';
    break;

  case 'D':
    Result = 'd';
    break;

  case 'E':
    Result = 'e';
    break;

  case 'F':
    Result = 'f';
    break;

  case 'G':
    Result = 'g';
    break;

  case 'H':
    Result = 'h';
    break;

  case 'I':
    Result = 'i';
    break;

  case 'J':
    Result = 'j';
    break;

  case 'K':
    Result = 'k';
    break;

  case 'L':
    Result = 'l';
    break;

  case 'M':
    Result = 'm';
    break;

  case 'N':
    Result = 'n';
    break;

  case 'O':
    Result = 'o';
    break;

  case 'P':
    Result = 'p';
    break;

  case 'Q':
    Result = 'q';
    break;

  case 'R':
    Result = 'r';
    break;

  case 'S':
    Result = 's';
    break;

  case 'T':
    Result = 't';
    break;

  case 'U':
    Result = 'u';
    break;

  case 'V':
    Result = 'v';
    break;

  case 'W':
    Result = 'w';
    break;

  case 'X':
    Result = 'x';
    break;

  case 'Y':
    Result = 'y';
    break;

  case 'Z':
    Result = 'z';
    break;
  }
  return Result;
}  /* TOLOWER */


/* END MODULE TOLOWER */

/* BEGIN MODULE GETCHAR */
/* Read a character from the console and check */
/*  for correct response.                      */
Static Void GETCHAR(CH, PNAME, ALLOWED)
Char *CH;
Char *PNAME;
long *ALLOWED;
{
  /* BEGIN MODULE TOUPPER */
  /* END MODULE TOUPPER */
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


/* END MODULE GETCHAR */

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


/* END MODULE GETREAL */

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


/* END MODULE GETINTEGER */

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


/* END MODULE NUC         VERSION= 'CSAPMODS     Version  4/13/88'; */

/* BEGIN MODULE AA */
/*****************************************************************/
/*  Convert a character to the appropriate amino acid symbol.    */
/*****************************************************************/
Static AMINOACID AA(CH)
Char CH;
{
  AMINOACID Result;

  switch (CH) {

  case 'A':
    Result = ALA;
    break;

  case 'C':
    Result = CYS;
    break;

  case 'D':
    Result = ASP;
    break;

  case 'E':
    Result = GLU;
    break;

  case 'F':
    Result = PHE;
    break;

  case 'G':
    Result = GLY;
    break;

  case 'H':
    Result = HIS;
    break;

  case 'I':
    Result = ILE;
    break;

  case 'K':
    Result = LYS;
    break;

  case 'L':
    Result = LEU;
    break;

  case 'M':
    Result = MET;
    break;

  case 'N':
    Result = ASN;
    break;

  case 'P':
    Result = PRO;
    break;

  case 'Q':
    Result = GLN;
    break;

  case 'R':
    Result = ARG;
    break;

  case 'S':
    Result = SER;
    break;

  case 'T':
    Result = THR;
    break;

  case 'V':
    Result = VALINE;
    break;

  case 'W':
    Result = TRP;
    break;

  case 'Y':
    Result = TYR;
    break;

  case 'B':
    Result = ASX;
    break;

  case 'Z':
    Result = GLX;
    break;

  case 'X':
    Result = UNKX;
    break;

  case '*':
    Result = TERM;
    break;
  }
  return Result;
}  /* AA */


/* END MODULE AA */

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


/* END MODULE SAMEWORD  9/12/91 */

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
/* p2c: sunmods.p, line 833: 
 * Note: Line breaker spent 0.0 seconds, 5000 tries on line 1488 [251] */
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


/* END MODULE READSEQ */

/* BEGIN MODULE READPRO */
/*  ******************************************* */
/*  Read a protein    sequence from SFILE       */
/*   and store it in PR.                        */
/*  ******************************************* */
Static Void READPRO(SFILE, PR, LEN, NAME)
_TEXT *SFILE;
AMINOACID *PR;
long *LEN;
WORD *NAME;
{
  Char CH, ANSWER, FILETYPE;
  long J;
  boolean TOOLONG = false;
  long FORLIM;

  /* Prompt user for sequence file type */
  printf("The following file formats can be read:\n");
  printf("  F:free format    N:NBRF\n");
  do {
    printf("Type letter of format (F|N)\n");
    scanf("%c%*[^\n]", &ANSWER);
    getchar();
    if (ANSWER == '\n')
      ANSWER = ' ';
  } while (ANSWER != 'n' && ANSWER != 'N' && ANSWER != 'f' && ANSWER != 'F');
  switch (ANSWER) {

  case 'F':
  case 'f':
    FILETYPE = 'F';
    break;

  case 'N':
  case 'n':
    FILETYPE = 'N';
    break;
  }
  printf("Reading input file...\n");

  /* NBRF: read sequence name and title line */
  /* 'Official' NBRF files have a name line followed by a title line.
      In the name line, four data characters preceed the name itself.
      These are deleted by the program.  Short NBRF files have a title
      following the name, on the same line. The name is not preceeded
      by data characters. */
  NAME->LEN = 0;
  if (FILETYPE == 'N') {
    if (P_peek(SFILE->f) == '>') {
      CH = getc(SFILE->f);
      if (CH == '\n')
	CH = ' ';
    }
    READWORD(SFILE, NAME);
    if (NAME->STR[2] == ';') {
      FORLIM = NAME->LEN;
      /* with NAME */
      for (J = 4; J <= FORLIM; J++)
	NAME->STR[J-4] = NAME->STR[J-1];
      NAME->LEN -= 3;
      fscanf(SFILE->f, "%*[^\n]");
      getc(SFILE->f);
    }
    fscanf(SFILE->f, "%*[^\n]");
    getc(SFILE->f);
  }  /* NBRF */

  /* Read in the sequence */
  J = 0;
  while (!BUFEOF(SFILE->f)) {
    while (!P_eoln(SFILE->f)) {
      CH = getc(SFILE->f);
      if (CH == '\n')
	CH = ' ';
      if (CH != 'X' && CH != '*' && CH != 'Z' && CH != 'B' && CH != 'R' &&
	  CH != 'K' && CH != 'H' && CH != 'E' && CH != 'D' && CH != 'W' &&
	  CH != 'Y' && CH != 'Q' && CH != 'N' && CH != 'C' && CH != 'T' &&
	  CH != 'S' && CH != 'P' && CH != 'F' && CH != 'M' && CH != 'I' &&
	  CH != 'L' && CH != 'V' && CH != 'A' && CH != 'G') {
	if (CH == ';') {   /* comment in input file */
	  fscanf(SFILE->f, "%*[^\n]");
	  getc(SFILE->f);
	}
	continue;
      }
      if (J >= MAXSEQ - MAXRANGE) {
	TOOLONG = true;
	continue;
      }
      J++;
      switch (CH) {

      case 'G':
	PR[J-1] = GLY;
	break;

      case 'A':
	PR[J-1] = ALA;
	break;

      case 'V':
	PR[J-1] = VALINE;
	break;

      case 'L':
	PR[J-1] = LEU;
	break;

      case 'I':
	PR[J-1] = ILE;
	break;

      case 'M':
	PR[J-1] = MET;
	break;

      case 'F':
	PR[J-1] = PHE;
	break;

      case 'P':
	PR[J-1] = PRO;
	break;

      case 'S':
	PR[J-1] = SER;
	break;

      case 'T':
	PR[J-1] = THR;
	break;

      case 'C':
	PR[J-1] = CYS;
	break;

      case 'N':
	PR[J-1] = ASN;
	break;

      case 'Q':
	PR[J-1] = GLN;
	break;

      case 'Y':
	PR[J-1] = TYR;
	break;

      case 'W':
	PR[J-1] = TRP;
	break;

      case 'D':
	PR[J-1] = ASP;
	break;

      case 'E':
	PR[J-1] = GLU;
	break;

      case 'H':
	PR[J-1] = HIS;
	break;

      case 'K':
	PR[J-1] = LYS;
	break;

      case 'R':
	PR[J-1] = ARG;
	break;

      case 'X':
	PR[J-1] = UNKX;
	break;

      case 'B':
	PR[J-1] = ASX;
	break;

      case 'Z':
	PR[J-1] = GLX;
	break;

      case '*':
	PR[J-1] = TERM;
	if (FILETYPE == 'N')
	  goto _L86;   /*ignore rest of file*/
	break;
      }
    }
    fscanf(SFILE->f, "%*[^\n]");
    getc(SFILE->f);
  }

  /* !!! */
_L86:
  /*CLOSE(SFILE); */
  /* branch destination for end of NBRF seq. */
  if (TOOLONG)
    printf(">>> WARNING! Sequence exceeds MAXSEQ-MAXRANGE. Sequence truncated.\n");
  *LEN = J;
}  /* READPRO */


/* END MODULE READPRO */

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


Local Void SWAP(FIRST, SECOND)
FRAGMENT **FIRST, **SECOND;
{
  FRAGMENT *TEMP;

  TEMP = *FIRST;
  *FIRST = *SECOND;
  *SECOND = TEMP;
}


/* END MODULE LINKED */

/* BEGIN MODULE SORT */
/*  Invariant:  The array elements > TOP are sorted.*/
/*    TOP >= unsorted elements.                     */
Static Void BUBBLESORT(TOP, BOTTOM)
long TOP, BOTTOM;
{
  long SMALLEST, NEXT;

  while (TOP >= BOTTOM) {
    /*bubble smallest unsorted number to the top of sorted list*/
    SMALLEST = BOTTOM;
    NEXT = BOTTOM + 1;
    while (NEXT <= TOP) {
      if (ORDER[SMALLEST-1]->SIZE < ORDER[NEXT-1]->SIZE)
	SWAP(&ORDER[SMALLEST-1], &ORDER[NEXT-1]);
      SMALLEST = NEXT;
      NEXT = SMALLEST + 1;
    }
    TOP--;
  }
}  /* BUBBLESORT */


/* END MODULE SORT */

/* BEGIN MODULE MV */
/*****************************************************************/
/*  Calculate the mean and variance of first N elements          */
/*  in array X.  The type VECTOR must be declared in a higher    */
/*  scope as   array[1..<positive constant>] of <real|integer>   */
/*****************************************************************/
Static Void MV(X, N_, MEAN, VARIANCE)
long *X;
long N_;
double *MEAN, *VARIANCE;
{
  double SX = 0.0, SSX = 0.0;
  long I, TEMP;

  for (I = 0; I < N_; I++) {
    SX += X[I];
    TEMP = X[I];
    SSX += TEMP * TEMP;
  }
  *MEAN = SX / N_;
  *VARIANCE = (SSX - SX * SX / N_) / N_ - 1;
}  /* MV */


/* END MODULE MV */

/* ----------------------------------------------------------  */
/* ----------------- MAIN  PROCEDURE  -----------------------  */
main(argc, argv)
int argc;
Char *argv[];
{
  long SET[4];
  _TEXT TEMP;

  /* BEGIN MODULE STARTUP */
  /* Peform operations which must be done at beginning of main
     procedure. */
  /*!!!   TERMIN(input);    Open input for interactive use */
  /*!!!   TERMOUT(output);   "   output "      "        "  */
  PASCAL_MAIN(argc, argv);
  INFILE.f = NULL;
  *INFILE.name = '\0';
  printf("%50s\n\n", VERSION);
  /* END MODULE STARTUP */

  printf("This program tests the module procedures in CSAPMOD.\n\n");
  printf("TEST OF INTERACTIVE INPUT\n");
  CH = ' ';
  printf("When prompted, type a letter between A and Z\n");
  GETCHAR(&CH, "CHARACTER:", P_addsetr(P_expset(SET, 0L), 'A', 'Z'));
  printf("The character is: %c\n", CH);
  printf("Type a single word and press RETURN\n");
  INPWORD(&WO);
  scanf("%*[^\n]");
  getchar();
  printf("The word you typed was: ");
  TEMP.f = stdout;
  *TEMP.name = '\0';
  WRITEWORD(&TEMP, WO, WO.LEN);
  printf("\nType a line of text and press RETURN\n");
  INPLINE(&LI);
  printf("The line you typed was:\n");
  TEMP.f = stdout;
  *TEMP.name = '\0';
  WRITELINE(&TEMP, LI, LI.LEN);
  printf("\nTEST OF NUMERICAL INPUT\n");
  printf("Type an integer between %12ld and %12ldand press RETURN\n",
	 LONG_MIN, LONG_MAX);
  GETINTEGER(&INT, LONG_MIN, LONG_MAX);
  scanf("%*[^\n]");
  getchar();
  printf("The integer you typed was: %12ld\n", INT);
  printf("Type a real number between % .5E and % .5E\n", -MAXREAL, MAXREAL);
  printf("and press RETURN\n");
  GETREAL(&RE, -MAXREAL, MAXREAL);
  scanf("%*[^\n]");
  getchar();
  printf("The real number you typed was: % .5E\n\n", RE);
  printf("TEST OF FILE I/O\n");
  printf("Enter sequence filename:\n");
  GETFILE(&INFILE, 'I', &LI);
  READSEQ(&INFILE, SEQ, &SEQLEN, &WO, &CIRCULAR);
  printf("The sequence in\n");
  TEMP.f = stdout;
  *TEMP.name = '\0';
  WRITELINE(&TEMP, LI, LI.LEN);
  printf("\nis %12ld bases long.\n", SEQLEN);
  if (INFILE.f != NULL)
    fclose(INFILE.f);
  exit(EXIT_SUCCESS);
}  /* CSAPMOD */




/* End. */
