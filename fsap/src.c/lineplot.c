/* Output from p2c 1.21alpha-07.Dec.93, the Pascal-to-C translator */
/* From input file "lineplot.p" */


/*$DEBUG- */
/* ************************************************************** */
/*                                                                */
/*   LINEPLOT    Version  4/27/92, Standard Pascal                */
/*               Brian Fristensky                                 */
/*               Dept. of Plant Science                           */
/*               University of Manitoba                           */
/*               Winnipeg, MB R3T 2N2  CANADA                     */
/*                                                                */
/*  Copyright (c) 1984-1990 by Brian Fristensky                   */
/* !!! in comment indicates feature which may need change         */
/* ************************************************************** */


#include <p2c.h>


/*,INFILE,OUTFILE*/
/*!!! Some Pascals require file parameters in program heading */

#define MAXABCISSA      120
#define MAXORDINATE     100
/* BEGIN MODULE REALLIMITS */
/*!!!*/

#define MAXREAL         1.7e38
/* END MODULE REALLIMITS         VERSION= 'SUNMODS     Version  6/26/01'; */

#define MAXPOINTS       12000   /* for DOS use 6000 */
#define MAXWORD         25
#define MAXLINE         150

#define VERSION         "LINEPLOT       Version 4/27/92"


typedef struct ROW {
  Char STR[MAXABCISSA + 1];
  long LEN;
} ROW;

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
typedef double POINTS[MAXPOINTS];


Static ROW GRAPH[MAXORDINATE + 1];   /*GRAPH CREATED SO FAR*/
Static _TEXT INFILE, OUTFILE;
Static LINE IFN, OFN;
Static Char SCALES, AXES, GRAPHTYPE;
Static ROW HTIT, VTIT;
Static long I, J, ABCISSA, ORDINATE;   /* Width and height of graph (char) */
Static double MINHSCAL, MAXHSCAL;
    /* Max & min values for horizontal scale */
Static double MINVSCAL, MAXVSCAL;   /* Max & min values for vertical scale */
Static double FV;   /*VERTICAL INTERPOLATION FACTOR*/
Static double FH;   /*HORIZONTAL INTERPOLATION FACTOR*/
Static double HAX, VAX;
Static double HSCALE[11], VSCALE[11];


/* BEGIN MODULE INPLINE */
/*  Read a line from the console.    */
Static Void INPLINE(W)
LINE *W;
{
  long I;
  Char CH;
  Char BACKSPACE = '\b';

  /*!!!   if eoln(input) then readln; */
  /*!!!   get(input); */
  W->LEN = 0;
  while (!P_eoln(stdin)) {
    CH = getchar();
    if (CH == '\n')
      CH = ' ';
    /*!!!       if not eoln then*/
    if (CH == BACKSPACE && W->LEN > 0) {
      putchar(BACKSPACE);
      W->LEN--;
    } else if (W->LEN < MAXLINE) {
      W->LEN++;
      W->STR[W->LEN - 1] = CH;
    }
  }
  scanf("%*[^\n]");
  getchar();
  for (I = W->LEN; I < MAXLINE; I++)
    W->STR[I] = ' ';
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


Local Void READTITLE(INFILE, TITLE)
_TEXT *INFILE;
ROW *TITLE;
{
  TITLE->LEN = 0;
  while (!P_eoln(INFILE->f) && TITLE->LEN < MAXABCISSA) {
    TITLE->LEN++;
    TITLE->STR[TITLE->LEN] = getc(INFILE->f);
    if (TITLE->STR[TITLE->LEN] == '\n')
      TITLE->STR[TITLE->LEN] = ' ';
  }
  fscanf(INFILE->f, "%*[^\n]");
  getc(INFILE->f);
}  /* READTITLE */


/* END MODULE GETFILE         VERSION= 'SUNMODS     Version  6/26/01'; */

/* ************************************************************* */
/*  Initialize the graph to blanks.                              */
/* ************************************************************* */
Static Void INITIALIZE()
{
  /*  Read a title from the input file.              */
  printf("Initializing...\n");
  ABCISSA = 50;
  ORDINATE = 40;
  /* Read parameters used in graph from INFILE. */
  READTITLE(&INFILE, &HTIT);
  READTITLE(&INFILE, &VTIT);
  fscanf(INFILE.f, "%lg%lg%lg%lg%lg%lg%*[^\n]", &MINHSCAL, &MAXHSCAL,
	 &MINVSCAL, &MAXVSCAL, &HAX, &VAX);
  getc(INFILE.f);

  SCALES = 'Y';
  AXES = 'Y';
}  /* INITIALIZE */


/***********************************************************/
/*  Skip to the next non-blank character or end-of-file    */
/***********************************************************/
Static Void SKIPBLANKS(F)
_TEXT *F;
{
  while ((P_peek(F->f) == ' ') & (!BUFEOF(F->f)))
    getc(F->f);
}


/*******************************************************/
/* Check to see whether V is within scale coordinates. */
/*******************************************************/
Static boolean INRANGE(V, LOWERLIM, UPPERLIM)
double V, LOWERLIM, UPPERLIM;
{
  if (V < LOWERLIM || V > UPPERLIM)
    return false;
  else
    return true;
}  /* INRANGE */


typedef Char LETTERS[10];
typedef long CHSET[4];


/*!!!*/

#define HIGHEXP         38
/*!!!*/
#define LOWEXP          (-38)


typedef char EXPONENT;



/* BEGIN MODULE INPWORD */
/* Read a WORD from the terminal. */
Local Void INPWORD(W)
WORD *W;
{
  long I;
  Char CH;
  Char BACKSPACE = '\b';

  W->LEN = 0;
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
    if (CH == BACKSPACE && W->LEN > 0) {
      putchar(BACKSPACE);
      W->LEN--;
    } else if (W->LEN < MAXWORD) {
      W->LEN++;
      W->STR[W->LEN - 1] = CH;
    }
  }
  /*readln;*/
  for (I = W->LEN; I < MAXWORD; I++)
    W->STR[I] = ' ';
}  /* INPWORD */


Local double TEN(E)
EXPONENT E;
{
  /* = 10**E */
  long I = 0;
  double T = 1.0;

  do {
    if (E & 1) {
      switch (I) {

      case 0:
	T *= 1.0e1;
	break;

      case 1:
	T *= 1.0e2;
	break;

      case 2:
	T *= 1.0e4;
	break;

      case 3:
	T *= 1.0e8;
	break;

      case 4:
	T *= 1.0e16;
	break;

      case 5:
	T *= 1.0e32;
	break;
	/*!!!   6: T:= T * 1.0E64; Max. exponent is 38
	        7: T:= T * 1.0E128;
	        8: T:= T * 1.0E256; */
      }
    }
    E /= 2;
    I++;
  } while (E != 0);
  return T;
}  /* TEN */

/* END MODULE INPWORD         VERSION= 'SUNMODS     Version  6/26/01'; */

/* BEGIN MODULE GETREAL */
/* The procedure GETREAL has been adapted from the procedure 'rdr', from
   Jensen,K., and Wirth,N. (1974) Pascal User Manual and Report, 2nd
   Ed., Springer-Verlag, pp122-123.  The scope of real numbers considered
   legal by GETREAL includes all possible real numbers as defined on p111
   with two additional allowances: (i) numbers may begin with a decimal
   point (ii) 'E' or 'e' may be used to indicate exponentiation. */

Local Void GETREAL(VAL, LOW, HIGH)
double *VAL, LOW, HIGH;
{
  Char CH;
  double LIMIT = MAXREAL / 10;
  double Y;
  long ORDZERO = '0';
  long I, E;
  boolean LEGAL, INRANGE;
  WORD NUMWORD;
  boolean NEG1, NEG2;
  double A, S;
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
      A = 0.0;
      E = 0;
      LCOUNT = 1;
      while (P_inset(CH, DIGITS)) {
	if (A < LIMIT)
	  A = 10 * A + CH - ORDZERO;
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
	  if (A < LIMIT) {
	    A = 10 * A + CH - ORDZERO;
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
	S = 0.0;
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
	S = CH - ORDZERO;
	I++;
	CH = NUMWORD.STR[I-1];
	while (P_inset(CH, DIGITS)) {
	  if (S < LIMIT)
	    S = 10 * S + CH - ORDZERO;
	  I++;
	  CH = NUMWORD.STR[I-1];
	}
	if (NEG2)
	  E -= (long)floor(S + 0.5);
	else
	  E += (long)floor(S + 0.5);
      }

      /* Check for errors  */
      if (I <= NUMWORD.LEN)   /*illegal char*/
	goto _L777;
      if (E - RCOUNT < LOWEXP || E + LCOUNT > HIGHEXP)
	goto _L776;

      /* Calculate final value */
      Y = A;
      if (NEG1)
	Y = -Y;
      if (E < 0)
	*VAL = Y / TEN((int)(-E));
      else if (E > 0)
	*VAL = Y * TEN((int)E);
      else
	*VAL = Y;
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


/*  Read an integer parameter from the console and check */
/*    that it is in range.                               */
Local Void GETNUMBER(P, PNAME, LOW, HIGH)
double *P;
Char *PNAME;
double LOW, HIGH;
{
  printf("\nType new value for %.10s  (CURRENT VALUE: % .5E)\n", PNAME, *P);
  GETREAL(P, LOW, HIGH);
  scanf("%*[^\n]");
  getchar();
  putchar('\n');
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
  fprintf(stdout, "\f");
  printf("Parameter   Description/Response                 Value\n");
  printf("------------------------------------------------------\n");
  printf(" 1)ABCISSA  width  (in characters) of graph     %10ld\n", ABCISSA);
  printf(" 2)ORDINATE height (in characters) of graph     %10ld\n", ORDINATE);
  printf(" 3)MINHSCAL minimum horizontal scale coordinate % .3E\n", MINHSCAL);
  printf(" 4)MAXHSCAL maximum horizontal scale coordinate % .3E\n", MAXHSCAL);
  printf(" 5)MINVSCAL minimum vertical   scale coordinate % .3E\n", MINVSCAL);
  printf(" 6)MAXVSCAL maximum vertical   scale coordinate % .3E\n", MAXVSCAL);
  printf(" 7)HAX      Y-coordinate of horizontal axis     % .3E\n", HAX);
  printf(" 8)VAX      X-coordinate of vertical   axis     % .3E\n", VAX);
  printf(" 9)SCALES   Print scales               (Y/N)    %10c\n", SCALES);
  printf("10)AXES     Print axes                 (Y/N)    %10c\n\n", AXES);
  printf("Type number of parameter you wish to change\n");
  printf("    (0 to continue)\n");
}  /* DISPLAY */


/* **************************************************** */
/* Prompt user for parameters used by program.          */
/* **************************************************** */
Static Void PARAMETERS()
{
  double TEMP;
  long RESPONSE;
  long SET[4];

  /* Set default values for parameters */
  SCALES = 'Y';
  AXES = 'Y';

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
	TEMP = ABCISSA;
	GETNUMBER(&TEMP, "ABCISSA   ", 50.0, (double)MAXABCISSA);
	ABCISSA = (long)floor(TEMP + 0.5);
	break;

      case 2:
	TEMP = ORDINATE;
	GETNUMBER(&TEMP, "ORDINATE  ", 40.0, (double)MAXORDINATE);
	ORDINATE = (long)floor(TEMP + 0.5);
	break;

      case 3:
	GETNUMBER(&MINHSCAL, "MINHSCAL  ", -MAXREAL, MAXHSCAL);
	break;

      case 4:
	GETNUMBER(&MAXHSCAL, "MAXHSCAL  ", MINHSCAL, MAXREAL);
	break;

      case 5:
	GETNUMBER(&MINVSCAL, "MINVSCAL  ", -MAXREAL, MAXVSCAL);
	break;

      case 6:
	GETNUMBER(&MAXVSCAL, "MAXVSCAL  ", MINVSCAL, MAXREAL);
	break;

      case 7:
	GETNUMBER(&HAX, "HAX       ", -MAXREAL, MAXREAL);
	break;

      case 8:
	GETNUMBER(&VAX, "VAX       ", -MAXREAL, MAXREAL);
	break;

      case 9:
	P_addset(P_expset(SET, 0L), 'Y');
	GETCHAR(&SCALES, "SCALES    ", P_addset(SET, 'N'));
	break;

      case 10:
	P_addset(P_expset(SET, 0L), 'Y');
	GETCHAR(&AXES, "AXES      ", P_addset(SET, 'N'));
	break;
      }
    }
  } while (RESPONSE != 0);   /* PARAMETERS */
}


/* *********************************************************** */
/*  Set the horizontal scale.                                  */
/* *********************************************************** */
Static Void SCALEH(MINHSCAL, MAXHSCAL)
double MINHSCAL, MAXHSCAL;
{
  long J;
  double INCREMENT;

  FH = ABCISSA / (MAXHSCAL - MINHSCAL);
  INCREMENT = (MAXHSCAL - MINHSCAL) / 10;
  HSCALE[0] = MINHSCAL;
  for (J = 1; J <= 10; J++)
    HSCALE[J] = HSCALE[J-1] + INCREMENT;
}  /* SCALEH */


/*  ***********************************************************  */
/*  Set the vertical scale.                                      */
/*  ***********************************************************  */
Static Void SCALEV(MINVSCAL, MAXVSCAL)
double MINVSCAL, MAXVSCAL;
{
  double INCREMENT;
  long J;

  FV = ORDINATE / (MAXVSCAL - MINVSCAL);
  INCREMENT = (MAXVSCAL - MINVSCAL) / 10;
  VSCALE[0] = MINVSCAL;
  for (J = 1; J <= 10; J++)
    VSCALE[J] = VSCALE[J-1] + INCREMENT;
}  /* SCALEV */


/* ************************************************************ */
/* Set the horizontal axis.                                     */
/* ************************************************************ */
Static Void AXISH(HAX)
double HAX;
{
  long I = 0;
  long J;
  long K = 0;
  long INTERVAL;
  ROW *WITH;

  if (!INRANGE(HAX, MINVSCAL, MAXVSCAL))
    return;
  J = (long)floor(ORDINATE + FV * (HAX - MAXVSCAL) + 0.5);
  WITH = &GRAPH[J];
  WITH->STR[0] = '-';
  INTERVAL = ABCISSA / 10;
  while (I < ABCISSA) {
    I++;
    K++;
    if (K == INTERVAL) {
      WITH->STR[I] = '+';
      K = 0;
    } else
      WITH->STR[I] = '-';
  }
  WITH->LEN = ABCISSA;
}  /* AXISH */


/* ************************************************************** */
/* Set the vertical axis.                                         */
/* ************************************************************** */
Static Void AXISV(VAX)
double VAX;
{
  long I, J, FORLIM;

  if (!INRANGE(VAX, MINHSCAL, MAXHSCAL))
    return;
  I = (long)floor(ABCISSA + FH * (VAX - MAXHSCAL) + 0.5);
  FORLIM = ORDINATE;
  for (J = 1; J <= FORLIM; J++) {
    GRAPH[J].STR[I] = '|';
    if (I > GRAPH[J].LEN)
      GRAPH[J].LEN = I;
  }
}  /* AXISV */


#define CIRCUMFERENCE   6.283185307   /* 2 * Pi */


/* Local variables for DRAW: */
struct LOC_DRAW {
  long N;
} ;

/* ************************************************************** */
/* Calculate slope M and constant B by method of least squares.   */
/* ************************************************************** */
Local Void LEASTSQ(X, Y, M, B, N, LINK)
double *X, *Y;
double *M, *B;
long N;
struct LOC_DRAW *LINK;
{
  long I;
  double SX = 0.0, SY = 0.0, SXY = 0.0, SXSQ = 0.0;
  double TEMP;

  /* Calculate sums */
  for (I = 0; I < N; I++) {
    SX += X[I];
    SY += Y[I];
    SXY += X[I] * Y[I];
    TEMP = X[I];
    SXSQ += TEMP * TEMP;
  }

  *M = (SXY - SX * SY / N) / (SXSQ - SX * SX / N);
  *B = (SY - *M * SX) / N;
}  /* LEASTSQ */

/* ********************************************************** */
/* Find the limits of the coordinates plotted.                */
/*    (i.e. the minimum and maximum x and y values)           */
/* ********************************************************** */
Local Void COORD(X, Y, X1, Y1, X2, Y2, M, B, N, LINK)
double *X, *Y;
double *X1, *Y1, *X2, *Y2, M, B;
long N;
struct LOC_DRAW *LINK;
{
  long I;

  /* Find X1 and X2 */
  *X1 = X[0];
  *X2 = X[1];
  for (I = 1; I < N; I++) {
    if (*X1 > X[I])
      *X1 = X[I];
    if (*X2 < X[I])
      *X2 = X[I];
  }
  /* Calculate Y1 and Y2. */
  *Y1 = M * *X1 + B;
  *Y2 = M * *X2 + B;
}  /* COORD */

Local Void SWAP(F, S)
double *F, *S;
{
  double T;

  T = *F;
  *F = *S;
  *S = T;
}  /* SWAP */

/***********************************************************/
/* Draw a line between points X1,Y1 and X2,Y2              */
/* Only the part of the line within scale limits is drawn. */
/***********************************************************/
Local Void MAKELINE(X1, Y1, X2, Y2, C, LINK)
double X1, Y1, X2, Y2;
Char C;
struct LOC_DRAW *LINK;
{
  long W, Z, W1, Z1;   /*FIRST POINT*/
  long W2, Z2;   /*SECOND POINT*/
  double M;   /*SLOPE*/


  /* Set X and Y coordinates to fall within scale boundaries */

  /* VERTICAL LINE - a special case */
  if (X1 == X2) {
    if (Y1 > Y2) {
      SWAP(&X1, &X2);
      SWAP(&Y1, &Y2);
    }
    if (Y1 < MINVSCAL)
      Y1 = MINVSCAL;
    if (Y2 > MAXVSCAL)
      Y2 = MAXVSCAL;
  }

  else {
    /* Make sure X1 < X2 */
    if (X1 > X2) {
      SWAP(&X1, &X2);
      SWAP(&Y1, &Y2);
    }
    M = (Y2 - Y1) / (X2 - X1);   /*slope*/
  }
  /* X1 <> X2 */

  /* Make sure X values are in range */
  if (X1 > MAXHSCAL || X2 < MINHSCAL)
    return;
  if (X1 < MINHSCAL) {
    Y1 = M * (MINHSCAL - X1) + Y1;
    X1 = MINHSCAL;
  }
  if (X2 > MAXHSCAL) {
    Y2 += M * (MAXHSCAL - X2);
    X2 = MAXHSCAL;
  }

  /* Make sure Y values are in range */
  if (Y1 > Y2) {
    SWAP(&X1, &X2);
    SWAP(&Y1, &Y2);
  }
  if (Y1 > MAXVSCAL || Y2 < MINVSCAL)
    return;
  if (Y1 < MINVSCAL) {
    X1 = (MINVSCAL - Y1) / M + X1;
    Y1 = MINVSCAL;
  }
  if (Y2 > MAXVSCAL) {
    X2 += (MAXVSCAL - Y2) / M;
    Y2 = MAXVSCAL;
  }
  /* Make sure X1 < X2 */
  if (X1 > X2) {
    SWAP(&X1, &X2);
    SWAP(&Y1, &Y2);
  }

  /*CONVERT SCALE COORDINATES TO GRAPH INDICES*/
  W1 = (long)floor(ABCISSA + FH * (X1 - MAXHSCAL) + 0.5);
  W2 = (long)floor(ABCISSA + FH * (X2 - MAXHSCAL) + 0.5);
  Z1 = (long)floor(ORDINATE + FV * (Y1 - MAXVSCAL) + 0.5);
  Z2 = (long)floor(ORDINATE + FV * (Y2 - MAXVSCAL) + 0.5);

  /*DETERMINE POINTS NEEDED TO APPROXIMATE A LINE*/
  /*BETWEEN (W1,Z1) AND (W2,Z2)*/
  if (W1 == W2) {
    if (Z1 == Z2) {
      GRAPH[Z1].STR[W1] = C;   /*SINGLE POINT*/
      if (W1 > GRAPH[Z1].LEN)
	GRAPH[Z1].LEN = W1;
      return;
    }
    for (Z = Z1; Z <= Z2; Z++) {  /*VERTICAL LINE*/
      GRAPH[Z].STR[W1] = C;
      if (W1 > GRAPH[Z].LEN)
	GRAPH[Z].LEN = W1;
    }
    return;
  }
  M = (double)(Z2 - Z1) / (W2 - W1);
  if (fabs(M) <= 0.5) {
    for (W = W1; W <= W2; W++) {
      Z = (long)floor(M * (W - W1) + Z1 + 0.5);
      if (INRANGE((double)Z, 0.0, (double)ORDINATE)) {
	GRAPH[Z].STR[W] = C;
	if (W > GRAPH[Z].LEN)
	  GRAPH[Z].LEN = W;
      }
    }
    return;
  }
  if (Z1 > Z2) {
    W = W1;
    W1 = W2;
    W2 = W;
    Z = Z1;
    Z1 = Z2;
    Z2 = Z;
  }
  for (Z = Z1; Z <= Z2; Z++) {
    W = (long)floor((Z - Z1) / M + W1 + 0.5);
    if (INRANGE((double)W, 0.0, (double)ABCISSA)) {
      GRAPH[Z].STR[W] = C;
      if (W > GRAPH[Z].LEN)
	GRAPH[Z].LEN = W;
    }
  }
}  /* MAKELINE */

/********************************************/
/* Store the coordinates of a circle in X,Y */
/********************************************/
Local Void CIRCLE(X, Y, RADIUS, X0, Y0, LINK)
double *X, *Y;
double RADIUS, X0, Y0;
struct LOC_DRAW *LINK;
{
  double THETA = 0.0, DELTATHETA = CIRCUMFERENCE / 360;

  LINK->N = 0;
  while (THETA < CIRCUMFERENCE) {
    LINK->N++;
    X[LINK->N-1] = X0 + RADIUS * cos(THETA);
    Y[LINK->N-1] = Y0 + RADIUS * sin(THETA);
    THETA += DELTATHETA;
  }
}  /* CIRCLE */

#undef CIRCUMFERENCE

/* ************************************************************** */
/* Put X,Y coordinates into graph.                                */
/* ************************************************************** */
Local Void PLOT(X, Y, C, N, LINK)
double *X, *Y;
Char C;
long N;
struct LOC_DRAW *LINK;
{
  long I, J, P;
  ROW *WITH;

  for (P = 0; P < N; P++) {
    if (INRANGE(X[P], MINHSCAL, MAXHSCAL) & INRANGE(Y[P], MINVSCAL, MAXVSCAL)) {
      I = (long)floor(ABCISSA + FH * (X[P] - MAXHSCAL) + 0.5);
      J = (long)floor(ORDINATE + FV * (Y[P] - MAXVSCAL) + 0.5);
      WITH = &GRAPH[J];
      WITH->STR[I] = C;
      if (I > WITH->LEN)
	WITH->LEN = I;
    }
  }
}  /* PLOT */


/* **************************************************************** */
/*  Draw the graph as specified by GRAPHTYPE from input.            */
/*    L- draw a straight line through data points using the         */
/*       method of lease squares                                    */
/*    C- connect each adjacent pair of points with a line           */
/*    O- circle with specified RADIUS centered on X0,Y0             */
/*    P- plot the points only                                       */
/* **************************************************************** */
Static Void DRAW()
{
  struct LOC_DRAW V;
  POINTS X, Y;
  double X1, Y1, X2, Y2, M, B, RADIUS;
  long I;
  Char C;
  long FORLIM;


  fscanf(INFILE.f, "%ld", &V.N);   /* Number of points*/
  SKIPBLANKS(&INFILE);
  C = getc(INFILE.f);
  if (C == '\n')
    C = ' ';
  FORLIM = V.N;
  /* Read in Cartesian coordinates. */
  for (I = 0; I < FORLIM; I++)
    fscanf(INFILE.f, "%lg%lg", &X[I], &Y[I]);

  /* L- LINE */
  if (GRAPHTYPE == 'L') {
    if (V.N > 2) {  /* Best line through N points */
      LEASTSQ(X, Y, &M, &B, V.N, &V);
      COORD(X, Y, &X1, &Y1, &X2, &Y2, M, B, V.N, &V);
      MAKELINE(X1, Y1, X2, Y2, '.', &V);
    } else if (V.N == 2) {   /*Two points*/
      MAKELINE(X[0], Y[0], X[1], Y[1], '.', &V);
      /* 0 or 1 points */
    }
  }

  /* C- CURVE */
  if (GRAPHTYPE == 'C') {
    I = 1;
    while (I < V.N) {
      MAKELINE(X[I-1], Y[I-1], X[I], Y[I], '.', &V);
      I++;
    }
  }

  /* O- CIRCLE */
  if (GRAPHTYPE == 'O') {
    fscanf(INFILE.f, "%lg%*[^\n]", &RADIUS);
    getc(INFILE.f);
    CIRCLE(X, Y, RADIUS, X[0], Y[0], &V);
  }

  /* P- POINTS (always plotted) */

  PLOT(X, Y, C, V.N, &V);
}  /* DRAW */


/*  Write a line to OUTFILE.*/
Local Void WRITESTRING(STRING)
ROW STRING;
{
  long I;

  for (I = 0; I <= STRING.LEN; I++)
    putc(STRING.STR[I], OUTFILE.f);
  putc('\n', OUTFILE.f);
}  /*WRITESTRING */


/* ********************************* */
/* Print the graph a line at a time. */
/* ********************************* */
Static Void PRINTGRAPH()
{
  long LINENUM;
  long CHARCOUNT = 1, K = 1, SCALECOUNT = 10;
  long FORLIM;

  for (LINENUM = ORDINATE; LINENUM >= 0; LINENUM--) {
    if (CHARCOUNT <= VTIT.LEN) {
      putc(VTIT.STR[CHARCOUNT], OUTFILE.f);
      CHARCOUNT++;
    } else
      putc(' ', OUTFILE.f);
    if (K == 1 && SCALES == 'Y') {
      fprintf(OUTFILE.f, "% .4E", VSCALE[SCALECOUNT]);
      SCALECOUNT--;
    } else
      fprintf(OUTFILE.f, "%11c", ' ');
    K++;
    if (K > (long)floor(ORDINATE / 10.0 + 0.5))
      K = 1;
    WRITESTRING(GRAPH[LINENUM]);
  }
  if (SCALES == 'Y') {
    fprintf(OUTFILE.f, "%13.3f", HSCALE[0]);
    for (SCALECOUNT = 1; SCALECOUNT <= 10; SCALECOUNT++) {
      if (SCALECOUNT & 1) {
	fprintf(OUTFILE.f, "%*c", (int)(ABCISSA / 100 * 10 + 2), ' ');
	/*!!!   else write(OUTFILE,HSCALE[SCALECOUNT]:8:3)   !!! UCSD Pascal only*/
	/* !!!*/
      } else {
	fprintf(OUTFILE.f, "% .1E", HSCALE[SCALECOUNT]);
	/* standard Pascal */
      }
    }
  }
  putc('\n', OUTFILE.f);
  FORLIM = HTIT.LEN;
  for (CHARCOUNT = 1; CHARCOUNT <= FORLIM; CHARCOUNT++)
    putc(HTIT.STR[CHARCOUNT], OUTFILE.f);
  putc('\n', OUTFILE.f);
  /*!!!*/
  if (OUTFILE.f != NULL)
    fclose(OUTFILE.f);
  OUTFILE.f = NULL;
}  /* PRINTGRAPH */


/*------------------------------------------------------------------*/
/*----------------------- MAIN PROCEDURE ---------------------------*/
main(argc, argv)
int argc;
Char *argv[];
{
  ROW *WITH;
  long FORLIM;

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

  /* Read in parameters and set graph accordingly */
  printf("Enter input data filename:\n");
  GETFILE(&INFILE, 'I', &IFN);
  printf("Enter output filename:\n");
  GETFILE(&OUTFILE, 'O', &OFN);
  INITIALIZE();
  PARAMETERS();

  /* Initialize graph */
  WITH = GRAPH;
  FORLIM = ABCISSA;
  for (J = 0; J <= FORLIM; J++)
    WITH->STR[J] = ' ';
  WITH->LEN = 0;
  FORLIM = ORDINATE;
  for (I = 1; I <= FORLIM; I++)
    GRAPH[I] = GRAPH[0];
  SCALEH(MINHSCAL, MAXHSCAL);
  SCALEV(MINVSCAL, MAXVSCAL);
  if (AXES == 'Y') {
    AXISH(HAX);
    AXISV(VAX);
  }

  /*Process a graph*/
  SKIPBLANKS(&INFILE);
  while (!BUFEOF(INFILE.f)) {
    GRAPHTYPE = getc(INFILE.f);
    if (GRAPHTYPE == '\n')
      GRAPHTYPE = ' ';
    if (GRAPHTYPE == 'P' || GRAPHTYPE == 'O' || GRAPHTYPE == 'L' ||
	GRAPHTYPE == 'C')
      DRAW();
    SKIPBLANKS(&INFILE);
  }
  /*!!!*/
  if (INFILE.f != NULL)
    fclose(INFILE.f);
  INFILE.f = NULL;
  PRINTGRAPH();
  if (INFILE.f != NULL)
    fclose(INFILE.f);
  if (OUTFILE.f != NULL)
    fclose(OUTFILE.f);
  exit(EXIT_SUCCESS);
}  /* LINEPLOT */



/* End. */
