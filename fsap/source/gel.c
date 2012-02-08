/* Output from p2c 1.21alpha-07.Dec.93, the Pascal-to-C translator */
/* From input file "gel.p" */


/*$DEBUG-*/
/* ********************************************************  */
/*                                                           */
/*   GEL       Version  5/13/91  Standard Pascal             */
/*             Brian Fristensky                              */
/*             Dept. of Plant Science                        */
/*             University of Manitoba                        */
/*             Winnipeg, MB R3T 2N2  CANADA                  */
/*                                                           */
/*  Copyright (c) l984-1990           by  Brian Fristensky   */
/*                                                           */
/*  !!! indicates feature which may need change              */
/*  *******************************************************  */
/*  Given the lengths and mobilities of a set of molecular   */
/*  weight standards on a gel, this program will estimate    */
/*  the sizes of unknown fragments whose mobilities are      */
/*  known, using the least squares approach found in:        */
/*                                                           */
/*  Schaffer & Sederoff, ANAL.BIOCHEM 115,p113-122 (1981)    */
/*************************************************************/

/*!!!*/

#include <p2c.h>


/*,OUTFILE*/
/*!!! Some Pascals require file parameters in program heading*/

#define MAXFRAGS        50
#define MAXLINE         150
#define MAXWORD         20
/* BEGIN MODULE REALLIMITS */
/*!!!*/

#define MAXREAL         1.7e38
/* END MODULE REALLIMITS         VERSION= 'SUNMODS     Version  9/12/91'; */

#define VERSION         "GEL            Version   5/13/91"


/* BEGIN MODULE TYPE.WORD */
/*   <word>::= <non-blank char>[<non-blank char>] */

typedef struct WORD {
  long LEN;
  Char STR[MAXWORD];
} WORD;

/* END MODULE TYPE.WORD         VERSION= 'SUNMODS     Version  9/12/91'; */

/* BEGIN MODULE TYPE.LINE */
typedef Char CHARARRAY[MAXLINE];

typedef struct LINE {
  CHARARRAY STR;
  long LEN;
} LINE;

/* END MODULE TYPE.LINE         VERSION= 'SUNMODS     Version  9/12/91'; */


Static _TEXT OUTFILE;
Static LINE OFN, HLINE;
Static Char ANSWER;

/* GLOBAL PARAMETERS */
Static double L0, M0, CCAP, SD, SC;
Static double L[MAXFRAGS], M[MAXFRAGS], PREDLEN[MAXFRAGS], D[MAXFRAGS],
	      C[MAXFRAGS];
Static long I, N;


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


/* END MODULE INPLINE         VERSION= 'SUNMODS     Version  9/12/91'; */

/* BEGIN MODULE WRITELINE */
/*  Write a line to a file using L char, left-justified.  */
Static Void WRITELINE(F, W, L)
_TEXT *F;
LINE W;
long L;
{
  long I;

  for (I = 1; I <= L; I++) {
    if (I <= W.LEN)
      putc(W.STR[I-1], F->f);
    else
      putc(' ', F->f);
  }
}  /* WRITELINE */


/* END MODULE WRITELINE         VERSION= 'SUNMODS     Version  9/12/91'; */

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


/* END MODULE GETFILE         VERSION= 'SUNMODS     Version  9/12/91'; */

/* BEGIN MODULE INPWORD */
/* Read a WORD from the terminal. */
Static Void INPWORD(W)
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


/* END MODULE INPWORD         VERSION= 'SUNMODS     Version  9/12/91'; */

/* BEGIN MODULE WRITEWORD */
/*  Write a word to a file using L char, left-justified.  */
Static Void WRITEWORD(F, W, L)
_TEXT *F;
WORD W;
long L;
{
  long I;

  for (I = 1; I <= L; I++) {
    if (I <= W.LEN)
      putc(W.STR[I-1], F->f);
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


/* END MODULE WRITEWORD         VERSION= 'SUNMODS     Version  9/12/91'; */

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


/* END MODULE GETREAL         VERSION= 'SUNMODS     Version  9/12/91'; */

/*********************************************************/
/*  Read in standard fragments from console.             */
/*********************************************************/
Static Void READSTANDARDS()
{
  double LENGTH, MOBILITY;

  fprintf(stdout, "\f");
  printf("Type in the size of each fragment and distances migrated,\n");
  printf("one fragment at a time. Press <RETURN> after each fragment.\n");
  printf("Type 0 0  <RETURN> when finished\n");
  printf(" Length    Mobility\n");
  printf("--------+------------\n");
  GETREAL(&LENGTH, 0.0, MAXREAL);
  GETREAL(&MOBILITY, 0.0, MAXREAL);
  scanf("%*[^\n]");
  getchar();
  N = 0;
  while (LENGTH != 0 && MOBILITY != 0 && N < MAXFRAGS) {
    N++;
    L[N-1] = LENGTH;
    M[N-1] = MOBILITY;
    GETREAL(&LENGTH, 0.0, MAXREAL);
    GETREAL(&MOBILITY, 0.0, MAXREAL);
    scanf("%*[^\n]");
    getchar();
  }
}  /* READSTANDARDS */


/*****************************************************/
/* Display fragment lengths and mobilities of screen */
/*****************************************************/
Static Void DISPLAYFRAGS()
{
  long I;
  _TEXT TEMP;
  long FORLIM;

  fprintf(stdout, "\f");
  TEMP.f = stdout;
  *TEMP.name = '\0';
  WRITELINE(&TEMP, HLINE, 80L);
  printf("\nStandard fragments:      Length  Mobility\n");
  FORLIM = N;
  for (I = 1; I <= FORLIM; I++)
    printf("%20ld)%10.3f%10.3f\n", I, L[I-1], M[I-1]);
  TEMP.f = stdout;
  *TEMP.name = '\0';
  WRITELINE(&TEMP, HLINE, 80L);
  putchar('\n');
}  /* DISPLAYFRAGS */


/****************************************************/
/* Let user edit the standards.                     */
/****************************************************/
Static Void EDIT()
{
  long P;
  double TEMP;

  if (N <= 0) {
    printf("There are  0 fragments.  Type in new standards.\n");
    return;
  }
  do {
    DISPLAYFRAGS();
    printf("Type number of a fragment you wish to change,\n");
    printf("or 0 if all are correct.\n");
    GETREAL(&TEMP, 0.0, (double)N);
    scanf("%*[^\n]");
    getchar();
    P = (long)floor(TEMP + 0.5);

    if (P > 0) {
      do {
	printf("Enter the correct length and mobility:\n");
	GETREAL(&L[P-1], 0.0, MAXREAL);
	GETREAL(&M[P-1], 0.0, MAXREAL);
	scanf("%*[^\n]");
	getchar();
      } while (L[P-1] <= 0 || M[P-1] <= 0);
    }
  } while (P != 0);
}  /* EDIT */


/*************************************/
/* Let user add standard fragments.  */
/*************************************/
Static Void ADDFRAGS()
{
  long I, POSITION;
  double TEMP;

  if (N <= 0) {
    printf("There are  0 fragments.  Type in new standards.\n");
    return;
  }
  do {
    if (N < MAXFRAGS) {
      DISPLAYFRAGS();
      printf("Add a fragment before which number? (0 to quit) \n");
      GETREAL(&TEMP, 0.0, N + 1.0);
      scanf("%*[^\n]");
      getchar();
      POSITION = (long)floor(TEMP + 0.5);
      if (POSITION > 0) {
	/* Push up the fragment stack to make room for a new fragment */
	for (I = N + 1; I > POSITION; I--) {
	  L[I-1] = L[I-2];
	  M[I-1] = M[I-2];
	}
	/* Add a new fragment */
	do {
	  printf("Type new fragment length and mobility\n");
	  GETREAL(&L[POSITION-1], 0.0, MAXREAL);
	  GETREAL(&M[POSITION-1], 0.0, MAXREAL);
	  scanf("%*[^\n]");
	  getchar();
	} while (L[POSITION-1] <= 0 || M[POSITION-1] <= 0);
	N++;
      }
    } else {
      printf("!!! Too many fragments !!!\n");
      POSITION = 0;
    }
  } while (POSITION != 0);
}  /*ADDFRAGS*/


/****************************************/
/* Let user delete standard fragments.  */
/****************************************/
Static Void DELETEFRAGS()
{
  long I, POSITION;
  double TEMP;
  long FORLIM;

  if (N <= 0) {
    printf("!!! There are  0 fragments.  Type in new standards!!! \n");
    return;
  }
  do {
    DISPLAYFRAGS();
    printf("Delete which fragment? (0 to quit)\n");
    GETREAL(&TEMP, 0.0, (double)N);
    scanf("%*[^\n]");
    getchar();
    POSITION = (long)floor(TEMP + 0.5);
    if (POSITION > 0) {
      FORLIM = N;
      for (I = POSITION; I < FORLIM; I++) {
	L[I-1] = L[I];
	M[I-1] = M[I];
      }
      N--;
    }
  } while (POSITION != 0);
}  /* DELETEFRAGS*/


/**********************************************************************/
/* Calculate parameters L0,M0,CCAP,SD,& SC  as described in the ref.  */
/**********************************************************************/
Static Void CALCULATE(L0, M0, CCAP, SD, SC)
double *L0, *M0, *CCAP, *SD, *SC;
{
  double MBAR = 0.0, LBAR = 0.0, MLBAR = 0.0;   /* means of M,L, & M*L */
  double CSSM = 0.0, CSSL = 0.0, CSCPML = 0.0;
      /* corrected sums of squares of M,L, & M*L*/
  double CSPMLL = 0.0, CSPMLM = 0.0;
  double DELTA;
  double SCi = 0.0, SCiS = 0.0, SDi = 0.0, SDiS = 0.0;
      /* Sums & sums of sqrs of Ci & Di */
  double PROD[MAXFRAGS], DWT[MAXFRAGS], DDIST[MAXFRAGS], DPROD[MAXFRAGS];
  long I, FORLIM;
  double TEMP;

  /* Calculate means */
  FORLIM = N;
  for (I = 0; I < FORLIM; I++) {
    MBAR += M[I];
    LBAR += L[I];
    PROD[I] = M[I] * L[I];
    MLBAR += PROD[I];
  }
  MBAR /= N;
  LBAR /= N;
  MLBAR /= N;

  FORLIM = N;
  /* Calculate deviations */
  for (I = 0; I < FORLIM; I++) {
    DWT[I] = L[I] - LBAR;
    DDIST[I] = M[I] - MBAR;
    DPROD[I] = PROD[I] - MLBAR;
    /* Calculate intermediate values */
  }

  FORLIM = N;
  for (I = 0; I < FORLIM; I++) {
    TEMP = DWT[I];
    CSSL += TEMP * TEMP;
    TEMP = DDIST[I];
    CSSM += TEMP * TEMP;
    CSCPML += DWT[I] * DDIST[I];
    CSPMLL += DPROD[I] * DWT[I];
    CSPMLM += DPROD[I] * DDIST[I];
  }

  /* Calculate L0, M0, and CCAP */
  DELTA = CSSM * CSSL - CSCPML * CSCPML;
  *M0 = (CSSM * CSPMLL - CSCPML * CSPMLM) / DELTA;
  *L0 = (CSSL * CSPMLM - CSCPML * CSPMLL) / DELTA;

  /* Standard deviation of Ci */
  *CCAP = 0.0;
  FORLIM = N;
  for (I = 0; I < FORLIM; I++) {
    C[I] = (M[I] - *M0) * (L[I] - *L0);
    SCi += C[I] - C[0];
    TEMP = C[I] - C[0];
    SCiS += TEMP * TEMP;
  }
  *CCAP = SCi / N + C[0];
  *SC = sqrt((SCiS - SCi * SCi / N) / (N - 1));

  /* Calculate predicted fragment sizes and deviations from input sizes */
  FORLIM = N;
  for (I = 0; I < FORLIM; I++) {
    PREDLEN[I] = *CCAP / (M[I] - *M0) + *L0;
    D[I] = L[I] - PREDLEN[I];
    SDi += D[I];
    TEMP = D[I];
    SDiS += TEMP * TEMP;
  }

  /* Standard deviation of fragment sizes */
  *SD = sqrt((SDiS - SDi * SDi / N) / (N - 3));

}  /* CALCULATE */


/******************************************************************/
/* Print a report giving the standards and statistics.            */
/******************************************************************/
Static Void REPORT(OUTFILE)
_TEXT *OUTFILE;
{
  long I;
  LINE TITLE;
  long FORLIM;

  printf("Type title to appear on output (<RETURN> for blank):\n");
  INPLINE(&TITLE);
  WRITELINE(OUTFILE, TITLE, TITLE.LEN);
  if (TITLE.LEN > 0)
    putc('\n', OUTFILE->f);

  fprintf(OUTFILE->f, "%10s%10s%10s%10s%10s%10s\n",
	  "STD LEN", "DIST", "PRED LEN", "DEVIATION", "%DEV", "C[I]");
  FORLIM = N;
  for (I = 0; I < FORLIM; I++)
    fprintf(OUTFILE->f, "%10.2f%10.3f%10.2f%10.3f%10.3f%10.3f\n",
	    L[I], M[I], PREDLEN[I], D[I], D[I] / L[I] * 100, C[I]);
  fprintf(OUTFILE->f, "M0= %10.4f  L0= %10.4f  CCAP= %10.4f\n", M0, L0, CCAP);
  fprintf(OUTFILE->f, "SC= %10.4f  SD= %10.4f\n\n", SC, SD);
}  /* REPORT */


/****************************************************************/
/* Prompt user for migration distances of unknown fragments and */
/*  calculate the sizes of unknowns.                            */
/****************************************************************/
Static Void UNKNOWNS(OUTFILE)
_TEXT *OUTFILE;
{
  double DISTANCE, LENGTH;
  LINE NAME;

  fprintf(OUTFILE->f, "UNKNOWN FRAGMENTS:\n");
  fprintf(OUTFILE->f, "%10s%10s PREDICTED LENGTH\n", "FRAGMENT", "DISTANCE");
  printf("Type an identifier (<=10 letters) for an unknown fragment\n");
  printf("  (<RETURN> to quit)\n");
  INPLINE(&NAME);
  while (NAME.LEN > 0) {
    printf("Distance migrated?\n");
    GETREAL(&DISTANCE, 0.0, MAXREAL);
    scanf("%*[^\n]");
    getchar();
    LENGTH = CCAP / (DISTANCE - M0) + L0;
    WRITELINE(OUTFILE, NAME, 10L);
    fprintf(OUTFILE->f, "%10.2f%17.3f\n", DISTANCE, LENGTH);
    printf("Type an identifier (<=10 letters) for an unknown fragment\n");
    printf("  (<RETURN> to quit)\n");
    INPLINE(&NAME);
  }
  putc('\n', OUTFILE->f);
}  /* UNKNOWNS */


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
  printf("%50s\n\n", VERSION);
  /* END MODULE STARTUP         VERSION= 'SUNMODS     Version  9/12/91'; */

  /* Initialize horizontal output line */

  /* Open output file */
  for (I = 1; I <= MAXLINE; I++)
    HLINE.STR[I-1] = '_';
  HLINE.LEN = MAXLINE;
  printf("Type output filename:\n");
  GETFILE(&OUTFILE, 'O', &OFN);
  N = 0;

  /* MAIN LOOP */
  do {
    putchar('\n');
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, HLINE, 80L);
    printf("\nGEL%40s\n", "MAIN MENU");
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, HLINE, 80L);
    printf("\n%20s", "Output file: ");
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, OFN, OFN.LEN);
    putchar('\n');
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, HLINE, 80L);
    printf("\n\n%20c1) Type in a set of standard fragments\n", ' ');
    printf("%20c2) Edit values of standard fragments\n", ' ');
    printf("%20c3) Add fragments\n", ' ');
    printf("%20c4) Delete fragments\n", ' ');
    printf("%20c5) Calculate sizes of unknowns (output to screen)\n", ' ');
    printf("%20c6) Calculate sizes of unknowns (output to file)\n", ' ');
    printf("%20c7) Open a new output file\n", ' ');
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, HLINE, 80L);
    printf("\nType number of your choice  (0 to quit)\n");
    scanf("%c%*[^\n]", &ANSWER);
    getchar();
    if (ANSWER == '\n')
      ANSWER = ' ';
    if (ANSWER >= '1' && ANSWER <= '7') {
      switch (ANSWER) {

      case '1':
	READSTANDARDS();
	break;

      case '2':
	EDIT();
	break;

      case '3':
	ADDFRAGS();
	break;

      case '4':
	DELETEFRAGS();
	break;

      case '5':
	if (N >= 4) {
	  CALCULATE(&L0, &M0, &CCAP, &SD, &SC);
	  TEMP.f = stdout;
	  *TEMP.name = '\0';
	  REPORT(&TEMP);
	  TEMP.f = stdout;
	  *TEMP.name = '\0';
	  UNKNOWNS(&TEMP);
	} else
	  printf(">>> There must be at least 4 standards. <<<\n");
	break;

      case '6':
	if (N >= 4) {
	  CALCULATE(&L0, &M0, &CCAP, &SD, &SC);
	  REPORT(&OUTFILE);
	  UNKNOWNS(&OUTFILE);
	} else
	  printf(">>> There must be at least 4 standards. <<<\n");
	break;

      case '7':
	if (OUTFILE.f != NULL)
	  fclose(OUTFILE.f);
	OUTFILE.f = NULL;
	printf("Type new output filename:\n");
	GETFILE(&OUTFILE, 'O', &OFN);
	break;
	/*!!!*/
      }
    }
  } while (ANSWER != '0');
  /*!!!*/
  if (OUTFILE.f != NULL)
    fclose(OUTFILE.f);
  OUTFILE.f = NULL;
  if (OUTFILE.f != NULL)
    fclose(OUTFILE.f);
  exit(EXIT_SUCCESS);
}  /* GEL */




/* End. */
