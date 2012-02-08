/* Output from p2c 1.21alpha-07.Dec.93, the Pascal-to-C translator */
/* From input file "p1hom.p" */


/*$DEBUG- */
/* ********************************************************  */
/*                                                           */
/*   P1HOM     Version  5/13/91, Standard Pascal             */
/*             Brian Fristensky                              */
/*             Dept. of Plant Science                        */
/*             University of Manitoba                        */
/*             Winnipeg, MB R3T 2N2  CANADA                  */
/*                                                           */
/* Copyright (c) 1984, 1986, 1987, 1991 by Brian Fristensky. */
/* !!! in comment indicates feature which may need change.   */
/*  *******************************************************  */


#include <p2c.h>


/*,SFILEX,SFILEY,OUTFILE*/
/*!!! Some Pascals require file parameters in program heading */

#define MAXSEQ          9000
/* BEGIN MODULE REALLIMITS */
/*!!!*/

#define MAXREAL         1.7e38
/* END MODULE REALLIMITS         VERSION= 'SUNMODS     Version  9/12/91'; */

#define MAXWORD         25
#define MAXRANGE        30
#define MAXLINE         150

#define VERSION         "P1HOM         Version  5/13/91"


typedef enum {
  GLY, ALA, VALINE, LEU, ILE, MET, PHE, PRO, SER, THR, CYS, ASN, GLN, TYR,
  TRP, ASP, GLU, HIS, LYS, ARG, ASX, GLX, TERM, UNKX, UNKY
} AMINOACID;

typedef AMINOACID PROTEIN[MAXSEQ + MAXRANGE + 1];

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

/* END MODULE TYPE.WORD         VERSION= 'SUNMODS     Version  9/12/91'; */

/* BEGIN MODULE TYPE.LINE */
typedef Char CHARARRAY[MAXLINE];

typedef struct LINE {
  CHARARRAY STR;
  long LEN;
} LINE;

/* END MODULE TYPE.LINE         VERSION= 'SUNMODS     Version  9/12/91'; */


Static _TEXT SFILEX, SFILEY, OUTFILE;
Static LINE XFN, YFN, OFN, HLINE;
Static PROTEIN SEQX, SEQY;
Static long LENX, LENY;
Static WORD NAMEX, NAMEY;
Static NODE *PEPTIDE[(long)ARG - (long)GLY + 1];
Static NODE *FREELIST;
Static Char AACHAR[25];
Static long STARTX, FINISHX, STARTY, FINISHY, HOMRANGE, MINPER, COMPFACT,
	    LINESIZE, MAXSCORE, ONEMATCH, MINSCORE;
Static double SCALEFAC, MAXFACTOR, TEMP;
Static long SUBSCORE[MAXRANGE + 1];
Static Char PC[102];   /* Printing characters for graph */
Static long CHOICE, J;


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

/**************************************************/
/*  WORD   I/O  PROCEDURES                        */
/**************************************************/
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

/* BEGIN MODULE READWORD */
/*  Read a word from a textfile           */
Static Void READWORD(F, W)
_TEXT *F;
WORD *W;
{
  long I;
  Char CH;

  W->LEN = 0;
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
    if (W->LEN < MAXWORD) {
      W->LEN++;
      W->STR[W->LEN - 1] = getc(F->f);
      if (W->STR[W->LEN - 1] == '\n')
	W->STR[W->LEN - 1] = ' ';
    } else {
      CH = getc(F->f);
      if (CH == '\n')
	CH = ' ';
    }
  }
  for (I = W->LEN; I < MAXWORD; I++)
    W->STR[I] = ' ';
}  /* READWORD */


/* END MODULE READWORD         VERSION= 'SUNMODS     Version  9/12/91'; */

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
	PR[J + MAXRANGE] = GLY;
	break;

      case 'A':
	PR[J + MAXRANGE] = ALA;
	break;

      case 'V':
	PR[J + MAXRANGE] = VALINE;
	break;

      case 'L':
	PR[J + MAXRANGE] = LEU;
	break;

      case 'I':
	PR[J + MAXRANGE] = ILE;
	break;

      case 'M':
	PR[J + MAXRANGE] = MET;
	break;

      case 'F':
	PR[J + MAXRANGE] = PHE;
	break;

      case 'P':
	PR[J + MAXRANGE] = PRO;
	break;

      case 'S':
	PR[J + MAXRANGE] = SER;
	break;

      case 'T':
	PR[J + MAXRANGE] = THR;
	break;

      case 'C':
	PR[J + MAXRANGE] = CYS;
	break;

      case 'N':
	PR[J + MAXRANGE] = ASN;
	break;

      case 'Q':
	PR[J + MAXRANGE] = GLN;
	break;

      case 'Y':
	PR[J + MAXRANGE] = TYR;
	break;

      case 'W':
	PR[J + MAXRANGE] = TRP;
	break;

      case 'D':
	PR[J + MAXRANGE] = ASP;
	break;

      case 'E':
	PR[J + MAXRANGE] = GLU;
	break;

      case 'H':
	PR[J + MAXRANGE] = HIS;
	break;

      case 'K':
	PR[J + MAXRANGE] = LYS;
	break;

      case 'R':
	PR[J + MAXRANGE] = ARG;
	break;

      case 'X':
	PR[J + MAXRANGE] = UNKX;
	break;

      case 'B':
	PR[J + MAXRANGE] = ASX;
	break;

      case 'Z':
	PR[J + MAXRANGE] = GLX;
	break;

      case '*':
	PR[J + MAXRANGE] = TERM;
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


/* END MODULE READPRO         VERSION= 'SUNMODS     Version  9/12/91'; */

/* **************************************************** */
/*              INITIALIZATION  PROCEDURES              */
/* **************************************************** */
Static Void INITIALIZE()
{
  AMINOACID A1;

  /* Initialize the PC array which holds the       */
  /* characers which symbolize percent identity.   */
  PC[0] = ' ';
  PC[1] = '.';
  PC[6] = 'w';
  PC[7] = 'v';
  PC[8] = 'v';
  PC[9] = 'u';
  PC[10] = 'u';
  PC[11] = 't';
  PC[12] = 't';
  PC[13] = 's';
  PC[14] = 's';
  PC[15] = 'r';
  PC[16] = 'r';
  PC[17] = 'q';
  PC[18] = 'q';
  PC[19] = 'p';
  PC[20] = 'p';
  PC[21] = 'o';
  PC[22] = 'o';
  PC[23] = 'n';
  PC[24] = 'n';
  PC[25] = 'm';
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

  /*   Initialize AACHAR array which holds the            */
  /*   character values of the AMINOACIDs.                */
  AACHAR[(long)GLY] = 'G';
  AACHAR[(long)ALA] = 'A';
  AACHAR[(long)VALINE] = 'V';
  AACHAR[(long)LEU] = 'L';
  AACHAR[(long)ILE] = 'I';
  AACHAR[(long)MET] = 'M';
  AACHAR[(long)PHE] = 'F';
  AACHAR[(long)PRO] = 'P';
  AACHAR[(long)SER] = 'S';
  AACHAR[(long)THR] = 'T';
  AACHAR[(long)CYS] = 'C';
  AACHAR[(long)ASN] = 'N';
  AACHAR[(long)GLN] = 'Q';
  AACHAR[(long)TYR] = 'Y';
  AACHAR[(long)TRP] = 'W';
  AACHAR[(long)ASP] = 'D';
  AACHAR[(long)GLU] = 'E';
  AACHAR[(long)HIS] = 'H';
  AACHAR[(long)LYS] = 'K';
  AACHAR[(long)ARG] = 'R';
  AACHAR[(long)ASX] = 'B';
  AACHAR[(long)GLX] = 'Z';
  AACHAR[(long)TERM] = '*';
  AACHAR[(long)UNKX] = 'X';
  AACHAR[(long)UNKY] = 'X';

  /*Set default values for parameters */
  HOMRANGE = 10;
  SCALEFAC = 0.90;
  MINPER = 50;
  COMPFACT = 10;
  LINESIZE = 70;

  /* Initialize PEPTIDE table to nil */
  for (A1 = GLY; (long)A1 <= (long)ARG; A1 = (AMINOACID)((long)A1 + 1)) {
    PEPTIDE[(long)A1 - (long)GLY] = NULL;

  }
}  /* INITIALIZE */


Local Void SEQENDS(PR, LEN, UNKNOWN)
AMINOACID *PR;
long LEN;
AMINOACID UNKNOWN;
{
  long J;

  /* Initialize negative and last parts of the sequence*/
  for (J = -MAXRANGE; J <= 0; J++)
    PR[J + MAXRANGE] = UNKNOWN;
  for (J = LEN + 1; J <= LEN + MAXRANGE; J++)
    PR[J + MAXRANGE] = UNKNOWN;
}  /* SEQENDS */


/********************************************************/
/* Open a sequence file and read in sequence.           */
/********************************************************/
Static Void NEWSEQ(SFILE, FN, PR, LEN, START, FINISH, NAME, AXIS)
_TEXT *SFILE;
LINE *FN;
AMINOACID *PR;
long *LEN, *START, *FINISH;
WORD *NAME;
Char AXIS;
{
  long J, FORLIM;

  printf("Enter filename for sequence on %c-axis:\n", AXIS);
  GETFILE(SFILE, 'I', FN);
  NAME->LEN = 0;
  READPRO(SFILE, PR, LEN, NAME);
  switch (AXIS) {

  case 'X':
    SEQENDS(PR, *LEN, UNKX);
    break;

  case 'Y':
    SEQENDS(PR, *LEN, UNKY);
    FORLIM = *LEN;
    for (J = 1; J <= FORLIM; J++) {
      if (PR[J + MAXRANGE] == UNKX)
	PR[J + MAXRANGE] = UNKY;
    }
    break;
  }
  if (NAME->LEN == 0) {
    printf("Type name for %c-axis sequence to appear on output:\n", AXIS);
    INPWORD(NAME);
    scanf("%*[^\n]");
    getchar();
  }
  *START = 1;
  *FINISH = *LEN;
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

/* Display  parameters on screen */
Local Void DISPLAY()
{
  _TEXT TEMP1;

  fprintf(stdout, "\f");
  printf("X-axis sequence: ");
  TEMP1.f = stdout;
  *TEMP1.name = '\0';
  WRITEWORD(&TEMP1, NAMEX, 20L);
  printf("%21s%10ld aa\n", " Length: ", LENX);
  printf("Y-axis sequence: ");
  TEMP1.f = stdout;
  *TEMP1.name = '\0';
  WRITEWORD(&TEMP1, NAMEY, 20L);
  printf("%21s%10ld aa\n", " Length: ", LENY);
  TEMP1.f = stdout;
  *TEMP1.name = '\0';
  WRITELINE(&TEMP1, HLINE, 80L);
  printf("\n%12cParameter   Description/Response                 Value\n",
	 ' ');
  TEMP1.f = stdout;
  *TEMP1.name = '\0';
  WRITELINE(&TEMP1, HLINE, 80L);
  printf("\n%12c 1)STARTX   first amino acid position in SEQX   %6ld\n",
	 ' ', STARTX);
  printf("%12c 2)FINISHX  last  amino acid position in SEQX   %6ld\n",
	 ' ', FINISHX);
  printf("%12c 3)STARTY   first amino acid position in SEQY   %6ld\n",
	 ' ', STARTY);
  printf("%12c 4)FINISHY  last  amino acid position in SEQY   %6ld\n",
	 ' ', FINISHY);
  printf("%12c 5)HOMRANGE dist.from central a.a. in a match.  %6ld\n",
	 ' ', HOMRANGE);
  printf("%12c 6)SCALEFAC scale factor for exponential curve  %6.2f\n",
	 ' ', SCALEFAC);
  printf("%12c 7)MINPER   minimum percent similarity printed  %6ld\n",
	 ' ', MINPER);
  printf("%12c 8)COMPFACT graph compression factor            %6ld\n",
	 ' ', COMPFACT);
  printf("%12c 9)LINESIZE width of output line (ex. 70,120)   %6ld\n",
	 ' ', LINESIZE);
  TEMP1.f = stdout;
  *TEMP1.name = '\0';
  WRITELINE(&TEMP1, HLINE, 80L);
  putchar('\n');
}  /* DISPLAY */


/* **************************************************** */
/* Prompt user for parameters used by program.          */
/* **************************************************** */
Static Void PARAMETERS()
{
  long RESPONSE;
  double TEMP;

  /* Prompt user for new parameter values */
  do {
    DISPLAY();
    printf("Type number of parameter you wish to change (0 to continue)\n");
    GETREAL(&TEMP, 0.0, 9.0);
    scanf("%*[^\n]");
    getchar();
    RESPONSE = (long)floor(TEMP + 0.5);
    if ((unsigned long)RESPONSE < 32 && ((1L << RESPONSE) & 0x3fe) != 0) {
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
	GETNUMBER(&MINPER, "MINPER    ", 5L, 100L);
	break;

      case 8:
	GETNUMBER(&COMPFACT, "COMPFACT  ", 1L, 100L);
	break;

      case 9:
	GETNUMBER(&LINESIZE, "LINESIZE  ", 40L, (long)MAXLINE);
	break;
      }
    }
  } while (RESPONSE != 0);   /* PARAMETERS */
}


/****************************************************************/
/* Calculate a table of scores for matches at each distance from*/
/* the central   peptide in a local homology.                   */
/****************************************************************/
Static Void CALCSCORES()
{
  long I;
  long V = 100;
  double S;
  long FORLIM;

  SUBSCORE[0] = V;
  MAXSCORE = V;
  S = SCALEFAC;
  FORLIM = HOMRANGE;
  for (I = 1; I <= FORLIM; I++) {
    SUBSCORE[I] = (long)floor(V * S + 0.5);
    MAXSCORE += SUBSCORE[I] * 2;
    S *= SCALEFAC;
  }
  ONEMATCH = SUBSCORE[0];
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
      putc(AACHAR[(long)SEQX[I + MAXRANGE]], LINK->OUTFILE->f);
  }
  putc('\n', LINK->OUTFILE->f);
}  /*HORAXIS*/

/* Local variables for MAKETABLE: */
struct LOC_MAKETABLE {
  struct LOC_QUICKSEARCH *LINK;
  long I;
} ;

Local Void ADDNODE(N, LINK)
NODE **N;
struct LOC_MAKETABLE *LINK;
{
  NODE *TEMP;

  TEMP = *N;
  if (FREELIST == NULL)
    *N = (NODE *)Malloc(sizeof(NODE));
  else {
    *N = FREELIST;
    FREELIST = FREELIST->NEXT;
  }
  (*N)->NEXT = TEMP;
  (*N)->POS = LINK->I;
}

/*  Make a table of locations of   peptides in SEQX.  Each   peptide  */
/*  has a stack of nodes, each of which holds a location in SEQX      */
/*  at which the   peptide occurs.                                    */
Local Void MAKETABLE(LEFT, RIGHT, LINK)
long LEFT, RIGHT;
struct LOC_QUICKSEARCH *LINK;
{
  struct LOC_MAKETABLE V;
  AMINOACID A1;

  V.LINK = LINK;
  for (V.I = RIGHT; V.I >= LEFT; V.I--) {
    A1 = SEQX[V.I + MAXRANGE];
    if ((long)A1 <= (long)ARG)
      ADDNODE(&PEPTIDE[(long)A1 - (long)GLY], &V);
  }
}  /*MAKETABLE*/

/* At each occurrence of the   peptide in SEQX, compare a region  */
/* HOMRANGE bases on either side.   If the match is good enough,  */
/* print a character at the corresponding point in the matrix.    */
Local Void SEARCHFOR(N, LINK)
NODE *N;
struct LOC_QUICKSEARCH *LINK;
{
  long LX, LY, RX, RY, DISTANCE, SCORE, PERCENT, X;

  while (N != NULL) {
    LINK->POSX = N->POS;
    SCORE = ONEMATCH;
    LX = LINK->POSX - 1;
    RX = LINK->POSX + 1;
    LY = LINK->POSY - 1;
    RY = LINK->POSY + 1;
    DISTANCE = 1;
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
    N = N->NEXT;
  }
}  /* SEARCHFOR */

/* Push the linked-list (if there is one) for each PEPTIDE to */
/* the top of the FREELIST and reset PEPTIDE to nil.          */
Local Void RIDOF(LINK)
struct LOC_QUICKSEARCH *LINK;
{
  AMINOACID A1;
  NODE *HEAD, *TAIL;

  for (A1 = GLY; (long)A1 <= (long)ARG; A1 = (AMINOACID)((long)A1 + 1)) {
    if (PEPTIDE[(long)A1 - (long)GLY] != NULL) {
      HEAD = PEPTIDE[(long)A1 - (long)GLY];
      TAIL = HEAD;
      while (TAIL->NEXT != NULL)
	TAIL = TAIL->NEXT;
      TAIL->NEXT = FREELIST;
      FREELIST = HEAD;
      PEPTIDE[(long)A1 - (long)GLY] = NULL;
    }
  }
}  /*RIDOF*/


/* ***************************************************************** */
/*  Find homologies (2*HOMRANGE)+1 long with MINPER or better match. */
/* ***************************************************************** */
Static Void QUICKSEARCH(OUTFILE_)
_TEXT *OUTFILE_;
{
  struct LOC_QUICKSEARCH V;
  long RIGHT, RIGHTLIM, I;
  long J = 1;
  long K, INDEX, NEXTLINE;
  short ONELINE[MAXLINE], TENLINE[MAXLINE];
  AMINOACID A1;
  long FORLIM, FORLIM1;

  V.OUTFILE = OUTFILE_;
  /* Initialize LEFT & RIGHT, which define the part of SEQX */
  /* to be compared with SEQY.                              */
  V.LEFT = STARTX;
  RIGHTLIM = FINISHX;
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
    K = STARTY % 10;
/* p2c: p1hom.p, line 790:
 * Note: Using % for possibly-negative arguments [317] */
    if (K == 0)
      memcpy(V.THISLINE, TENLINE, MAXLINE * sizeof(short));
    else
      memcpy(V.THISLINE, ONELINE, MAXLINE * sizeof(short));
    FORLIM = FINISHY;
    for (V.POSY = STARTY; V.POSY <= FORLIM; V.POSY++) {
      A1 = SEQY[V.POSY + MAXRANGE];
      if ((long)A1 <= (long)ARG)
	SEARCHFOR(PEPTIDE[(long)A1 - (long)GLY], &V);
      INDEX++;
      if (INDEX == COMPFACT) {
	if (K < 9) {
	  fprintf(V.OUTFILE->f, "%6c", ' ');
	  NEXTLINE = 1;
	} else if (K == 9) {
	  fprintf(V.OUTFILE->f, "%6c", ' ');
	  NEXTLINE = 10;
	} else {
	  fprintf(V.OUTFILE->f, "%6ld", V.POSY);
	  NEXTLINE = 1;
	  K = 0;
	}
	if (COMPFACT == 1)
	  putc(AACHAR[(long)SEQY[V.POSY + MAXRANGE]], V.OUTFILE->f);
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
	K++;
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
  /* END MODULE STARTUP         VERSION= 'SUNMODS     Version  9/12/91'; */

  /* Read in X-axis sequence */
  NEWSEQ(&SFILEX, &XFN, SEQX, &LENX, &STARTX, &FINISHX, &NAMEX, 'X');

  /* Read in Y-axis sequence */
  NEWSEQ(&SFILEY, &YFN, SEQY, &LENY, &STARTY, &FINISHY, &NAMEY, 'Y');

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
    printf("\nP1HOM   %30s\n", "MAIN MENU");
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
    printf("%20c5) Compare sequences and write output to screen\n", ' ');
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
      NEWSEQ(&SFILEX, &XFN, SEQX, &LENX, &STARTX, &FINISHX, &NAMEX, 'X');
      break;

    case 2:
      NEWSEQ(&SFILEY, &YFN, SEQY, &LENY, &STARTY, &FINISHY, &NAMEY, 'Y');
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
}  /* P1HOM  */



/* End. */
