/* Output from p2c 1.21alpha-07.Dec.93, the Pascal-to-C translator */
/* From input file "reform.p" */


/***********************************************************/
/*                                                         */
/*  REFORM    VERSION   4/ 1/2000  Standard Pascal           */
/*            Brian Fristensky                             */
/*            Dept. of Plant Science                       */
/*            University of Manitoba                       */
/*            Winnipeg, MB R3T 2N2  CANADA                 */
/*                                                         */
/*  SYNOPSIS                                               */
/*  reform [-gpcnm] [-ln] [-sn] file                       */
/*                           or                            */
/*  ralign file window wordsize mm in del |reform [options]*/
/*                                                         */
/*  DESCRIPTION                                            */
/*  Reformats multiple alignment output                    */
/*                                                         */
/*      g    Gaps are to be represented by dashes (-).     */
/*      p    Bases which agree with the consensus are      */
/*           represented by periods (.).                   */
/*      c    Positions at which all sequences agree are    */
/*           capitalized in the consensus.                 */
/*      n    Sequence data is nucleic acid. Defult protein.*/
/*      fx   Specify input file format, where x is:        */
/*           r:RALIGN  p:PEARSON  i:INTELLIGENETICS (MASE) */
/*      m    each input sequence may run over many lines.  */
/*           Sequences begin with >name, after Pearson.    */
/*           Obsolete, equivalent to -fp                   */
/*      ln   The output linelength is set to n.            */
/*           Default is 70.                                */
/*      sn   numbering starts with n (0 default).          */
/*                                                         */
/*    file   Sequence file as described in ralign docu-    */
/*           mentation.  reform needs to re-read the       */
/*           sequence file read by ralign to get the       */
/*           names of the sequences, which ralign ignores. */
/*                                                         */
/*  Copyright (c) l988-2000  by Brian Fristensky           */
/*  !!! in comment indicates feature which may need change */
/***********************************************************/

#include <p2c.h>


/*,SEQFILE*/
/*!!!  Some Pascals require file parameters in program heading */

#define MAXLINE         80   /* Maximum length of input line */
#define MAXWORD         25   /* Maximum length of word eg. seq. name */
#define MAXSEQ          250   /* Maximum number of sequences in alignment */
#define MAXLEN          25000   /* Maximum length of sequences in alignment */
/* BEGIN MODULE STARTARGNUM */
#define STARTARGNUM     1
    /* SUN Pascal: ARG(1) is 1st command line argument*/


/*STARTARGNUM=2;       ATT Pascal: ARG(2) is 1st command line argument*/
/* END MODULE STARTARGNUM         VERSION= 'SUNMODS     Version  8/ 9/94'; */

typedef enum {
  RALIGN, IG, PEARSON
} ALIGNTYPE;   /* alignment type */

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


/* Names, sequences, and related parameters */
Static WORD NAMES[MAXSEQ + 1];
Static Char ALIGNMENT[MAXSEQ + 1][MAXLEN];
Static long SEQLEN[MAXSEQ + 1];
Static long NUMSEQ;

/* Variables associated with command line options */
Static _TEXT SEQFILE;   /* sequence file, also read by ralign */
Static Char ARGUMENT[132];   /* command line argument */
Static boolean FILEARGUMENT;   /*=false if argument preceeded by '-' */
Static boolean CAPS, PERIODS;
Static ALIGNTYPE FORMAT;
Static Char GAPS, UNKNOWN;
Static long I, ARGNUM, TOTALARGS, STARTNUM, SIGN, LINELENGTH;
Static LINE DUMMY;


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


/* END MODULE FILEARGS         VERSION= 'SUNMODS     Version  8/ 9/94'; */

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


/* END MODULE WRITELINE         VERSION= 'SUNMODS     Version  8/ 9/94'; */

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


/* END MODULE READWORD         VERSION= 'SUNMODS     Version  8/ 9/94'; */

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


/* END MODULE WRITEWORD         VERSION= 'SUNMODS     Version  8/ 9/94'; */

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


/* END MODULE TOLOWER         VERSION= 'SUNMODS     Version  8/ 9/94'; */

/***************************************************************/
/* Read aligned sequences into arrays.                         */
/***************************************************************/
Static Void READALIGNMENT(INFILE, OUTFILE)
_TEXT *INFILE, *OUTFILE;
{
  LINE CLINE;
  long LONGEST = 0;   /* length of longest sequence */
  long I;
  boolean ENDOFSEQUENCE;
  Char CH;

  if (FORMAT == RALIGN) {  /*ie. RALIGN output */
    /* Copy consensus word lines of RALIGN output to OUTFILE. */
    if (!BUFEOF(INFILE->f)) {
      fscanf(INFILE->f, "%*[^\n]");
      getc(INFILE->f);
    }
    /*ignore 1st line of pralign
                                             output */
    while (isdigit(P_peek(INFILE->f))) {
      READLINE(INFILE, &CLINE);
      WRITELINE(OUTFILE, CLINE, CLINE.LEN);
      putc('\n', OUTFILE->f);
    }
    putc('\n', OUTFILE->f);   /* ignore consensus line */
    fscanf(INFILE->f, "%*[^\n]");
    getc(INFILE->f);
  }  /* RALIGN */

  /* Read in each sequence. */
  NUMSEQ = 0;
  for (I = 0; I <= MAXSEQ; I++)
    SEQLEN[I] = 0;

  while (!BUFEOF(INFILE->f)) {
    /* IG,PEARSON: move past comments to name */
    if (((1L << ((long)FORMAT)) &
	 ((1L << ((long)IG)) | (1L << ((long)PEARSON)))) != 0) {
      while (P_peek(INFILE->f) == ';') {   /* ignore comment lines */
	fscanf(INFILE->f, "%*[^\n]");
	getc(INFILE->f);
      }
      if (!BUFEOF(INFILE->f)) {
	if (P_peek(INFILE->f) == '>') {   /*read past '>' */
	  CH = getc(INFILE->f);
	  if (CH == '\n')
	    CH = ' ';
	}
	NUMSEQ++;
	READWORD(INFILE, &NAMES[NUMSEQ]);   /* read name */
	if (!BUFEOF(INFILE->f)) {
	  fscanf(INFILE->f, "%*[^\n]");
	  getc(INFILE->f);
	}
      }
    }  /* FORMAT in [IG,PEARSON] */
    else
      NUMSEQ++;

    /* Read sequence */
    ENDOFSEQUENCE = false;
    while (!ENDOFSEQUENCE) {
      while (!P_eoln(INFILE->f)) {
	CH = getc(INFILE->f);
	if (CH == '\n')
	  CH = ' ';
	/* Truncate alignment if > MAXLEN */
	if (SEQLEN[NUMSEQ] < MAXLEN) {
	  SEQLEN[NUMSEQ]++;
	  if (CH == '-' || CH == ' ')
	    CH = GAPS;
	  ALIGNMENT[NUMSEQ][SEQLEN[NUMSEQ] - 1] = TOLOWER(CH);
	}
      }
      if (!BUFEOF(INFILE->f)) {
	fscanf(INFILE->f, "%*[^\n]");
	getc(INFILE->f);
      }
      if (BUFEOF(INFILE->f) || FORMAT == RALIGN)
	ENDOFSEQUENCE = true;
      else if (P_peek(INFILE->f) == '>' || P_peek(INFILE->f) == ';')
	ENDOFSEQUENCE = true;
      if (ENDOFSEQUENCE) {
	if (SEQLEN[NUMSEQ] > LONGEST)
	  LONGEST = SEQLEN[NUMSEQ];
      }
    }  /* not ENDOFSEQUENCE */
  }  /* not eof(INFILE) */
  SEQLEN[0] = LONGEST;
}  /* READALIGNMENT */


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


/* Return the majority character position J, or if there is no majority,
   "n" for DNA and "x" for proteins */
Local Char CONSENSUS(J)
long J;
{
  long I, K;
  long SETSIZE = 1, KNOWNCOUNT = 0;
  long MAXNUMBER;
  Char LETTER[MAXSEQ];
  long NUMBER[MAXSEQ];
  long FORLIM;

  /* Count the number of each different type of character at position
     J. SETSIZE=number of different characters examined so far. */
  LETTER[0] = ALIGNMENT[1][J-1];
  NUMBER[0] = 1;
  FORLIM = NUMSEQ;
  for (I = 2; I <= FORLIM; I++) {
    K = 1;
    while (ALIGNMENT[I][J-1] != LETTER[K-1] && K <= SETSIZE)
      K++;
    if (K > SETSIZE) {
      LETTER[K-1] = ALIGNMENT[I][J-1];
      NUMBER[K-1] = 1;
      SETSIZE++;
    } else
      NUMBER[K-1]++;
  }  /* for I */

  /* Set CONSENSUS equal to the majority character, if there is one */
  if (LETTER[0] == UNKNOWN)
    I = 2;
  else
    I = 1;
  MAXNUMBER = I;
  for (K = I; K <= SETSIZE; K++) {
    if (LETTER[K-1] != UNKNOWN) {
      if (NUMBER[K-1] > NUMBER[MAXNUMBER-1])
	MAXNUMBER = K;
      KNOWNCOUNT += NUMBER[K-1];
    }
  }
  if (NUMBER[MAXNUMBER-1] > KNOWNCOUNT / 2)
    return (LETTER[MAXNUMBER-1]);
  else
    return UNKNOWN;
}  /* CONSENSUS */


/* END MODULE TOUPPER         VERSION= 'SUNMODS     Version  8/ 9/94'; */

/***************************************************************/
/* Modify the alignment according to the command line options. */
/***************************************************************/
Static Void MODIFY()
{
  long I, J;
  Char CONCHAR;
  boolean AGREE;
  long FORLIM, FORLIM1;

  FORLIM = NUMSEQ;
  /* The consensus will be the length of the longest sequence. Any
     sequence shorter than the consensus will be padded with gaps. */
  for (I = 1; I <= FORLIM; I++) {
    FORLIM1 = SEQLEN[0];
    for (J = SEQLEN[I]; J < FORLIM1; J++)
      ALIGNMENT[I][J] = GAPS;
  }

  FORLIM = SEQLEN[0];
  /* For each position, calculate the consensus */
  for (J = 0; J < FORLIM; J++) {
    CONCHAR = CONSENSUS(J + 1);
    ALIGNMENT[0][J] = CONCHAR;
    AGREE = true;
    FORLIM1 = NUMSEQ;
    for (I = 1; I <= FORLIM1; I++) {
      if (ALIGNMENT[I][J] == CONCHAR) {
	if (PERIODS && CONCHAR != UNKNOWN)
	  ALIGNMENT[I][J] = '.';
      } else {
	AGREE = false;
	if (ALIGNMENT[I][J] == '-' || ALIGNMENT[I][J] == ' ')
	  ALIGNMENT[I][J] = GAPS;
      }
    }  /* for I */
    if (AGREE && CAPS)
      ALIGNMENT[0][J] = TOUPPER(ALIGNMENT[0][J]);
  }  /* for J */
}  /* MODIFY */


/***************************************************************/
/* Write the alignment to OUTFILE.                             */
/***************************************************************/
Static Void WRITEALIGNMENT(OUTFILE)
_TEXT *OUTFILE;
{
  long NUCSPRINTED[MAXSEQ + 1];
  boolean DONE = false;   /* true when all of each sequence has been printed*/
  long I, J, THISLINE;   /* index of last nucleotide printed on a line */
  long NUMBER, FORLIM;

  FORLIM = NUMSEQ;
  for (I = 0; I <= FORLIM; I++)
    NUCSPRINTED[I] = 0;
  NUMBER = STARTNUM - 1;   /* number used for printing */
  while (!DONE) {
    DONE = true;

    /* Write a line of numbers */
    fprintf(OUTFILE->f, "%10c", ' ');
    THISLINE = NUCSPRINTED[0] + LINELENGTH;
    if (THISLINE > SEQLEN[0])
      THISLINE = SEQLEN[0];
    FORLIM = (THISLINE - NUCSPRINTED[0]) / 10;
    for (J = 1; J <= FORLIM; J++) {
      /* The coordinate 0 is not skipped. */
      if (NUMBER >= -10 && NUMBER < 0)
	NUMBER += 11;
      else
	NUMBER += 10;
      fprintf(OUTFILE->f, "%10ld", NUMBER);
    }
    putc('\n', OUTFILE->f);

    FORLIM = NUMSEQ;
    /* Write LINELENGTH nucleotides for each position */
    for (I = 0; I <= FORLIM; I++) {
      WRITEWORD(OUTFILE, NAMES[I], 10L);
      THISLINE = NUCSPRINTED[I] + LINELENGTH;
      if (THISLINE > SEQLEN[I])
	THISLINE = SEQLEN[I];
      for (J = NUCSPRINTED[I]; J < THISLINE; J++)
	putc(ALIGNMENT[I][J], OUTFILE->f);
      NUCSPRINTED[I] = THISLINE;
      putc('\n', OUTFILE->f);
      if (NUCSPRINTED[I] < SEQLEN[I])
	DONE = false;
    }  /* for I */
    putc('\n', OUTFILE->f);
  }  /* not DONE */
}  /* WRITEALIGNMENT */


/* ----------------------------------------------------------  */
/* ----------------- MAIN  PROCEDURE  -----------------------  */
main(argc, argv)
int argc;
Char *argv[];
{
  _TEXT TEMP, TEMP1;

  PASCAL_MAIN(argc, argv);
  SEQFILE.f = NULL;
  *SEQFILE.name = '\0';
  /* Read options from command line */
  GAPS = ' ';
  CAPS = false;
  PERIODS = false;
  UNKNOWN = 'x';
  FORMAT = RALIGN;
  LINELENGTH = 70;
  STARTNUM = 1;
  ARGNUM = STARTARGNUM;
  FILEARGUMENT = false;
  TOTALARGS = P_argc - 1;

  while (ARGNUM <= TOTALARGS && !FILEARGUMENT) {
    P_sun_argv(ARGUMENT, 132, (int)ARGNUM);
    if (ARGUMENT[0] != '-') {
      FILEARGUMENT = true;
      break;
    }  /* ARGUMENT[1]='-' */
    if (ARGUMENT[1] == 's' || ARGUMENT[1] == 'l' || ARGUMENT[1] == 'm' ||
	ARGUMENT[1] == 'p' || ARGUMENT[1] == 'f' || ARGUMENT[1] == 'n' ||
	ARGUMENT[1] == 'c' || ARGUMENT[1] == 'g') {
      switch (ARGUMENT[1]) {

      case 'g':
	GAPS = '-';
	break;

      case 'p':
	PERIODS = true;
	break;

      case 'c':
	CAPS = true;
	break;

      case 'n':
	UNKNOWN = 'n';
	break;

      case 'f':
	if (ARGUMENT[2] == 'p' || ARGUMENT[2] == 'i' || ARGUMENT[2] == 'r') {
	  switch (ARGUMENT[2]) {

	  case 'r':
	    FORMAT = RALIGN;
	    break;

	  case 'i':
	    FORMAT = IG;
	    break;

	  case 'p':
	    FORMAT = PEARSON;
	    break;
	  }
	}
	break;

      case 'm':
	FORMAT = PEARSON;
	break;

      case 's':
	STARTNUM = 0;
	I = 3;
	if (ARGUMENT[I-1] == '-') {
	  SIGN = -1;
	  I++;
	} else
	  SIGN = 1;
	while (isdigit(ARGUMENT[I-1])) {
	  STARTNUM = STARTNUM * 10 + ARGUMENT[I-1] - '0';
	  I++;
	}
	STARTNUM *= SIGN;
	break;
	/* s */

      case 'l':
	LINELENGTH = 0;
	I = 3;
	while (isdigit(ARGUMENT[I-1])) {
	  LINELENGTH = LINELENGTH * 10 + ARGUMENT[I-1] - '0';
	  I++;
	}
	if (LINELENGTH < 1)
	  LINELENGTH = 70;
	break;
	/* l */
      }
    }
    ARGNUM++;
  }

  /* Read sequence names from SEQFILE. (RALIGN output only) */
  if (FORMAT == RALIGN) {
    FILEARGS(&SEQFILE, 'I', &ARGNUM);
    I = 0;
    NAMES[0].LEN = 0;   /* consensus name is blank */
    while (!BUFEOF(SEQFILE.f)) {
      /* skip comment lines */
      while (P_peek(SEQFILE.f) == ';') {
	fscanf(SEQFILE.f, "%*[^\n]");
	getc(SEQFILE.f);
      }
      if (BUFEOF(SEQFILE.f))
	break;
      /* read the name */
      I++;
      READWORD(&SEQFILE, &NAMES[I]);
      if (!BUFEOF(SEQFILE.f)) {
	fscanf(SEQFILE.f, "%*[^\n]");
	getc(SEQFILE.f);
      }
      /* skip the sequence */
      do {
	READLINE(&SEQFILE, &DUMMY);
      } while (!((DUMMY.LEN == 0) | BUFEOF(SEQFILE.f)));
	  /* if not eof(SEQFILE) */
    }
  }  /* RALIGN */

  TEMP.f = stdin;
  *TEMP.name = '\0';
  TEMP1.f = stdout;
  *TEMP1.name = '\0';
  /* Read in the aligned sequences produced by ralign.*/
  READALIGNMENT(&TEMP, &TEMP1);

  /* Modify the output according to command line options.*/
  MODIFY();

  TEMP.f = stdout;
  *TEMP.name = '\0';
  /* Write the re-formatted alignment with numbering.*/
  WRITEALIGNMENT(&TEMP);
  if (SEQFILE.f != NULL)
    fclose(SEQFILE.f);
  exit(EXIT_SUCCESS);
}  /* REFORM  */




/* End. */
