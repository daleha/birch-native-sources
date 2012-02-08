/* Output from p2c 1.21alpha-07.Dec.93, the Pascal-to-C translator */
/* From input file "shuffle.p" */


/* ********************************************************  */
/*                                                           */
/*   SHUFFLE   VERSION 10/23/94, Standard Pascal             */
/*             Brian Fristensky                              */
/*             Dept. of Plant Science                        */
/*             University of Manitoba                        */
/*             Winnipeg, MB  R3T 2N2 CANADA                  */
/*                                                           */
/* SYNOPSIS                                                  */
/*    shuffle -s<n> [-w<n> -o<n>]                            */
/*                                                           */
/* DESCRIPTION                                               */
/*    Shuffles nucleici acid or protein sequences in a       */
/*    Pearson (.wrp) format file. see Lipman et al           */
/*                                                           */
/*    -s<n>  n is a random integer between 0 and 32767.      */
/*           This number must be provided for each run.      */
/*                                                           */
/*    -w<n>  n is an integer, indicating the width of the    */
/*           window for random localization. It should       */
/*           never exceed the length of the shortest input   */
/*           sequence.                                       */
/*                                                           */
/*    -o<n>  n is an integer, indicating the number of nuc-  */
/*           leotides overlap between adjacent windows. It   */
/*           should never exceed the window size.  o def-    */
/*           aults to 0 of not specified.                    */
/*                                                           */
/*    If w and o are specified, overlapping windows of w     */
/*    nucleotides are shuffled, thus preserving the local    */
/*    characteristic base composition. Windows overlap by    */
/*    o nucleotides.                                         */
/*                                                           */
/*    If w and o are not specified, each sequence is shuf-   */
/*    fled globally, thus preserving the overall base        */
/*    composition, but not the local variations in comp.     */
/*                                                           */
/* Copyright (c) 1988, 1990  by Brian Fristensky.            */
/* !!! in comment indicates feature which may need change.   */
/*  *******************************************************  */

/*!!!*/

#include <p2c.h>


/*!!! Some Pascals require file parameters in heading */

#define MAXSEQ          750000L

#define MAXWORD         25
#define MAXLINE         120
#define STARTARGNUM     1
/* BEGIN MODULE REALLIMITS */
/*!!!*/

#define MAXREAL         1.7e38
/* END MODULE REALLIMITS         VERSION= 'SUNMODS     Version  8/ 9/94'; */

#define VERSION         "SHUFFLE       Version 10/23/94"


typedef enum {
  T, U, C, A, G, N, R, Y, M, W, S, K, D, H, V, B, L, Z, F, P, E, I, Q,
  ASTERISK, X, GAP
} SEQCHAR;
typedef SEQCHAR SEQUENCE[MAXSEQ];

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


Static SEQUENCE SEQ;
Static long SEQLEN;
Static Char SCHAR[26];
/* global processing parameters */
Static long ARGNUM, POS;
Static boolean UNREADOPTIONS;
Static CHARARRAY ARGUMENT;
Static long RSEED, OVERLAP, WINDOW, J;


/* BEGIN MODULE NUMBER */
/* Extract an integer from a CHARARRAY (see TYPE.LINE), starting
   at the current position. */
Static long NUMBER(TARGET, POS)
Char *TARGET;
long *POS;
{
  long N_ = 0, ORDZERO = '0';

  /* evaluate characteristic */
  while (isdigit(TARGET[*POS - 1])) {
    N_ = N_ * 10 + TARGET[*POS - 1] - ORDZERO;
    (*POS)++;
  }
  return N_;
}  /* NUMBER */


/* END MODULE NUMBER         VERSION= 'SUNMODS     Version  8/ 9/94'; */

/* BEGIN MODULE READLINE */
/* Read a line from a file, omitting trailing blanks */
Static Void READLINE(F_, L_)
_TEXT *F_;
LINE *L_;
{
  long LASTNONBLANK = 0;
  long I_;
  Char CH;

  /* with L*/
  L_->LEN = 0;
  while (!P_eoln(F_->f)) {
    CH = getc(F_->f);
    if (CH == '\n')
      CH = ' ';
    if (L_->LEN < MAXLINE) {
      L_->LEN++;
      L_->STR[L_->LEN - 1] = CH;
      if (CH != ' ')
	LASTNONBLANK = L_->LEN;
    }
  }
  if (!BUFEOF(F_->f)) {
    fscanf(F_->f, "%*[^\n]");
    getc(F_->f);
  }
  L_->LEN = LASTNONBLANK;
  for (I_ = L_->LEN; I_ < MAXLINE; I_++)
    L_->STR[I_] = ' ';
}  /* READLINE */


/* END MODULE READLINE         VERSION= 'SUNMODS     Version  8/ 9/94'; */

/* BEGIN MODULE WRITELINE */
/*  Write a line to a file using L char, left-justified.  */
Static Void WRITELINE(F_, W_, L_)
_TEXT *F_;
LINE W_;
long L_;
{
  long I_;

  for (I_ = 1; I_ <= L_; I_++) {
    if (I_ <= W_.LEN)
      putc(W_.STR[I_-1], F_->f);
    else
      putc(' ', F_->f);
  }
}  /* WRITELINE */


/* END MODULE WRITELINE         VERSION= 'SUNMODS     Version  8/ 9/94'; */

/* BEGIN MODULE NUCPROT */
/*****************************************************************/
/*  Convert a character to the appropriate sequence symbol.      */
/*****************************************************************/
Static SEQCHAR NUCPROT(CH)
Char CH;
{
  SEQCHAR Result;

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

  case 'U':
  case 'u':
    Result = U;
    break;

  case 'L':
  case 'l':
    Result = L;
    break;

  case 'T':
  case 't':
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

  case 'Z':
  case 'z':
    Result = Z;
    break;

  case 'F':
  case 'f':
    Result = F;
    break;

  case 'P':
  case 'p':
    Result = P;
    break;

  case 'E':
  case 'e':
    Result = E;
    break;

  case 'I':
  case 'i':
    Result = I;
    break;

  case 'Q':
  case 'q':
    Result = Q;
    break;

  case '*':
    Result = ASTERISK;
    break;

  case 'X':
  case 'x':
    Result = X;
    break;

  case '-':
    Result = GAP;
    break;
  }
  return Result;
}


/* END MODULE NUCPROT */

/*****************************************************/
/* Initialization procedures.                        */
/*****************************************************/
Static Void INITIALIZE()
{
  /*Set default values for parameters */
  WINDOW = 0;
  OVERLAP = 0;
  RSEED = 0;

  /*   Initialize SCHAR array which holds the          */
  /*   character values of the SEQCHARS.               */
  SCHAR[(long)A] = 'A';
  SCHAR[(long)C] = 'C';
  SCHAR[(long)G] = 'G';
  SCHAR[(long)T] = 'T';
  SCHAR[(long)U] = 'U';
  SCHAR[(long)L] = 'L';
  SCHAR[(long)R] = 'R';
  SCHAR[(long)Y] = 'Y';
  SCHAR[(long)M] = 'M';
  SCHAR[(long)W] = 'W';
  SCHAR[(long)S] = 'S';
  SCHAR[(long)K] = 'K';
  SCHAR[(long)D] = 'D';
  SCHAR[(long)H] = 'H';
  SCHAR[(long)V] = 'V';
  SCHAR[(long)B] = 'B';
  SCHAR[(long)N] = 'N';
  SCHAR[(long)Z] = 'Z';
  SCHAR[(long)F] = 'F';
  SCHAR[(long)P] = 'P';
  SCHAR[(long)E] = 'E';
  SCHAR[(long)I] = 'I';
  SCHAR[(long)Q] = 'Q';
  SCHAR[(long)X] = 'X';
  SCHAR[(long)ASTERISK] = '*';
  SCHAR[(long)GAP] = '-';

}  /* INITIALIZE */


/* BEGIN MODULE RANDSEED */
Static Void RANDSEED(X_)
long X_;
{
  srand(X_);
/* p2c: shuffle.p, line 190: Warning: Symbol 'srand' is not defined [221] */
}  /* RANDINT */


/* END MODULE RANDSEED */

/* BEGIN MODULE RANDREAL */

/* It's real hard to make random number generation portable, because the
   C rand() function generates random integers in an architecture dependent
   range, usually in a range between 0 and (2**16)-1 or (2**32)-1.*/
Static double RANDREAL()
{
  long X_;

  X_=rand();
  /* if  random integer is a  31bit integer, perform an arithmetic shift
     right 16 bits to make it a 16 bit integer. */
  if (X_ > 65536L)
    X_ /= 65536L;

  /* Convert random integer into random real between 0 and 1.*/
  return (X_ / 32768.0);
}  /* RANDREAL */


/* Local variables for SCRAMBLE: */
struct LOC_SCRAMBLE {
  double RNUM;
} ;

/* - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/* Read in the next sequence                            */
Local Void READSEQUENCE(SEQFILE, SEQ, LINK)
_TEXT *SEQFILE;
SEQCHAR *SEQ;
struct LOC_SCRAMBLE *LINK;
{
  Char CH;
  boolean TOOLONG = false;

  SEQLEN = 0;
  while ((!BUFEOF(SEQFILE->f)) & (P_peek(SEQFILE->f) != ';' &&
				  P_peek(SEQFILE->f) != '>')) {
    while (!P_eoln(SEQFILE->f)) {
      CH = getc(SEQFILE->f);
      if (CH == '\n')
	CH = ' ';
      if (CH == '-' || CH == 'x' || CH == 'X' || CH == '*' || CH == 'q' ||
	  CH == 'Q' || CH == 'i' || CH == 'I' || CH == 'e' || CH == 'E' ||
	  CH == 'p' || CH == 'P' || CH == 'f' || CH == 'F' || CH == 'z' ||
	  CH == 'Z' || CH == 'b' || CH == 'B' || CH == 'v' || CH == 'V' ||
	  CH == 'h' || CH == 'H' || CH == 'd' || CH == 'D' || CH == 'l' ||
	  CH == 'L' || CH == 'k' || CH == 'K' || CH == 's' || CH == 'S' ||
	  CH == 'w' || CH == 'W' || CH == 'm' || CH == 'M' || CH == 'y' ||
	  CH == 'Y' || CH == 'r' || CH == 'R' || CH == 'n' || CH == 'N' ||
	  CH == 'u' || CH == 'U' || CH == 't' || CH == 'T' || CH == 'g' ||
	  CH == 'G' || CH == 'c' || CH == 'C' || CH == 'a' || CH == 'A') {
/* p2c: shuffle.p, line 350: Note:
 * Line breaker spent 1.8+0.64 seconds, 5000 tries on line 428 [251] */
	if (SEQLEN < MAXSEQ - 2) {
	  SEQLEN++;
	  SEQ[SEQLEN-1] = NUCPROT(CH);
	} else
	  TOOLONG = true;
      }
    }  /* eoln */
    fscanf(SEQFILE->f, "%*[^\n]");
    getc(SEQFILE->f);
  }  /* eof */
  if (TOOLONG)
    printf(">>> WARNING! Sequence length exceeds MAXSEQ-2. Seq. truncated.\n");
}  /* READSEQUENCE */

/* Local variables for RANDOMIZE: */
struct LOC_RANDOMIZE {
  struct LOC_SCRAMBLE *LINK;
  SEQCHAR *SEQ;
  long WINDOW;
} ;

/* Swap nucleotides or amino acids between two positions. */
Local Void SWAP(N1, N2)
SEQCHAR *N1, *N2;
{
  SEQCHAR TEMP;

  TEMP = *N1;
  *N1 = *N2;
  *N2 = TEMP;
}  /* SWAP */

/* Randomize nucleotides/aa's within a window by swapping. */
Local Void SCRAMBLE_WINDOW(LEFT, RIGHT, LINK)
long LEFT, RIGHT;
struct LOC_RANDOMIZE *LINK;
{
  long I_, J;

  for (I_ = LEFT - 1; I_ < RIGHT; I_++) {
    LINK->LINK->RNUM = RANDREAL();
    J = LEFT + (long)floor(LINK->LINK->RNUM * LINK->WINDOW + 0.5);
    if (J > RIGHT)   /* Make sure to stay within bounds */
      J = RIGHT;
    SWAP(&LINK->SEQ[I_], &LINK->SEQ[J-1]);
  }
}  /* SCRAMBLE_WINDOW */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
Local Void RANDOMIZE(SEQ_, WINDOW_, LINK)
SEQCHAR *SEQ_;
long WINDOW_;
struct LOC_SCRAMBLE *LINK;
{
  struct LOC_RANDOMIZE V;
  long LEFT = 1;
  long RIGHT, OFFSET;

  V.LINK = LINK;
  V.SEQ = SEQ_;
  V.WINDOW = WINDOW_;
  if (V.WINDOW < 1)
    V.WINDOW = SEQLEN;
  else if (V.WINDOW > SEQLEN)
    V.WINDOW = SEQLEN;
  RIGHT = V.WINDOW;
  OFFSET = V.WINDOW - OVERLAP;

  while (RIGHT <= SEQLEN) {
    SCRAMBLE_WINDOW(LEFT, RIGHT, &V);
    LEFT += OFFSET;
    RIGHT += OFFSET;
    if (RIGHT <= SEQLEN) {  /* end of sequence, a special case */
      continue;
    }  /* RIGHT > SEQLEN */
    RIGHT = SEQLEN;
    V.WINDOW = RIGHT - LEFT + 1;
    if (V.WINDOW > 1)
      SCRAMBLE_WINDOW(LEFT, RIGHT, &V);
    RIGHT += OFFSET;
  }  /* RIGHT <= SEQLEN */
}  /* RANDOMIZE */

/* - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
Local Void WRITESEQ(F_, SEQ, LINK)
_TEXT *F_;
SEQCHAR *SEQ;
struct LOC_SCRAMBLE *LINK;
{
  long SEQPRINTED = 0;
  long THISLINE;

  while (SEQPRINTED < SEQLEN) {
    THISLINE = SEQPRINTED + 50;
    if (THISLINE > SEQLEN)
      THISLINE = SEQLEN;
    while (SEQPRINTED < THISLINE) {
      SEQPRINTED++;
      putc(SCHAR[(long)SEQ[SEQPRINTED-1]], F_->f);
    }
    putc('\n', F_->f);
  }  /* SEQPRINTED < SEQLEN */
}  /* WRITESEQ */


/* END MODULE RANDREAL */

/* ******************************************************************* */
/*  Randomize each sequence by shuffling groups of WINDOW bases .      */
/*  Each adjecent window overlaps by OVERLAP positions.                */
/* ******************************************************************* */
Static Void SCRAMBLE(SEQFILE, OUTFILE)
_TEXT *SEQFILE, *OUTFILE;
{  /* ----------------- SCRAMBLE ----------------------- */
  struct LOC_SCRAMBLE V;
  LINE COMMENTLINE;

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
  /* The next two output lines appear as comments in the first two
     lines of the output.*/
  fprintf(OUTFILE->f, ">%s\n", VERSION);
  fprintf(OUTFILE->f,
	  "> RANDOM SEED: %10ld%10cWINDOW: %5ld%10cOVERLAP: %5ld\n",
	  RSEED, ' ', WINDOW, ' ', OVERLAP);

  RANDSEED(RSEED);
  while (!BUFEOF(SEQFILE->f)) {
    /* Read in comment lines */
    while (P_peek(SEQFILE->f) == ';' || P_peek(SEQFILE->f) == '>') {
      READLINE(SEQFILE, &COMMENTLINE);
      WRITELINE(OUTFILE, COMMENTLINE, COMMENTLINE.LEN);
      putc('\n', OUTFILE->f);
    }  /* SEQFILE^ in ['>',';'] */

    /* Process a sequence */
    if (BUFEOF(SEQFILE->f))
      break;
    READSEQUENCE(SEQFILE, SEQ, &V);
    RANDOMIZE(SEQ, WINDOW, &V);
    WRITESEQ(OUTFILE, SEQ, &V);
  }  /* not eof(SEQFILE) */

}  /* SCRAMBLE */


/* ----------------------------------------------------------  */
/* ----------------- MAIN  PROCEDURE  -----------------------  */
main(argc, argv)
int argc;
Char *argv[];
{
  _TEXT TEMP, TEMP1;

  PASCAL_MAIN(argc, argv);
  INITIALIZE();

  /* Read options from the command line. Note that -s is required */
  ARGNUM = STARTARGNUM;
  UNREADOPTIONS = true;
  while (UNREADOPTIONS) {
    if (ARGNUM >= P_argc) {
      UNREADOPTIONS = false;
      break;
    }  /* if ARGNUM <= argc */
    P_sun_argv(ARGUMENT, sizeof(CHARARRAY), (int)ARGNUM);
    if (ARGUMENT[0] == '-') {
      if (ARGUMENT[1] == 'o' || ARGUMENT[1] == 'w' || ARGUMENT[1] == 's') {
	POS = 3;
	switch (ARGUMENT[1]) {

	case 's':
	  RSEED = NUMBER(ARGUMENT, &POS);
	  break;

	case 'w':
	  WINDOW = NUMBER(ARGUMENT, &POS);
	  break;

	case 'o':
	  OVERLAP = NUMBER(ARGUMENT, &POS);
	  break;
	}
      }  /* s,w,o */
    }  /* '-' */
    else
      UNREADOPTIONS = false;
    ARGNUM++;
  }

  /* Scramble the sequence(s) as specified. */
  if (RSEED > 0) {
    TEMP.f = stdin;
    *TEMP.name = '\0';
    TEMP1.f = stdout;
    *TEMP1.name = '\0';
    SCRAMBLE(&TEMP, &TEMP1);
  } else
    printf("ERROR >>> Random seed (-s) must be set in command line.\n");
  exit(EXIT_SUCCESS);
}  /* SHUFFLE */



/* End. */
