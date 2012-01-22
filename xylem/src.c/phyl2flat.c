/* Output from p2c 1.21alpha-07.Dec.93, the Pascal-to-C translator */
/* From input file "phyl2flat.p" */


/* ********************************************************  */
/*                                                           */
/*  PHYL2FLAT  VERSION  3/12/97, Standard Pascal             */
/*             Brian Fristensky                              */
/*             Dept. of Plant Science                        */
/*             University of Manitoba                        */
/*             Winnipeg, MB  R3T 2N2 CANADA                  */
/*                                                           */
/* SYNOPSIS                                                  */
/*    phyl2flat [-i] < Phylip_file > GDE_flatfile            */
/*                                                           */
/* DESCRIPTION                                               */
/*    Convert a phylip discrete character file into a        */
/*    GDE flatfile.                                          */
/*                                                           */
/*             -i  invert characters, so that 0 -> 1 and     */
/*                 1 -> 0.                                   */
/*             -n  read non-interleaved format               */
/*                                                           */
/* Copyright (c) 1996        by Brian Fristensky.            */
/* !!! in comment indicates feature which may need change.   */
/*  *******************************************************  */

/*!!!*/

#include <p2c.h>


/*!!! Some Pascals require file parameters in heading */

#define MAXCHAR         10000
#define MAXSEQ          500
#define MAXWORD         10
#define STARTARGNUM     1
/* BEGIN MODULE REALLIMITS */
/*!!!*/

#define MAXREAL         1.7e38
/* END MODULE REALLIMITS         VERSION= 'SUNMODS     Version  8/ 9/94'; */

#define VERSION         "PHYL2FLAT     Version 3/12/97"


/* BEGIN MODULE TYPE.WORD */
/*   <word>::= <non-blank char>[<non-blank char>] */

typedef struct WORD {
  long LEN;
  Char STR[MAXWORD];
} WORD;

/* END MODULE TYPE.WORD         VERSION= 'SUNMODS     Version  8/ 9/94'; */

typedef WORD NAMEARRAY[MAXSEQ];
typedef Char SEQARRAY[MAXSEQ][MAXCHAR];
typedef long INTARRAY[MAXSEQ];


Static SEQARRAY SEQ;
Static NAMEARRAY NAME;
Static INTARRAY LENGTH;
Static long NUMSEQ, SEQLEN;
/* global processing parameters */
Static long ARGNUM, POS;
Static boolean OKAY, INTERLEAVED, INVERT, UNREADOPTIONS;
Static Char ARGUMENT[132];   /* command line argument */


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


/* --------------------------------------------*/
Local Void READNAME(F, NAME, I, OKAY)
_TEXT *F;
WORD *NAME;
long I;
boolean *OKAY;
{
  long J = 1;
  Char CH;

  while (J <= MAXWORD && *OKAY) {
    if (P_eoln(F->f)) {
      *OKAY = false;
      printf("PHYL2FLAT: >>>Truncated sequence line.\n");
      continue;
    }
    CH = getc(F->f);
    if (CH == '\n')
      CH = ' ';
    if (CH != ' ') {
      NAME[I-1].STR[J-1] = CH;
      NAME[I-1].LEN++;
    }  /* CH <> ' ' */
    J++;
  }
}  /* READNAME */

/* --------------------------------------------*/
Local Void READSEQLINE(F, SEQ, LENGTH, I)
_TEXT *F;
Char (*SEQ)[MAXCHAR];
long *LENGTH;
long I;
{
  long K;
  Char CH;

  /* Read in the sequence */
  K = LENGTH[I-1];
  while (!P_eoln(F->f) && K < MAXCHAR) {
    if (P_peek(F->f) == 'B' || P_peek(F->f) == 'P' || P_peek(F->f) == '?' ||
	P_peek(F->f) == '-' || P_peek(F->f) == '+' || P_peek(F->f) == '1' ||
	P_peek(F->f) == '0') {
      K++;
      SEQ[I-1][K-1] = getc(F->f);
      if (SEQ[I-1][K-1] == '\n')
	SEQ[I-1][K-1] = ' ';
    } else {
      CH = getc(F->f);
      if (CH == '\n')
	CH = ' ';
    }
  }
  if (!BUFEOF(F->f)) {
    fscanf(F->f, "%*[^\n]");
    getc(F->f);
  }
  LENGTH[I-1] = K;
}  /* READSEQLINE */

/* --------------------------------------------*/
/* Calculate the minimum and maximum sequence lengths read so far. */
Local Void MINMAXREAD(LENGTH, MINREAD, MAXREAD)
long *LENGTH;
long *MINREAD, *MAXREAD;
{
  long I, FORLIM;

  *MINREAD = LENGTH[0];
  *MAXREAD = *MINREAD;
  FORLIM = NUMSEQ;
  for (I = 1; I < FORLIM; I++) {
    if (LENGTH[I] < *MINREAD)
      *MINREAD = LENGTH[I];
    if (LENGTH[I] > *MAXREAD)
      *MAXREAD = LENGTH[I];
  }
}  /* MINMAXREAD */


/* END MODULE WRITEWORD         VERSION= 'SUNMODS     Version  8/ 9/94'; */

/* ******************************************************************* */
/* Read in the  names and sequences.                                   */
/* ******************************************************************* */
Static Void NAMESANDSEQ(F, NAME, SEQ, SEQLEN, INTERLEAVED, OKAY)
_TEXT *F;
WORD *NAME;
Char (*SEQ)[MAXCHAR];
long *SEQLEN;
boolean INTERLEAVED, *OKAY;
{  /* NAMESANDSEQ -------------------------------------- */
  long I;
  long MINREAD = 0, MAXREAD = 0;
  long FORLIM;

  if (*F->name != '\0') {
    if (F->f != NULL)
      F->f = freopen(F->name, "r", F->f);
    else
      F->f = fopen(F->name, "r");
  } else
    rewind(F->f);
  if (F->f == NULL)
    _EscIO2(FileNotFound, F->name);
  RESETBUF(F->f, Char);

  /* Read line 1, telling number of 'sequences' and number of
     characters per sequence */
  fscanf(F->f, "%ld%ld%*[^\n]", &NUMSEQ, SEQLEN);
  getc(F->f);
  FORLIM = NUMSEQ;
  for (I = 0; I < FORLIM; I++) {
    NAME[I].LEN = 0;
    LENGTH[I] = 0;
  }

  /* Read in the names and sequences */
  *OKAY = true;

  if (INTERLEAVED) {
    while (MAXREAD < *SEQLEN && *OKAY) {
      I = 1;
      while (I <= NUMSEQ && *OKAY) {
	/* First set of seq. lines, read in the name.
	   In interleaved format, as many additional sets of
	   sequence lines, in the same order as the first,
	   will follow, without names. */
	if (NAME[I-1].LEN == 0)
	  READNAME(F, NAME, I, OKAY);
	READSEQLINE(F, SEQ, LENGTH, I);
	I++;
      }  /* while I <= NUMSEQ */
      MINMAXREAD(LENGTH, &MINREAD, &MAXREAD);
    }  /* while MAXREAD < SEQLEN */

    if (MINREAD >= MAXREAD)
      return;
    printf(">>> PHYL2FLAT: sequences not all same length\n");
    printf(">>> truncating to %12ld\n", MINREAD);
    *SEQLEN = MINREAD;
    return;
  }  /* if INTERLEAVED */
  I = 1;
  while (I <= NUMSEQ && *OKAY) {
    /* First set of seq. lines, read in the name.
       In interleaved format, as many additional sets of
       sequence lines, in the same order as the first,
       will follow, without names. */
    if (NAME[I-1].LEN == 0)
      READNAME(F, NAME, I, OKAY);
    while ((LENGTH[I-1] < *SEQLEN) & (!BUFEOF(F->f)))
      READSEQLINE(F, SEQ, LENGTH, I);
    I++;
    /* non-interleaved format */
    /* non-interleaved format */

  }  /* while (I <= NUMSEQ) and OKAY */
}  /* NAMESANDSEQ */


/* ******************************************************************* */
/* Invert the sequence ie. 1->0 and 0->1                               */
/* ******************************************************************* */
Static Void INVERTSEQ(SEQ)
Char (*SEQ)[MAXCHAR];
{
  long I, K, FORLIM, FORLIM1;

  FORLIM = NUMSEQ;
  for (I = 0; I < FORLIM; I++) {
    FORLIM1 = SEQLEN;
    for (K = 0; K < FORLIM1; K++) {
      if (SEQ[I][K] == '-' || SEQ[I][K] == '+' || SEQ[I][K] == '1' ||
	  SEQ[I][K] == '0') {
	switch (SEQ[I][K]) {

	case '0':
	  SEQ[I][K] = '1';
	  break;

	case '1':
	  SEQ[I][K] = '0';
	  break;

	case '-':
	  SEQ[I][K] = '+';
	  break;

	case '+':
	  SEQ[I][K] = '-';
	  break;
	}
      }
    }
  }
}  /* INVERTSEQ */


#define LINELEN         50


/* ******************************************************************* */
/* Write names and sequences in GDE flatfile format.                   */
/* ******************************************************************* */
Static Void WRITEFLAT(F, NUMSEQ, SEQLEN, NAME, SEQ)
_TEXT *F;
long NUMSEQ, SEQLEN;
WORD *NAME;
Char (*SEQ)[MAXCHAR];
{
  long I, K, THISLINE;

  for (I = 0; I < NUMSEQ; I++) {
    putc('"', F->f);   /* indicates text in GDE flatfile */
    WRITEWORD(F, NAME[I], NAME[I].LEN);
    putc('\n', F->f);
    K = 1;
    THISLINE = 0;
    while (K <= SEQLEN) {
      THISLINE += LINELEN;
      if (THISLINE > SEQLEN)
	THISLINE = SEQLEN;
      while (K <= THISLINE) {
	putc(SEQ[I][K-1], F->f);
	K++;
      }  /* K <= THISLINE */
      putc('\n', F->f);
    }  /* K <= SEQLEN*/
  }  /* for I */
}  /* WRITESEQ */

#undef LINELEN


/* ----------------------------------------------------------  */
/* ----------------- MAIN  PROCEDURE  -----------------------  */
main(argc, argv)
int argc;
Char *argv[];
{
  _TEXT TEMP;

  PASCAL_MAIN(argc, argv);
  /* Read options from the command line. Note that -s is required */
  INVERT = false;
  INTERLEAVED = true;
  ARGNUM = STARTARGNUM;
  UNREADOPTIONS = true;
  while (UNREADOPTIONS) {
    if (ARGNUM >= P_argc) {
      UNREADOPTIONS = false;
      break;
    }  /* if ARGNUM <= argc */
    P_sun_argv(ARGUMENT, 132, (int)ARGNUM);
    if (ARGUMENT[0] == '-') {
      if (ARGUMENT[1] == 'n' || ARGUMENT[1] == 'i') {
	POS = 3;
	switch (ARGUMENT[1]) {

	case 'i':
	  INVERT = true;
	  break;

	case 'n':
	  INTERLEAVED = false;
	  break;
	}
      }  /* i,n */
    }  /* '-' */
    else
      UNREADOPTIONS = false;
    ARGNUM++;
  }

  TEMP.f = stdin;
  *TEMP.name = '\0';
  /* Read in names and sequence */
  NAMESANDSEQ(&TEMP, NAME, SEQ, &SEQLEN, INTERLEAVED, &OKAY);

  if (!OKAY) {
    exit(EXIT_SUCCESS);
  }  /* if OKAY */


  /* If -i, then invert ie. 1->0 and 0->1 */
  if (INVERT)
    INVERTSEQ(SEQ);

  TEMP.f = stdout;
  *TEMP.name = '\0';
  /* Write names and sequence to output file */
  WRITEFLAT(&TEMP, NUMSEQ, SEQLEN, NAME, SEQ);
  exit(EXIT_SUCCESS);
}  /* PHYL2FLAT */



/* End. */
