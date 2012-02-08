/* Output from p2c 1.21alpha-07.Dec.93, the Pascal-to-C translator */
/* From input file "flat2phyl.p" */


/* ********************************************************  */
/*                                                           */
/*  FLAT2PHYL  VERSION  5/ 2/2000, Standard Pascal           */
/*             Brian Fristensky                              */
/*             Dept. of Plant Science                        */
/*             University of Manitoba                        */
/*             Winnipeg, MB  R3T 2N2 CANADA                  */
/*                                                           */
/* SYNOPSIS                                                  */
/*    flat2phyl [-i] < GDE_flatfile > Phylip_file            */
/*                                                           */
/* DESCRIPTION                                               */
/*    Convert a GDE flatfile into a                          */
/*    PHYLIP discrete character file.                        */
/*                                                           */
/*             -s  write in Phylip sequential format.        */
/*                 rather than the default, which is         */
/*                 Phylip interleaved format.                */
/*             -i  invert characters, so that 0 -> 1 and     */
/*                 1 -> 0.                                   */
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

#define VERSION         "FLAT2PHYL     Version 5/ 2/2000"


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
Static boolean SEQUENTIAL, INVERT, UNREADOPTIONS;
Static Char ARGUMENT[132];   /* command line argument */


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


/* Local variables for READFLAT: */
struct LOC_READFLAT {
  long *NUMSEQ;
} ;


/* --------------------------------------------*/
/* Calculate the minimum and maximum sequence lengths read so far. */
Local Void MINMAXREAD(LENGTH, MINREAD, MAXREAD, LINK)
long *LENGTH;
long *MINREAD, *MAXREAD;
struct LOC_READFLAT *LINK;
{
  long I, FORLIM;

  *MINREAD = LENGTH[0];
  *MAXREAD = *MINREAD;
  FORLIM = *LINK->NUMSEQ;
  for (I = 1; I < FORLIM; I++) {
    if (LENGTH[I] < *MINREAD)
      *MINREAD = LENGTH[I];
    if (LENGTH[I] > *MAXREAD)
      *MAXREAD = LENGTH[I];
  }
}  /* MINMAXREAD */


/* END MODULE WRITEWORD         VERSION= 'SUNMODS     Version  8/ 9/94'; */

/***************************************************************/
/* Read aligned sequences into arrays.                         */
/***************************************************************/
Static Void READFLAT(F, NUMSEQ_, NAME, SEQ, SEQLEN)
_TEXT *F;
long *NUMSEQ_;
WORD *NAME;
Char (*SEQ)[MAXCHAR];
long *SEQLEN;
{  /* READFLAT ----------------------------------------------- */
  struct LOC_READFLAT V;
  long I, MINREAD, MAXREAD;
  boolean ENDOFSEQUENCE;
  Char CH;

  V.NUMSEQ = NUMSEQ_;
  /* Read in each sequence. */
  *V.NUMSEQ = 0;
  for (I = 0; I < MAXSEQ; I++)
    LENGTH[I] = 0;

  while (!BUFEOF(F->f)) {
    /* Read a name */
    if (P_peek(F->f) == '"') {
      CH = getc(F->f);   /*read past '"' */
      if (CH == '\n')
	CH = ' ';
      (*V.NUMSEQ)++;
      READWORD(F, &NAME[*V.NUMSEQ - 1]);   /* read name */
      if (!BUFEOF(F->f)) {
	fscanf(F->f, "%*[^\n]");
	getc(F->f);
      }
    }

    /* Read sequence */
    ENDOFSEQUENCE = false;
    while (!ENDOFSEQUENCE) {
      while (!P_eoln(F->f)) {
	CH = getc(F->f);
	if (CH == '\n')
	  CH = ' ';
	if (CH == 'B' || CH == 'P' || CH == '?' || CH == '-' || CH == '+' ||
	    CH == '1' || CH == '0') {
	  LENGTH[*V.NUMSEQ - 1]++;
	  SEQ[*V.NUMSEQ - 1][LENGTH[*V.NUMSEQ - 1] - 1] = CH;
	}
      }  /* while not eoln(F) */
      if (!BUFEOF(F->f)) {
	fscanf(F->f, "%*[^\n]");
	getc(F->f);
      }
      if (BUFEOF(F->f))
	ENDOFSEQUENCE = true;
      else if (P_peek(F->f) == '"')
	ENDOFSEQUENCE = true;
    }  /* not ENDOFSEQUENCE */
  }  /* not eof(F) */

  MINMAXREAD(LENGTH, &MINREAD, &MAXREAD, &V);
  *SEQLEN = MINREAD;
}  /* READFLAT */


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
/* Write names and sequences in Phylip interleaved format.             */
/* ******************************************************************* */
Static Void WRITEPHYL(F, NUMSEQ, SEQLEN, NAME, SEQ, SEQUENTIAL)
_TEXT *F;
long NUMSEQ, SEQLEN;
WORD *NAME;
Char (*SEQ)[MAXCHAR];
boolean SEQUENTIAL;
{
  long I, K, THISLINE, SEQPRINTED;
  boolean FIRSTLINE = true;

  fprintf(F->f, "%10ld%10ld\n", NUMSEQ, SEQLEN);

  if (SEQUENTIAL) {  /* Phylip Sequential Format*/
    for (I = 0; I < NUMSEQ; I++) {
      /* some PHYLIP programs want 10 letter names,
         with at least 1 char of padding at the end. */
      WRITEWORD(F, NAME[I], 11L);
      putc('\n', F->f);
      SEQPRINTED = 0;
      while (SEQPRINTED < SEQLEN) {
	/* Write a line of sequence data */
	K = 1;
	THISLINE = SEQPRINTED + LINELEN;
	if (THISLINE > SEQLEN)
	  THISLINE = SEQLEN;
	for (K = SEQPRINTED; K < THISLINE; K++)
	  putc(SEQ[I][K], F->f);
	putc('\n', F->f);
	SEQPRINTED = THISLINE;
      }
    }

    return;
  }  /* Phylip Sequential Format*/

  SEQPRINTED = 0;
  while (SEQPRINTED < SEQLEN) {
    for (I = 0; I < NUMSEQ; I++) {
      /* Write out names on the first set of lines */
      if (FIRSTLINE)
	WRITEWORD(F, NAME[I], 10L);

      /* Write a line of sequence data */
      K = 1;
      THISLINE = SEQPRINTED + LINELEN;
      if (THISLINE > SEQLEN)
	THISLINE = SEQLEN;
      for (K = SEQPRINTED; K < THISLINE; K++)
	putc(SEQ[I][K], F->f);
      putc('\n', F->f);
    }  /* for I:= 1 to NUMSEQ */

    FIRSTLINE = false;
    SEQPRINTED = THISLINE;
    /* Phylip Interleaved Format*/
    /* Phylip Interleaved Format*/
  }  /* while SEQPRINTED < SEQLEN */
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
  SEQUENTIAL = false;
  ARGNUM = STARTARGNUM;
  UNREADOPTIONS = true;
  while (UNREADOPTIONS) {
    if (ARGNUM >= P_argc) {
      UNREADOPTIONS = false;
      break;
    }  /* if ARGNUM <= argc */
    P_sun_argv(ARGUMENT, 132, (int)ARGNUM);
    if (ARGUMENT[0] == '-') {
      if (ARGUMENT[1] == 's' || ARGUMENT[1] == 'i') {
	POS = 3;
	switch (ARGUMENT[1]) {

	case 'i':
	  INVERT = true;
	  break;

	case 's':
	  SEQUENTIAL = true;
	  break;
	}
      }  /* i */
    }  /* '-' */
    else
      UNREADOPTIONS = false;
    ARGNUM++;
  }

  TEMP.f = stdin;
  *TEMP.name = '\0';
  /* Read in names and sequence */
  READFLAT(&TEMP, &NUMSEQ, NAME, SEQ, &SEQLEN);

  /* If -i, then invert ie. 1->0 and 0->1 */
  if (INVERT)
    INVERTSEQ(SEQ);

  TEMP.f = stdout;
  *TEMP.name = '\0';
  /* Write names and sequence to output file */
  WRITEPHYL(&TEMP, NUMSEQ, SEQLEN, NAME, SEQ, SEQUENTIAL);

  exit(EXIT_SUCCESS);
}  /* FLAT2PHYL */



/* End. */
