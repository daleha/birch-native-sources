/* Output from p2c 1.21alpha-07.Dec.93, the Pascal-to-C translator */
/* From input file "getloc.p" */


/* Revision History:
     25 Jan 02 Removed support for LiMB, Vecbase
               PIR files (-p) will be written with the .pir extension
     30 May 95 CRUNCHOFFSET changed to 33
     29 May 95 Copied FINDDATA from GETOB. Now, NAMEFILE does not have
               be sorted, although sorted files will be searched faster.
     21 Dec 94 restores leading blanks to lines for which they have been
               crunched by splitdb -c */
/**************************************************************/
/*                                                            */
/*  GETLOC    VERSION   2/25/02  Standard Pascal              */
/*            Brian Fristensky                                */
/*            Dept. of Plant Science                          */
/*            University of Manitoba                          */
/*            Winnipeg, MB R3T 2N2  CANADA                    */
/*                                                            */
/* SYNOPSIS                                                   */
/* getloc -a -s -f -c namefile anofile seqfile indfile outfile*/
/*                                                            */
/* DESCRIPTION                                                */
/*  Gets GENBANK files from database and writes to OUTFILE    */
/*                                                            */
/*    -a        write annotation data only                    */
/*    -s        write sequence data only                      */
/*    -f        write each entry to a separate file.          */
/*              the filename consists of the locus name,      */
/*              followed by a 3-letter file extesion:         */
/*                   .gen (default, complete genbank entry    */
/*                   .ano (annotation only)                   */
/*                  .wrp (sequence only)                      */
/*    -c        namefile contains accession numbers           */
/*    -database    one of:                                    */
/*                 g  - GenBank (default)                     */
/*                 p  - PIR (NBRF)                            */
/*                 e  - EMBL                                  */
/*                                                            */
/*    NAMEFILE- contains names of entries to get              */
/*    ANOFILE - contains annotation parts of entries          */
/*    SEQFILE - contains sequence parts of entries            */
/*    INDFILE - contains linenumbers for beginning of each    */
/*              entry in ANOFILE and SEQFILE                  */
/*                                                            */
/*  Copyright (c)  1990-2002 by Brian Fristensky              */
/*  !!! in comment indicates feature which may need change    */
/**************************************************************/

#include <p2c.h>


/*!!!  Some Pascals require file parameters in program heading */

#define MAXWORD         25
#define MAXLINE         80
#define CRUNCHOFFSET    33   /* ASCII value of CRUNCHFLAG */
/* BEGIN MODULE STARTARGNUM */
#define STARTARGNUM     1
    /* SUN Pascal: ARG(1) is 1st command line argument*/


/*STARTARGNUM=2;       ATT Pascal: ARG(2) is 1st command line argument*/
/* END MODULE STARTARGNUM         VERSION= 'SUNMODS     Version  8/ 9/94'; */

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

typedef enum {
  GB, PIR, EMBL
} DBTYPE;

typedef struct LOC {
  WORD NAME;
  long ANOLINE, SEQLINE;
} LOC;


Static _TEXT NAMEFILE, ANOFILE, SEQFILE, INDFILE, OUTFILE;
Static WORD SEQNAME, DUMMY;
Static Char ARGUMENT[132];
Static LOC LOCATIONS;   /* locations of a locus in .ANO and .SEQ files */
Static DBTYPE DATABASE;
Static boolean PRINTSEQ, PRINTANO, FILES, ACNO, FOUND;
Static long ARGNUM, CA, CS;
    /* current lines in annotation and sequence files*/
Static long I;
Static LINE FILENAME;
Static Char CRUNCHFLAG;


/***************************************************************/
/* Read options from command line.                             */
/***************************************************************/
Static Void READOPTIONS(ARGNUM, PRINTANO, PRINTSEQ, FILES, ACNO, DATABASE)
long *ARGNUM;
boolean *PRINTANO, *PRINTSEQ, *FILES, *ACNO;
DBTYPE *DATABASE;
{
  Char ARGUMENT[132];

  /* Set defaults */
  *PRINTSEQ = true;
  *PRINTANO = true;
  *FILES = false;
  *ACNO = false;
  *DATABASE = GB;

  /* Read options.*/
  *ARGNUM = STARTARGNUM;
  do {
    P_sun_argv(ARGUMENT, 132, (int)(*ARGNUM));
    if (ARGUMENT[0] == '-') {
      if (ARGUMENT[1] == 'v' || ARGUMENT[1] == 'l' || ARGUMENT[1] == 'e' ||
	  ARGUMENT[1] == 'p' || ARGUMENT[1] == 'g' || ARGUMENT[1] == 'c' ||
	  ARGUMENT[1] == 'f' || ARGUMENT[1] == 's' || ARGUMENT[1] == 'a') {
	switch (ARGUMENT[1]) {

	case 'a':
	  *PRINTSEQ = false;
	  break;

	case 's':
	  *PRINTANO = false;
	  break;

	case 'f':
	  *FILES = true;
	  break;

	case 'c':
	  *ACNO = true;
	  break;

	case 'g':
	  *DATABASE = GB;
	  break;

	case 'p':
	  *DATABASE = PIR;
	  break;

	case 'e':
	  *DATABASE = EMBL;
	  break;
	}
      }
      (*ARGNUM)++;
    }
  } while (ARGUMENT[0] == '-');   /* READOPTIONS */
}


Local Void FILEARGS(F, FTYPE, ARGNUM)
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


/***************************************************************/
/* Open files.                                                 */
/***************************************************************/
Static Void OPENFILES(NAMEFILE, ANOFILE, SEQFILE, INDFILE, OUTFILE)
_TEXT *NAMEFILE, *ANOFILE, *SEQFILE, *INDFILE, *OUTFILE;
{

  /* BEGIN MODULE FILEARGS */
  /* This procedure overcomes one of the stupidest aspects of UNIX Pascal,
     namely the fact that filenames in the program statement are supposed to
     be actual UNIX filenames!  To overcome this, the 2-argument version of
     reset and rewrite must be used with string variables.  This module
     need only contain the reset and rewrite statements in any normal
     implementation of Pascal. */
  /* END MODULE FILEARGS         VERSION= 'SUNMODS     Version  8/ 9/94'; */

  FILEARGS(NAMEFILE, 'I', &ARGNUM);
  if (PRINTANO)
    FILEARGS(ANOFILE, 'I', &ARGNUM);
  if (PRINTSEQ)
    FILEARGS(SEQFILE, 'I', &ARGNUM);
  FILEARGS(INDFILE, 'I', &ARGNUM);
  if (!FILES)
    FILEARGS(OUTFILE, 'O', &ARGNUM);
}  /* OPENFILES */


/***************************************************************/
/* I/O ROUTINES                                                */
/***************************************************************/

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

/**************************************************************/
/* Read a LOCUS name for next sequence to pull from database. */
/**************************************************************/
Static Void READNAME(NAMEFILE, SEQNAME)
_TEXT *NAMEFILE;
WORD *SEQNAME;
{
  /* Skip comment lines */
  while ((P_peek(NAMEFILE->f) == ';') & (!BUFEOF(NAMEFILE->f))) {
    fscanf(NAMEFILE->f, "%*[^\n]");
    getc(NAMEFILE->f);
  }

  /* Read the first token in the next non-comment line */
  SEQNAME->LEN = 0;
  if (!BUFEOF(NAMEFILE->f))
    READWORD(NAMEFILE, SEQNAME);
  if (!BUFEOF(NAMEFILE->f)) {
    fscanf(NAMEFILE->f, "%*[^\n]");
    getc(NAMEFILE->f);
  }
}  /* READNAME */


/* Local variables for FINDDATA: */
struct LOC_FINDDATA {
  _TEXT *INDFILE;
  WORD ID;
  LOC *LOCATIONS;
  boolean *FOUND;
  Char FLAG;
  WORD DUMMY;
} ;

/* BEGIN MODULE SAMEWORD */
/* Compare two WORDS for equality */
Local boolean SAMEWORD(W1, W2, LINK)
WORD *W1, *W2;
struct LOC_FINDDATA *LINK;
{
  long I;
  boolean T;

  if (W1->LEN == W2->LEN) {
    T = true;
    I = 1;
    while (I <= W1->LEN && T) {
      if (W1->STR[I-1] == W2->STR[I-1])
	I++;
      else
	T = false;
    }
    return T;
  } else
    return false;
}  /* SAMEWORD */

/* END MODULE SAMEWORD         VERSION= 'SUNMODS     Version  8/ 9/94'; */

/* BEGIN MODULE TOUPPER */
/* Change a character from lower to uppercase */
Local Char TOUPPER(CH, LINK)
Char CH;
struct LOC_FINDDATA *LINK;
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

/* END MODULE TOUPPER         VERSION= 'SUNMODS     Version  8/ 9/94'; */

/* Advance to first line with first char of ID */
Local Void ZOOM(LINK)
struct LOC_FINDDATA *LINK;
{
  while ((!BUFEOF(LINK->INDFILE->f)) & (P_peek(LINK->INDFILE->f) != LINK->FLAG)) {
    fscanf(LINK->INDFILE->f, "%*[^\n]");
    getc(LINK->INDFILE->f);
  }
}  /* ZOOM */

/* Search for LOCUS corresponding to ID */
Local Void LCRAWL(LINK)
struct LOC_FINDDATA *LINK;
{
  LOC *WITH;

  while (P_peek(LINK->INDFILE->f) == LINK->FLAG && !*LINK->FOUND) {
    WITH = LINK->LOCATIONS;
    READWORD(LINK->INDFILE, &WITH->NAME);
    if (SAMEWORD(&LINK->ID, &WITH->NAME, LINK)) {
      *LINK->FOUND = true;
      READWORD(LINK->INDFILE, &LINK->DUMMY);   /* ignore ACCESSION number */
      fscanf(LINK->INDFILE->f, "%ld%ld", &WITH->ANOLINE, &WITH->SEQLINE);
    }
    if (!BUFEOF(LINK->INDFILE->f)) {
      fscanf(LINK->INDFILE->f, "%*[^\n]");
      getc(LINK->INDFILE->f);
    }
  }
}  /* LCRAWL */

/* Search for ACCESSION number corresponding to ID */
Local Void ACRAWL(LINK)
struct LOC_FINDDATA *LINK;
{
  LOC *WITH;

  while ((!*LINK->FOUND) & (!BUFEOF(LINK->INDFILE->f))) {
    WITH = LINK->LOCATIONS;
    READWORD(LINK->INDFILE, &LINK->DUMMY);   /* skip over LOCUS name */
    READWORD(LINK->INDFILE, &WITH->NAME);
    if (SAMEWORD(&LINK->ID, &WITH->NAME, LINK)) {
      *LINK->FOUND = true;
      fscanf(LINK->INDFILE->f, "%ld%ld", &WITH->ANOLINE, &WITH->SEQLINE);
    }
    if (!BUFEOF(LINK->INDFILE->f)) {
      fscanf(LINK->INDFILE->f, "%*[^\n]");
      getc(LINK->INDFILE->f);
    }
  }
}  /* ACRAWL */


/**********************************************************************/
/* Look for a given LOCUS in the index file.  Since names in NAMEFILE */
/* are supposed to be sorted, reaching the end of file means that the */
/* locus in question is not in the database.                          */
/**********************************************************************/
Static Void FINDDATA(INDFILE_, ID_, LOCATIONS_, FOUND_)
_TEXT *INDFILE_;
WORD ID_;
LOC *LOCATIONS_;
boolean *FOUND_;
{
  struct LOC_FINDDATA V;
  long FORLIM;

  V.INDFILE = INDFILE_;
  V.ID = ID_;
  V.LOCATIONS = LOCATIONS_;
  V.FOUND = FOUND_;
  if (BUFEOF(V.INDFILE->f)) {
    if (*V.INDFILE->name != '\0') {
      if (V.INDFILE->f != NULL)
	V.INDFILE->f = freopen(V.INDFILE->name, "r", V.INDFILE->f);
      else
	V.INDFILE->f = fopen(V.INDFILE->name, "r");
    } else
      rewind(V.INDFILE->f);
    if (V.INDFILE->f == NULL)
      _EscIO2(FileNotFound, V.INDFILE->name);
    RESETBUF(V.INDFILE->f, Char);
  }
  *V.FOUND = false;
  FORLIM = V.ID.LEN;
  /* Make sure that the name or accession number is capitalized. */
  for (I = 1; I <= FORLIM; I++)
    V.ID.STR[I-1] = TOUPPER(V.ID.STR[I-1], &V);
  V.FLAG = V.ID.STR[0];

  /* Search to end of file or until FOUND */
  if (ACNO)
    ACRAWL(&V);
  else {
    while ((!*V.FOUND) & (!BUFEOF(V.INDFILE->f))) {
      ZOOM(&V);
      LCRAWL(&V);
    }
  }

  /* If list isn't sorted, go back to the beginning and search from top */
  if (*V.FOUND)
    return;
  if (*V.INDFILE->name != '\0') {
    if (V.INDFILE->f != NULL)
      V.INDFILE->f = freopen(V.INDFILE->name, "r", V.INDFILE->f);
    else
      V.INDFILE->f = fopen(V.INDFILE->name, "r");
  } else
    rewind(V.INDFILE->f);
  if (V.INDFILE->f == NULL)
    _EscIO2(FileNotFound, V.INDFILE->name);
  RESETBUF(V.INDFILE->f, Char);
  if (ACNO) {
    ACRAWL(&V);
    return;
  }
  while ((!*V.FOUND) & (!BUFEOF(V.INDFILE->f))) {
    ZOOM(&V);
    LCRAWL(&V);
  }
}  /* FINDDATA */


/*****************************************************/
/* Create a filename using the locus.                */
/*****************************************************/
Static Void MAKEFN(NAME, FILENAME)
WORD NAME;
LINE *FILENAME;
{
  long I = 1;

  while (I <= NAME.LEN) {
    FILENAME->STR[I-1] = NAME.STR[I-1];
    I++;
  }
  FILENAME->STR[I-1] = '.';
  if (!PRINTSEQ && PRINTANO) {
    FILENAME->STR[I] = 'a';
    FILENAME->STR[I+1] = 'n';
    FILENAME->STR[I+2] = 'o';
  } else if (!PRINTANO && PRINTSEQ) {
    FILENAME->STR[I] = 'w';
    FILENAME->STR[I+1] = 'r';
    FILENAME->STR[I+2] = 'p';
  } else {
    if (DATABASE == PIR) {
      FILENAME->STR[I] = 'p';
      FILENAME->STR[I+1] = 'i';
      FILENAME->STR[I+2] = 'r';
    } else {
      FILENAME->STR[I] = 'g';
      FILENAME->STR[I+1] = 'e';
      FILENAME->STR[I+2] = 'n';
    }
  }
  FILENAME->LEN = I + 3;
  for (I = FILENAME->LEN; I < MAXLINE; I++)
    FILENAME->STR[I] = ' ';
}  /* MAKEFN */


/*********************************************************/
/* Read annoation data from ANOFILE and write to OUTFILE.*/
/*********************************************************/
Static Void GETANO(ANOLINE, CA)
long ANOLINE, *CA;
{
  Char CH;

  /* Advance to ANOLINE, the beginning of the entry.*/
  if (*CA > ANOLINE) {
    if (*ANOFILE.name != '\0') {
      if (ANOFILE.f != NULL)
	ANOFILE.f = freopen(ANOFILE.name, "r", ANOFILE.f);
      else
	ANOFILE.f = fopen(ANOFILE.name, "r");
    } else
      rewind(ANOFILE.f);
    if (ANOFILE.f == NULL)
      _EscIO2(FileNotFound, ANOFILE.name);
    RESETBUF(ANOFILE.f, Char);
    *CA = 1;
  }
  while (*CA < ANOLINE) {
    fscanf(ANOFILE.f, "%*[^\n]");
    getc(ANOFILE.f);
    (*CA)++;
  }

  /* Transfer data until the end of entry marker is reached.*/
  while (P_peek(ANOFILE.f) != '/') {
    /* If the second character holds a compression value. Expand
       the leading blanks so that there are
       ord(CH)-CRUNCHOFFSET leading blanks. */
    if (P_peek(ANOFILE.f) == CRUNCHFLAG) {
      CH = getc(ANOFILE.f);
      if (CH == '\n')
	CH = ' ';
      CH = getc(ANOFILE.f);
      if (CH == '\n')
	CH = ' ';
      fprintf(OUTFILE.f, "%*c", CH - CRUNCHOFFSET, ' ');
    }

    while (!P_eoln(ANOFILE.f)) {
      CH = getc(ANOFILE.f);
      if (CH == '\n')
	CH = ' ';
      putc(CH, OUTFILE.f);
    }  /* not eoln(ANOFILE) */
    fscanf(ANOFILE.f, "%*[^\n]");
    getc(ANOFILE.f);
    (*CA)++;
    putc('\n', OUTFILE.f);
  }  /* not(ANOFILE^='/') */


  /* if -a option is specified, divide entries with a line for easy reading. */
  if (!PRINTSEQ)
    fprintf(OUTFILE.f,
      "//----------------------------------------------------------------------\n");
  fscanf(ANOFILE.f, "%*[^\n]");
  getc(ANOFILE.f);
  (*CA)++;
}  /* GETANO */


/********************************************************/
/* Read sequence data from SEQFILE and write to OUTFILE.*/
/********************************************************/
Static Void GETSEQ(SEQLINE, CS)
long SEQLINE, *CS;
{
  Char CH;
  boolean MORESEQ;   /* =true if more sequence remains to be read. */
  long SEQLEN;   /* length of sequence printed so far */
  long NUMWIDTH;   /* width of numbers printed to left of sequence */
  long SPACERWIDTH;   /* width of white space between number and sequence */
  long GROUP, LINESIZE;
  /* print <= LINESIZE nucleotides per line in groups
                     of GROUP */
  long I;   /* width of current group */
  long WIDTH;   /* width of current line */
  LINE NAMELINE;

  /* Advance to SEQLINE, the beginning of the entry.*/
  if (*CS > SEQLINE) {
    if (*SEQFILE.name != '\0') {
      if (SEQFILE.f != NULL)
	SEQFILE.f = freopen(SEQFILE.name, "r", SEQFILE.f);
      else
	SEQFILE.f = fopen(SEQFILE.name, "r");
    } else
      rewind(SEQFILE.f);
    if (SEQFILE.f == NULL)
      _EscIO2(FileNotFound, SEQFILE.name);
    RESETBUF(SEQFILE.f, Char);
    *CS = 1;
  }
  while (*CS < SEQLINE) {
    fscanf(SEQFILE.f, "%*[^\n]");
    getc(SEQFILE.f);
    (*CS)++;
  }
  /* Read past name line */
  READLINE(&SEQFILE, &NAMELINE);
  (*CS)++;

  /* Transfer data until the next sequence or end of file is reached.*/
  if (PRINTANO) {   /* write in database-specific format */
    switch (DATABASE) {

    case GB:
    case EMBL:
      /* Set printing parameters according to DATABASE */
      switch (DATABASE) {

      case GB:
	NUMWIDTH = 9;
	SPACERWIDTH = 1;
	LINESIZE = 60;
	break;

      case EMBL:
	SPACERWIDTH = 5;
	LINESIZE = 60;
	break;
      }
      GROUP = 10;
      SEQLEN = 1;

      if ((!BUFEOF(SEQFILE.f)) & (P_peek(SEQFILE.f) != '>')) {
	MORESEQ = true;
	if (DATABASE != EMBL)
	  fprintf(OUTFILE.f, "%*ld", (int)NUMWIDTH, SEQLEN);
	fprintf(OUTFILE.f, "%*c", (int)SPACERWIDTH, ' ');
      } else
	MORESEQ = false;

      /* For each character in the sequence: */
      WIDTH = 0;
      I = 0;
      while (MORESEQ) {
	CH = getc(SEQFILE.f);
	if (CH == '\n')
	  CH = ' ';
	putc(CH, OUTFILE.f);
	if (P_eoln(SEQFILE.f) & (!BUFEOF(SEQFILE.f))) {
	  fscanf(SEQFILE.f, "%*[^\n]");
	  getc(SEQFILE.f);
	  (*CS)++;
	}
	if (BUFEOF(SEQFILE.f) | (P_peek(SEQFILE.f) == '>'))
	  MORESEQ = false;

	if (!MORESEQ) {  /*write spaces or /newline, if necessary*/
	  /* WIDTH=LINESIZE */
	  putc('\n', OUTFILE.f);
	  break;
	}  /* if MORESEQ */
	I++;
	WIDTH++;
	if (I != GROUP) {
	  continue;
	}  /* I=GROUP*/
	I = 0;
	if (WIDTH < LINESIZE) {
	  putc(' ', OUTFILE.f);
	  continue;
	}
	putc('\n', OUTFILE.f);
	SEQLEN += WIDTH;
	WIDTH = 0;
	if (DATABASE != EMBL)
	  fprintf(OUTFILE.f, "%*ld", (int)NUMWIDTH, SEQLEN);
	fprintf(OUTFILE.f, "%*c", (int)SPACERWIDTH, ' ');
      }  /* while MORESEQ */
      /* Read a character, write to OUTFILE, and determine if there
         if more sequence left to read. */
      fprintf(OUTFILE.f, "//\n");
      break;
      /* GB,EMBL,VECTOR */

    case PIR:
      fprintf(OUTFILE.f, "%7c", ' ');
      I = 5;
      while (I <= 30) {
	fprintf(OUTFILE.f, "%10ld", I);
	I += 5;
      }
      putc('\n', OUTFILE.f);

      /* write sequence in PIR format */
      /* punctuation characters occupy odd numbered columnes 1..61
         and amino acids occupy even numbered columns 2..60. This
         algorithm assumes that punctuation characters may not occur
         adjacent to one another. */
      NUMWIDTH = 7;
      LINESIZE = 60;
      SPACERWIDTH = 1;
      SEQLEN = 1;
      WIDTH = 0;

      if ((!BUFEOF(SEQFILE.f)) & (P_peek(SEQFILE.f) != '>')) {
	MORESEQ = true;
	fprintf(OUTFILE.f, "%*ld", (int)NUMWIDTH, SEQLEN);
      } else
	MORESEQ = false;

      while (MORESEQ) {
	CH = getc(SEQFILE.f);
	if (CH == '\n')
	  CH = ' ';
	if (P_eoln(SEQFILE.f) & (!BUFEOF(SEQFILE.f))) {
	  fscanf(SEQFILE.f, "%*[^\n]");
	  getc(SEQFILE.f);
	  (*CS)++;
	}
	if (BUFEOF(SEQFILE.f) | (P_peek(SEQFILE.f) == '>'))
	  MORESEQ = false;
	if (CH == ',' || CH == '.' || CH == '/' || CH == '=' || CH == ')' ||
	    CH == '(') {
	  putc(CH, OUTFILE.f);
	  WIDTH++;
	} else {
	  if (!(WIDTH & 1)) {
	    putc(' ', OUTFILE.f);
	    WIDTH++;
	  }
	  putc(CH, OUTFILE.f);
	  WIDTH++;
	}
	if (!MORESEQ) {
	  putc('\n', OUTFILE.f);
	  break;
	}
	if (WIDTH >= LINESIZE) {
	  putc('\n', OUTFILE.f);
	  SEQLEN += 30;
	  WIDTH = 0;
	  fprintf(OUTFILE.f, "%*ld", (int)NUMWIDTH, SEQLEN);
	}  /* WIDTH>=LINESIZE */
      }  /* while MORESEQ */
      fprintf(OUTFILE.f, "///\n");
      break;
      /* write number header */
      /*PIR*/
    }/* case DATABASE */

    return;
  }
  WRITELINE(&OUTFILE, NAMELINE, NAMELINE.LEN);
  putc('\n', OUTFILE.f);
  while (!((P_peek(SEQFILE.f) == '>') | BUFEOF(SEQFILE.f))) {
    while (!P_eoln(SEQFILE.f)) {
      CH = getc(SEQFILE.f);
      if (CH == '\n')
	CH = ' ';
      putc(CH, OUTFILE.f);
    }
    fscanf(SEQFILE.f, "%*[^\n]");
    getc(SEQFILE.f);
    (*CS)++;
    putc('\n', OUTFILE.f);
    /* write in Pearson format */
  }
}  /* GETSEQ */


/* ----------------------------------------------------------  */
/* ----------------- MAIN  PROCEDURE  -----------------------  */
main(argc, argv)
int argc;
Char *argv[];
{
  PASCAL_MAIN(argc, argv);
  OUTFILE.f = NULL;
  strcpy(OUTFILE.name, "OUTFILE");
  INDFILE.f = NULL;
  strcpy(INDFILE.name, "INDFILE");
  SEQFILE.f = NULL;
  strcpy(SEQFILE.name, "SEQFILE");
  ANOFILE.f = NULL;
  strcpy(ANOFILE.name, "ANOFILE");
  NAMEFILE.f = NULL;
  strcpy(NAMEFILE.name, "NAMEFILE");
  /* Read options from the command line */
  READOPTIONS(&ARGNUM, &PRINTANO, &PRINTSEQ, &FILES, &ACNO, &DATABASE);

  /* Open files.*/
  OPENFILES(&NAMEFILE, &ANOFILE, &SEQFILE, &INDFILE, &OUTFILE);


  /* For each locus in NAMEFILE, find the locations of the annotation and
     sequence parts in ANOFILE and SEQFILE, and then get the data
     and write to OUTFILE in GenBank tape format. */
  CRUNCHFLAG = (Char)CRUNCHOFFSET;
  CA = 1;
  CS = 1;
  while (!BUFEOF(NAMEFILE.f)) {
    READNAME(&NAMEFILE, &SEQNAME);
    if (SEQNAME.LEN > 0)
      FINDDATA(&INDFILE, SEQNAME, &LOCATIONS, &FOUND);
    if (!FOUND)
      continue;
    if (FILES) {
      MAKEFN(SEQNAME, &FILENAME);
      strcpy(OUTFILE.name, P_trimname(FILENAME.STR, sizeof(CHARARRAY)));
      if (OUTFILE.f != NULL)
	OUTFILE.f = freopen(OUTFILE.name, "w", OUTFILE.f);
      else
	OUTFILE.f = fopen(OUTFILE.name, "w");
      if (OUTFILE.f == NULL)
	_EscIO2(FileNotFound, OUTFILE.name);
      SETUPBUF(OUTFILE.f, Char);
    }
    if (PRINTANO)
      GETANO(LOCATIONS.ANOLINE, &CA);
    if (PRINTSEQ)
      GETSEQ(LOCATIONS.SEQLINE, &CS);
    /*!!!      if FILES then CLOSE(OUTFILE) */
  }  /* while not eof(NAMEFILE) */

  /*!!! if not FILES then CLOSE(OUTFILE); */
  if (NAMEFILE.f != NULL)
    fclose(NAMEFILE.f);
  if (ANOFILE.f != NULL)
    fclose(ANOFILE.f);
  if (SEQFILE.f != NULL)
    fclose(SEQFILE.f);
  if (INDFILE.f != NULL)
    fclose(INDFILE.f);
  if (OUTFILE.f != NULL)
    fclose(OUTFILE.f);
  exit(EXIT_SUCCESS);
}  /* GETLOC  */




/* End. */
