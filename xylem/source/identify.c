/* Output from p2c 1.21alpha-07.Dec.93, the Pascal-to-C translator */
/* From input file "identify.p" */


/***********************************************************/
/*                                                         */
/*  IDENTIFY  VERSION   8/19/90  UNIX Pascal               */
/*            Brian Fristensky                             */
/*            Dept. of Plant Science                       */
/*            University of Manitoba                       */
/*            Winnipeg, MB  R3T 2N2  CANADA                */
/*                                                         */
/*  Uses GREP output to identify locus names in .IND file  */
/*  that correspond to lines found in .ANO file.           */
/*                                                         */
/*    GREPFILE- contains lines found in .ANO file by GREP  */
/*    INDFILE - contains linenumbers for beginning of each */
/*              entry in ANOFILE and SEQFILE               */
/*    NAMEFILE- output file listing all loci found         */
/*    FINDFILE- output file with locus names and lines     */
/*              found for each locus                       */
/*                                                         */
/*  Copyright (c) l988, 1990 by Brian Fristensky           */
/*  !!! in comment indicates feature which may need change */
/***********************************************************/

#include <p2c.h>


/*!!!  Some Pascals require file parameters in program heading */

#define MAXWORD         25
#define MAXLINE         132
/*BEGIN MODULE STARTARGNUM */
#define STARTARGNUM     1


/*END MODULE STARTARGNUM */

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


Static _TEXT GREPFILE, INDFILE, NAMEFILE, FINDFILE;
Static LINE GREPLINE;
Static WORD CURNAME, NEXTNAME, DUMMY;
Static long ARGNUM, GREPNUM, CURNUM, NEXTNUM;


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


/* END MODULE FILEARGS         VERSION= 'SUNMODS     Version  9/12/91'; */

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


/* END MODULE READLINE         VERSION= 'SUNMODS     Version  9/12/91'; */

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


/* END MODULE WRITEWORD         VERSION= 'SUNMODS     Version  9/12/91'; */

/*********************************************************/
/* Read a line number and line from GREPFILE.            */
/*********************************************************/
Static Void GLINE(GREPFILE, GREPLINE, GREPNUM)
_TEXT *GREPFILE;
LINE *GREPLINE;
long *GREPNUM;
{
  Char CH;
  long ORDZERO = '0';

  *GREPNUM = 0;
  do {
    CH = getc(GREPFILE->f);
    if (CH == '\n')
      CH = ' ';
    if (CH != ':')
      *GREPNUM = *GREPNUM * 10 + CH - ORDZERO;
  } while (CH != ':');
  READLINE(GREPFILE, GREPLINE);
}  /* GLINE */


/*********************************************************/
/* Read INDFILE until CURNUM <= GREPNUM < NEXTNUM.       */
/*********************************************************/
Static Void SCAN(INDFILE, CURNAME, NEXTNAME, CURNUM, NEXTNUM)
_TEXT *INDFILE;
WORD *CURNAME, *NEXTNAME;
long *CURNUM, *NEXTNUM;
{
  WORD DUMMY;

  if (!BUFEOF(INDFILE->f)) {
    while ((*NEXTNUM <= GREPNUM) & (!BUFEOF(INDFILE->f))) {
      *CURNAME = *NEXTNAME;
      *CURNUM = *NEXTNUM;
      READWORD(INDFILE, NEXTNAME);
      READWORD(INDFILE, &DUMMY);   /* ignore ACCESSION number */
      fscanf(INDFILE->f, "%ld%*[^\n]", NEXTNUM);
      getc(INDFILE->f);
    }
    return;
  }
  *CURNAME = *NEXTNAME;
  *CURNUM = *NEXTNUM;
  *NEXTNUM = 10000000L;   /* GREPNUM can never be this large */
  /* Special case for last entry in INDFILE */
}  /* SCAN */


/* ----------------------------------------------------------  */
/* ----------------- MAIN  PROCEDURE  -----------------------  */
main(argc, argv)
int argc;
Char *argv[];
{
  PASCAL_MAIN(argc, argv);
  FINDFILE.f = NULL;
  strcpy(FINDFILE.name, "FINDFILE");
  NAMEFILE.f = NULL;
  strcpy(NAMEFILE.name, "NAMEFILE");
  INDFILE.f = NULL;
  strcpy(INDFILE.name, "INDFILE");
  GREPFILE.f = NULL;
  strcpy(GREPFILE.name, "GREPFILE");
  /* Open files.*/
  ARGNUM = STARTARGNUM;
  FILEARGS(&GREPFILE, 'I', &ARGNUM);
  FILEARGS(&INDFILE, 'I', &ARGNUM);
  FILEARGS(&NAMEFILE, 'O', &ARGNUM);
  FILEARGS(&FINDFILE, 'O', &ARGNUM);

  /* Read past comment lines at beginning of INDFILE, if any.*/
  /* Read the first two loci in INDFILE */
  while ((P_peek(INDFILE.f) == ';') & (!BUFEOF(INDFILE.f))) {
    fscanf(INDFILE.f, "%*[^\n]");
    getc(INDFILE.f);
  }
  READWORD(&INDFILE, &CURNAME);
  READWORD(&INDFILE, &DUMMY);   /* ignore ACCESSION number */
  fscanf(INDFILE.f, "%ld%*[^\n]", &CURNUM);
  getc(INDFILE.f);
  READWORD(&INDFILE, &NEXTNAME);
  READWORD(&INDFILE, &DUMMY);   /* ignore ACCESSION number */
  fscanf(INDFILE.f, "%ld%*[^\n]", &NEXTNUM);
  getc(INDFILE.f);

  /* Read the first line in GREPFILE */
  GLINE(&GREPFILE, &GREPLINE, &GREPNUM);

  /* Find the LOCUS corresponding to GREPNUM */
  SCAN(&INDFILE, &CURNAME, &NEXTNAME, &CURNUM, &NEXTNUM);

  /* Write the first LOCUS name to NAMEFILE and FINDFILE and the
     GREPLINE to FINDFILE */
  WRITEWORD(&FINDFILE, CURNAME, CURNAME.LEN);
  putc('\n', FINDFILE.f);
  WRITEWORD(&NAMEFILE, CURNAME, CURNAME.LEN);
  putc('\n', NAMEFILE.f);
  WRITELINE(&FINDFILE, GREPLINE, GREPLINE.LEN);
  putc('\n', FINDFILE.f);


  /* MAIN LOOP */
  /* Invariant: CURNUM <= GREPNUM < NEXTNUM */
  while (!BUFEOF(GREPFILE.f)) {
    /* Read a line from GREPFILE.  The first characters on the line tell
       the line number, the rest is from the .ANO file. */
    GLINE(&GREPFILE, &GREPLINE, &GREPNUM);

    /* Search INDFILE until the position of the next entry is greater
       than the line number from the GREPFILE.*/
    if (GREPNUM >= NEXTNUM) {
      fprintf(FINDFILE.f, "//\n");
      SCAN(&INDFILE, &CURNAME, &NEXTNAME, &CURNUM, &NEXTNUM);
      WRITEWORD(&FINDFILE, CURNAME, CURNAME.LEN);
      putc('\n', FINDFILE.f);
      WRITEWORD(&NAMEFILE, CURNAME, CURNAME.LEN);
      putc('\n', NAMEFILE.f);
    }

    /* Write the line into FINDFILE */
    WRITELINE(&FINDFILE, GREPLINE, GREPLINE.LEN);
    putc('\n', FINDFILE.f);
  }


  fprintf(FINDFILE.f, "//\n");
  /*!!!  CLOSE(NAMEFILE); */
  /*!!!  CLOSE(FINDFILE)  */
  if (GREPFILE.f != NULL)
    fclose(GREPFILE.f);
  if (INDFILE.f != NULL)
    fclose(INDFILE.f);
  if (NAMEFILE.f != NULL)
    fclose(NAMEFILE.f);
  if (FINDFILE.f != NULL)
    fclose(FINDFILE.f);
  exit(EXIT_SUCCESS);
}  /* IDENTIFY  */



/* End. */
