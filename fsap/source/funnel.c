/* Output from p2c 1.21alpha-07.Dec.93, the Pascal-to-C translator */
/* From input file "funnel.p" */


/*$DEBUG- */
/* ********************************************************* */
/*                                                           */
/* FUNNEL  Version  5/13/91   Standard Pascal                */
/*         Brian Fristensky                                  */
/*         Dept. of Plant Science                            */
/*         University of Manitoba                            */
/*         Winnipeg, MB R3T 2N2  CANADA                      */
/*                                                           */
/* Copyright (c) 1987, 1990 by Brian Fristensky              */
/* !!! in comment indicates feature which may need change    */
/* ********************************************************* */

#include <p2c.h>


/*,INFILE,OUTFILE*/
/*!!! Some Pascals require file parameters in program heading */

#define VERSION         "FUNNEL            Version  5/13/91"

#define MAXWORD         25
#define MAXLINE         10000


/* BEGIN MODULE INTLIMITS */
/*!!!  MAXINT =  2147483647; */
/*!!!  MININT = -2147483647; */
/* END MODULE INTLIMITS         VERSION= 'SUNMODS     Version  9/12/91'; */

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


Static _TEXT INFILE, OUTFILE;
Static LINE IFN, OFN;
Static Char CH;
Static long I, LINESIZE;
Static boolean COMMENT;


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


/* END MODULE GETINTEGER         VERSION= 'SUNMODS     Version  9/12/91'; */

main(argc, argv)
int argc;
Char *argv[];
{
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
  /* END MODULE STARTUP         VERSION= 'SUNMODS     Version  9/12/91'; */
  /*Ask user for INFILE and OUTFILE filenames*/
  printf("INPUT FILENAME?\n");
  GETFILE(&INFILE, 'I', &IFN);
  printf("OUTPUT FILENAME?\n");
  GETFILE(&OUTFILE, 'O', &OFN);
  printf("Enter LINESIZE:\n");
  GETINTEGER(&LINESIZE, 1L, (long)MAXLINE);
  scanf("%*[^\n]");
  getchar();
  /*READ INPUT FILE CHARACTER BY CHARACTER*/
  /*AND WRITE EACH CHARACTER INTO OUTPUT FILE*/
  I = 0;
  COMMENT = false;
  while (!BUFEOF(INFILE.f)) {
    if (COMMENT)
      I = 0;
    COMMENT = false;
    while (!P_eoln(INFILE.f)) {
      CH = getc(INFILE.f);
      if (CH == '\n')
	CH = ' ';
      if (CH == '>' || CH == ';') {
	COMMENT = true;
	if (I > 0)
	  putc('\n', OUTFILE.f);
      }
      switch (COMMENT) {

      case true:
	putc(CH, OUTFILE.f);
	break;

      case false:
	if (CH == '*' || isalpha(CH)) {
	  putc(CH, OUTFILE.f);
	  I++;
	  if (I == LINESIZE) {
	    I = 0;
	    putc('\n', OUTFILE.f);
	  }
	}
	break;
      }
    }
    fscanf(INFILE.f, "%*[^\n]");
    getc(INFILE.f);
    if (COMMENT)
      putc('\n', OUTFILE.f);
  }
  /*!!!*/
  if (OUTFILE.f != NULL)
    fclose(OUTFILE.f);
  OUTFILE.f = NULL;
  if (INFILE.f != NULL)
    fclose(INFILE.f);
  if (OUTFILE.f != NULL)
    fclose(OUTFILE.f);
  exit(EXIT_SUCCESS);
}  /* FUNNEL */




/* End. */
