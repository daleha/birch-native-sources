/* Output from p2c 1.21alpha-07.Dec.93, the Pascal-to-C translator */
/* From input file "splitdb.p" */


/* UPDATE HISTORY
      2 Sep 2001 GenBank has shortened LOCUS names to 16 char.
     15 Aug 2001 Revised to allow LOCUS names up to 18 char., to
     conform with upcomming changes in GenBank LOCUS line. Also
     removed support for the ancient LIMB database.
     28 Mar 98 added -t option to include genus & species in .ind file
     26 Aug 95 Writes ACCESSION number to INDFILE using 8 characters.
     30 May 95 Changed CRUNCHOFFSET to ASCII 33 = "!". VT100 emulators
     didn't like ASCII 29.
     21 Dec 94 -c option compresses leading blanks in .ano file.
     This version of splitdb now looks first for an ACCESSIONS line, and if
     none is found, it looks for the first #accession line.
  */

/***********************************************************/
/*                                                         */
/*  SPLITDB   VERSION    9/ 2/2001  Standard Pascal        */
/*            Brian Fristensky                             */
/*            Dept. of Plant Science                       */
/*            University of Manitoba                       */
/*            Winnipeg, MB Canada R3T 2N2                  */
/*                                                         */
/*  Splits GENBANK file into three files:                  */
/*    ANOFILE - contains annotation parts of entries       */
/*    SEQFILE - contains sequence parts of entries         */
/*    INDFILE - contains linenumbers for beginning of each */
/*              entry in ANOFILE and SEQFILE               */
/*                                                         */
/*           -c compress leading blanks in .ano file       */
/*           -t append the contents of the first ORGANISM  */
/*              line to each line in INDFILE               */
/*                                                         */
/*  Copyright (c) 1988 - 2001      by Brian Fristensky     */
/*  !!! in comment indicates feature which may need change */
/***********************************************************/

#include <p2c.h>


/*!!!  Some Pascals require file parameters in program heading */

#define MAXLINE         132
#define MAXWORD         25
#define SEQLINELEN      75   /* length of sequence line in SEQFILE */
#define CRUNCHOFFSET    33

/* BEGIN MODULE STARTARGNUM */
#define STARTARGNUM     1
    /* SUN Pascal: ARG(1) is 1st command line argument*/


/*STARTARGNUM=2;       ATT Pascal: ARG(2) is 1st command line argument*/
/* END MODULE STARTARGNUM         VERSION= 'SUNMODS     Version  8/ 9/94'; */

typedef enum {
  GB, PIR, EMBL, VECTOR
} DBTYPE;

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


Static _TEXT DBFILE, ANOFILE, SEQFILE, INDFILE;
Static WORD NAMEID, DEFID, ACID, SEQID, SOURCEID, ORGID;
    /* keyword identifiers */
Static WORD OLDACID, REFID;   /*can be eliminated after PIR changes */
Static WORD NAME, ACCESSION;   /* name and accession numbers */
Static LINE ORGANISM;
Static Char CRUNCHFLAG;   /* =chr(CRUNCHOFFSET), indicates -c compression */
Static Char CH;
Static DBTYPE DATABASE;
Static boolean TAXONOMY;   /* =true, include genus & sp. in index file */
Static boolean COMPRESSION;   /* =true, then compress leading blanks */
Static LINE CLINE, TITLE;   /* current line, title line for SEQFILE */
Static long POSITION;   /* next unprocessed position in CLINE */
Static long CANO, CSEQ;   /* current line pointers for ANOFILE & SEQFILE */
Static long FIRSTANO;   /* first annotation line for an entry in ANOFILE*/
Static long FIRSTSEQ;   /* title line for a sequence in SEQFILE */
Static long ARGNUM;   /* number of command line argument */
/* arg # of annotation file, used for rewrite */
Static long ANOARGNUM;

Static Char ANOFILENAME[132];


/***************************************************************/
/* Read options from command line and set DATABASE             */
/***************************************************************/
Static Void READOPTIONS(ARGNUM, DATABASE, COMPRESSION, TAXONOMY)
long *ARGNUM;
DBTYPE *DATABASE;
boolean *COMPRESSION, *TAXONOMY;
{
  Char ARGUMENT[132];
  boolean OPTIONSDONE = false;

  *ARGNUM = STARTARGNUM;
  *DATABASE = GB;
  *TAXONOMY = false;
  *COMPRESSION = false;
  do {
    P_sun_argv(ARGUMENT, 132, (int)(*ARGNUM));
    if (ARGUMENT[0] == '-') {
      if (ARGUMENT[1] == 'c' || ARGUMENT[1] == 't' || ARGUMENT[1] == 'v' ||
	  ARGUMENT[1] == 'e' || ARGUMENT[1] == 'p' || ARGUMENT[1] == 'g') {
	switch (ARGUMENT[1]) {

	case 'g':
	  *DATABASE = GB;
	  break;

	case 'p':
	  *DATABASE = PIR;
	  break;

	case 'e':
	  *DATABASE = EMBL;
	  break;

	case 'v':
	  *DATABASE = VECTOR;
	  break;

	case 't':
	  *TAXONOMY = true;
	  break;

	case 'c':
	  *COMPRESSION = true;
	  break;
	}
      }
      (*ARGNUM)++;
    } else
      OPTIONSDONE = true;
  } while (!OPTIONSDONE);   /* READOPTIONS */
}


/***************************************************************/
/* Initialize keywords. Keywords identify the line on which a  */
/* data item occurs.  NAMEID is the keyword that identifies the*/
/* line on which the name of the entry occurs (eg. 'LOCUS' for */
/* GenBank, 'ENTRY' for PIR and so forth.)  DEFID is the TITLE,*/
/* or DEFINITION line, to be written to SEQFILE.  ACID is the  */
/* ACCESSION number line. SEQID denotes where the sequence     */
/* begins. Not all databases have all types of keywords, and   */
/* the order of occurrence varies with database. SOURCEID and  */
/* ORGID are used for GenBank SOURCE and ORGANISM lines.       */
/***************************************************************/
Static Void INITKEYS(NAMEID, DEFID, ACID, SEQID, SOURCEID, ORGID)
WORD *NAMEID, *DEFID, *ACID, *SEQID, *SOURCEID, *ORGID;
{
  switch (DATABASE) {

  case GB:
  case VECTOR:
    NAMEID->STR[0] = 'L';
    NAMEID->STR[1] = 'O';
    NAMEID->STR[2] = 'C';
    NAMEID->STR[3] = 'U';
    NAMEID->STR[4] = 'S';
    NAMEID->LEN = 5;
    break;
    /* GB,VECTOR */

  case PIR:
    NAMEID->STR[0] = 'E';
    NAMEID->STR[1] = 'N';
    NAMEID->STR[2] = 'T';
    NAMEID->STR[3] = 'R';
    NAMEID->STR[4] = 'Y';
    NAMEID->LEN = 5;
    break;
    /* PIR */

  case EMBL:
    NAMEID->STR[0] = 'I';
    NAMEID->STR[1] = 'D';
    NAMEID->LEN = 2;
    break;
    /* EMBL */
  }

  switch (DATABASE) {

  case GB:
  case VECTOR:
    DEFID->STR[0] = 'D';
    DEFID->STR[1] = 'E';
    DEFID->STR[2] = 'F';
    DEFID->STR[3] = 'I';
    DEFID->STR[4] = 'N';
    DEFID->STR[5] = 'I';
    DEFID->STR[6] = 'T';
    DEFID->STR[7] = 'I';
    DEFID->STR[8] = 'O';
    DEFID->STR[9] = 'N';
    DEFID->LEN = 10;
    break;
    /* GB,VECTOR */

  case PIR:
    DEFID->STR[0] = 'T';
    DEFID->STR[1] = 'I';
    DEFID->STR[2] = 'T';
    DEFID->STR[3] = 'L';
    DEFID->STR[4] = 'E';
    DEFID->LEN = 5;
    break;
    /* PIR */

  case EMBL:
    DEFID->STR[0] = 'D';
    DEFID->STR[1] = 'E';
    DEFID->LEN = 2;
    break;
    /* EMBL */
  }

  switch (DATABASE) {

  case GB:
  case VECTOR:
    ACID->STR[0] = 'A';
    ACID->STR[1] = 'C';
    ACID->STR[2] = 'C';
    ACID->STR[3] = 'E';
    ACID->STR[4] = 'S';
    ACID->STR[5] = 'S';
    ACID->STR[6] = 'I';
    ACID->STR[7] = 'O';
    ACID->STR[8] = 'N';
    ACID->LEN = 9;
    break;
    /* GB,VECTOR */

  case PIR:
    ACID->STR[0] = '#';
    ACID->STR[1] = 'a';
    ACID->STR[2] = 'c';
    ACID->STR[3] = 'c';
    ACID->STR[4] = 'e';
    ACID->STR[5] = 's';
    ACID->STR[6] = 's';
    ACID->STR[7] = 'i';
    ACID->STR[8] = 'o';
    ACID->STR[9] = 'n';
    ACID->LEN = 10;
    break;
    /* PIR */

  case EMBL:
    ACID->STR[0] = 'A';
    ACID->STR[1] = 'C';
    ACID->LEN = 2;
    break;
    /* EMBL */
  }

  switch (DATABASE) {

  case GB:
  case VECTOR:
    SEQID->STR[0] = 'O';
    SEQID->STR[1] = 'R';
    SEQID->STR[2] = 'I';
    SEQID->STR[3] = 'G';
    SEQID->STR[4] = 'I';
    SEQID->STR[5] = 'N';
    SEQID->LEN = 6;
    break;
    /* GB,VECTOR */

  case PIR:
    SEQID->STR[0] = 'S';
    SEQID->STR[1] = 'E';
    SEQID->STR[2] = 'Q';
    SEQID->STR[3] = 'U';
    SEQID->STR[4] = 'E';
    SEQID->STR[5] = 'N';
    SEQID->STR[6] = 'C';
    SEQID->STR[7] = 'E';
    SEQID->LEN = 8;
    break;
    /* PIR */

  case EMBL:
    SEQID->STR[0] = 'S';
    SEQID->STR[1] = 'Q';
    SEQID->LEN = 2;
    break;
    /* EMBL */
  }
  /* Used for GenBank SOURCE and ORGANISM lines. */
  SOURCEID->STR[0] = 'S';
  SOURCEID->STR[1] = 'O';
  SOURCEID->STR[2] = 'U';
  SOURCEID->STR[3] = 'R';
  SOURCEID->STR[4] = 'C';
  SOURCEID->STR[5] = 'E';
  SOURCEID->LEN = 6;


  ORGID->STR[0] = 'O';
  ORGID->STR[1] = 'R';
  ORGID->STR[2] = 'G';
  ORGID->STR[3] = 'A';
  ORGID->STR[4] = 'N';
  ORGID->STR[5] = 'I';
  ORGID->STR[6] = 'S';
  ORGID->STR[7] = 'M';
  ORGID->LEN = 8;
  /* These keywords can be eliminated after PIR makes the switch over
     to the #accession protocol */
  OLDACID.STR[0] = 'A';
  OLDACID.STR[1] = 'C';
  OLDACID.STR[2] = 'C';
  OLDACID.STR[3] = 'E';
  OLDACID.STR[4] = 'S';
  OLDACID.STR[5] = 'S';
  OLDACID.STR[6] = 'I';
  OLDACID.STR[7] = 'O';
  OLDACID.STR[8] = 'N';
  OLDACID.STR[9] = 'S';
  OLDACID.LEN = 10;
  REFID.STR[0] = 'R';
  REFID.STR[1] = 'E';
  REFID.STR[2] = 'F';
  REFID.STR[3] = 'E';
  REFID.STR[4] = 'R';
  REFID.STR[5] = 'E';
  REFID.STR[6] = 'N';
  REFID.STR[7] = 'C';
  REFID.STR[8] = 'E';
  REFID.LEN = 9;
}  /* INITKEYS */


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
/* I/O procedures                                              */
/***************************************************************/
Static Void OPENFILES(DBFILE, ANOFILE, SEQFILE, INDFILE, ANOARGNUM)
_TEXT *DBFILE, *ANOFILE, *SEQFILE, *INDFILE;
long *ANOARGNUM;
{

  /* BEGIN MODULE FILEARGS */
  /* This procedure overcomes one of the stupidest aspects of UNIX Pascal,
     namely the fact that filenames in the program statement are supposed to
     be actual UNIX filenames!  To overcome this, the 2-argument version of
     reset and rewrite must be used with string variables.  This module
     need only contain the reset and rewrite statements in any normal
     implementation of Pascal. */
  /* END MODULE FILEARGS         VERSION= 'SUNMODS     Version  8/ 9/94'; */

  FILEARGS(DBFILE, 'I', &ARGNUM);

  /* We need to rewrite ANOFILE later, so save the ARGNUM as ANOARGNUM */
  *ANOARGNUM = ARGNUM;
  FILEARGS(ANOFILE, 'O', &ARGNUM);

  if (((1L << ((long)DATABASE)) & ((1L << ((long)GB)) | (1L << ((long)PIR)) |
	 (1L << ((long)EMBL)) | (1L << ((long)VECTOR)))) != 0)
    FILEARGS(SEQFILE, 'O', &ARGNUM);
  FILEARGS(INDFILE, 'O', &ARGNUM);
}  /* OPENFILES */


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

/* BEGIN MODULE GETWORD */
/*  Read a word from LINE L */
Static Void GETWORD(CLINE, W, P)
LINE CLINE;
WORD *W;
long *P;
{
  long I;

  W->LEN = 0;
  while (CLINE.STR[*P - 1] == ' ')
    (*P)++;
  while (CLINE.STR[*P - 1] != ' ') {
    if (W->LEN < MAXWORD) {
      W->LEN++;
      W->STR[W->LEN - 1] = CLINE.STR[*P - 1];
    }
    (*P)++;
  }
  for (I = W->LEN; I < MAXWORD; I++)
    W->STR[I] = ' ';
}  /* GETWORD */


/* END MODULE GETWORD */

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

/* BEGIN MODULE SAMEWORD */
/* Compare two WORDS for equality */
Static boolean SAMEWORD(W1, W2)
WORD *W1, *W2;
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


#define MINBLANKS       3


/* If the number of leading blanks is greater than 3, represent*/
/* the leading blanks by <CRUNCHFLAG><CRUNCHOFFSET+NUMBLANKS>  */
/* where CRUNCHFLAG is the character used to signal that the   */
/* line contains compressed blanks, and CRUNCHOFFSET+NUMBLANKS */
/* is a character whose ASCII value is equal to number of      */
/* leading blanks that the line should have plus CRUNCHOFFSET. */
/* CRUNCHOFFSET is added to prevent the use of the control     */
/* characters with low ASCII values. Since none of the data-   */
/* bases commonly used with XYLEM should have more then 80     */
/* leading blanks, ASCII 127 should never be exceeded.         */
Local Void WRITECRUNCH(F, CLINE)
_TEXT *F;
LINE *CLINE;
{
  long NUMBLANKS = 3, I = 4;
  long FORLIM;

  if (CLINE->LEN < MINBLANKS)
    return;
  if (CLINE->STR[0] != ' ' || CLINE->STR[1] != ' ' || CLINE->STR[2] != ' ') {
    WRITELINE(F, *CLINE, CLINE->LEN);
    return;
  }
  putc(CRUNCHFLAG, F->f);
  /* Set NUMBLANKS = the number of leading blanks. */
  /* Write the corresponding ASCII character. */
  while (CLINE->STR[I-1] == ' ' && I <= CLINE->LEN) {
    NUMBLANKS++;
    I++;
  }
  putc((Char)(CRUNCHOFFSET + NUMBLANKS), F->f);

  FORLIM = CLINE->LEN;
  /* Write the non-blank portion of the line. */
  for (I = NUMBLANKS; I < FORLIM; I++)
    putc(CLINE->STR[I], F->f);
}  /* WRITECRUNCH */

#undef MINBLANKS


/* END MODULE SAMEWORD         VERSION= 'SUNMODS     Version  8/ 9/94'; */


/***************************************************************/
/* Read in each line and write to output file, until KEYWORD   */
/* is found starting in the column specified by POSITION       */
/* When done, POSITION returns the column containing the next  */
/* unread character in CLINE, after KEYWORD is read.           */
/***************************************************************/
Static Void ADVANCE(CLINE, KEYWORD, POSITION, LINENUM)
LINE *CLINE;
WORD KEYWORD;
long *POSITION, *LINENUM;
{
  Char FLAG;
  long FLAGPOSN;
  WORD ID;

  FLAGPOSN = *POSITION;
  FLAG = KEYWORD.STR[0];
  do {
    do {
      READLINE(&DBFILE, CLINE);
      if (COMPRESSION)
	WRITECRUNCH(&ANOFILE, CLINE);
      else
	WRITELINE(&ANOFILE, *CLINE, CLINE->LEN);
      putc('\n', ANOFILE.f);
      (*LINENUM)++;
    } while (!(BUFEOF(DBFILE.f) || CLINE->STR[FLAGPOSN-1] == FLAG));
    *POSITION = FLAGPOSN;
    GETWORD(*CLINE, &ID, POSITION);
  } while (!(BUFEOF(DBFILE.f) | SAMEWORD(&KEYWORD, &ID)));   /* ADVANCE */
}


/******************************************************************/
/* Transfer lines until the DEFINITION (or equivalent) is reached.*/
/******************************************************************/
Static Void READDEF(CLINE, CANO, CSEQ)
LINE *CLINE;
long *CANO, *CSEQ;
{
  long I = 0, POSITION = 1;

  /* Get DEFINITION line */
  ADVANCE(CLINE, DEFID, &POSITION, CANO);

  /* Read in partial title from DEFINITION line */
  while (CLINE->STR[POSITION-1] == ' ' && POSITION < CLINE->LEN)
    POSITION++;
  while (I < 50 && POSITION <= CLINE->LEN) {
    I++;
    TITLE.STR[I-1] = CLINE->STR[POSITION-1];
    POSITION++;
  }
  TITLE.LEN = I;

  /* Title line:
     >name - parital title from DEFINITION line */
  putc('>', SEQFILE.f);
  WRITEWORD(&SEQFILE, NAME, NAME.LEN);
  fprintf(SEQFILE.f, " - ");
  WRITELINE(&SEQFILE, TITLE, TITLE.LEN);
  putc('\n', SEQFILE.f);
  (*CSEQ)++;
  FIRSTSEQ = *CSEQ;
}  /* READDEF */


/* Local variables for READAC: */
struct LOC_READAC {
  long *CANO;
} ;

/* Search for ACCESSION line */
Local Void FINDAC(CLINE, POSITION, LINENUM, LINK)
LINE *CLINE;
long *POSITION, *LINENUM;
struct LOC_READAC *LINK;
{
  boolean GOTIT = false;
  long FLAGPOSN = 1;
  WORD ID;

  /* Keep searching until either an ACCESSIONS or REFERENCE line is
     found. */
  do {
    do {
      READLINE(&DBFILE, CLINE);
      WRITELINE(&ANOFILE, *CLINE, CLINE->LEN);
      putc('\n', ANOFILE.f);
      (*LINENUM)++;
    } while (!(BUFEOF(DBFILE.f) || CLINE->STR[FLAGPOSN-1] == 'R' ||
	       CLINE->STR[FLAGPOSN-1] == 'A'));
    *POSITION = FLAGPOSN;
    GETWORD(*CLINE, &ID, POSITION);
    if (SAMEWORD(&OLDACID, &ID))
      GOTIT = true;
    else if (SAMEWORD(&REFID, &ID))
      GOTIT = true;
  } while (!(BUFEOF(DBFILE.f) || GOTIT));

  if (SAMEWORD(&REFID, &ID)) {
    /* We're on the REFERENCE line. No ACCESSIONS line was found. Move to
       the first #accession line. */
    *POSITION = 4;
    ADVANCE(CLINE, ACID, POSITION, LINK->CANO);
  }
}  /* FINDAC */


/***************************************************************/
/* Transfer lines until the ACCESSION line is found. Read      */
/* ACCESSION from CLINE.                                       */
/***************************************************************/
Static Void READAC(CLINE, ACCESSION, CANO_)
LINE *CLINE;
WORD *ACCESSION;
long *CANO_;
{
  struct LOC_READAC V;
  long POSITION;


  V.CANO = CANO_;
  /* Get ACCESSION line */
  if (DATABASE == PIR) {
    /* Since the primary accession number could be on an ACCESSIONS
       line (old format) or an #accession line (new format) we
       have be able to find either. In either case, FINDAC will
       set POSITION to the first character of the accession num. */
    FINDAC(CLINE, &POSITION, V.CANO, &V);
    GETWORD(*CLINE, ACCESSION, &POSITION);
  }  /* PIR */

  else {
    POSITION = 1;
    ADVANCE(CLINE, ACID, &POSITION, V.CANO);
    GETWORD(*CLINE, ACCESSION, &POSITION);
  }
  /* PIR - If there were secondary accession numbers, the
     primary accession number will be terminated by a semicolon (;)
     which must be removed */
  if (DATABASE == PIR) {
    if (ACCESSION->STR[ACCESSION->LEN - 1] == ';')
      ACCESSION->LEN--;
  }
}  /* READAC */


/******************************************************************/
/* Transfer lines until the SOURCE line is reached. Source may span*/
/* more than 1 line.                                               */
/* ORGANISM line may take the form                                 */
/* [genome] genus species subspecies                               */
/* so we just read in the rest of the line.                        */
/******************************************************************/
Static Void READSOURCE(CLINE, CANO, ORGANISM)
LINE *CLINE;
long *CANO;
LINE *ORGANISM;
{
  long I = 3, J = 1;

  /* Get SOURCE line */
  /*         I:=1;
           ADVANCE(CLINE,SOURCEID,I,CANO);
  */
  /* Read in ORGANISM line and write to .ano file. */
  ADVANCE(CLINE, ORGID, &I, CANO);

  /* Copy remaining part of line to ORGANISM */
  /* with CLINE */

  while (CLINE->STR[I-1] == ' ' && I <= CLINE->LEN)
    I++;
  while (I <= CLINE->LEN) {
    ORGANISM->STR[J-1] = CLINE->STR[I-1];
    I++;
    J++;
  }
  ORGANISM->LEN = J - 1;
}  /* READSOURCE */


/***************************************************************/
/* Transfer lines until sequence is found. Write sequence to   */
/* SEQFILE, stripping off blanks and numbers.                  */
/***************************************************************/
Static Void READSEQ(CLINE, CSEQ)
LINE *CLINE;
long *CSEQ;
{
  long WIDTH = 0, POSITION = 1;

  /* Advance to ORIGIN line. Sequence begins on next line */
  ADVANCE(CLINE, SEQID, &POSITION, &CANO);
  if (DATABASE == PIR)
    fprintf(ANOFILE.f, "///\n");
  else
    fprintf(ANOFILE.f, "//\n");
  CANO++;

  /* Read in sequence and write to output. */
  /* Omit extra blanks and numbers. */
  while (P_peek(DBFILE.f) != '/') {
    while (!P_eoln(DBFILE.f)) {
      CH = getc(DBFILE.f);
      if (CH == '\n')
	CH = ' ';
      if (CH == ' ' || isdigit(CH)) {
	continue;
      }  /* CH in [' ','0'..'9'] */
      putc(CH, SEQFILE.f);
      WIDTH++;
      if (WIDTH == SEQLINELEN) {
	putc('\n', SEQFILE.f);
	(*CSEQ)++;
	WIDTH = 0;
      }
    }  /* not eoln */
    fscanf(DBFILE.f, "%*[^\n]");
    getc(DBFILE.f);
  }  /* not (DBFILE^='/') */
  if (WIDTH > 0) {
    putc('\n', SEQFILE.f);
    (*CSEQ)++;
  }
  if (!BUFEOF(DBFILE.f)) {
    fscanf(DBFILE.f, "%*[^\n]");
    getc(DBFILE.f);
  }

  /* PIR files terminate with '\\\'. However, splitdb does not
     demand that this terminator line be there. It is not re-
     generated by getloc. */
  if (P_peek(DBFILE.f) == '\\') {
    fscanf(DBFILE.f, "%*[^\n]");
    getc(DBFILE.f);
  }

  /* Read past trailing blank lines, if any */
  while ((!BUFEOF(DBFILE.f)) & (P_peek(DBFILE.f) == ' ')) {
    fscanf(DBFILE.f, "%*[^\n]");
    getc(DBFILE.f);
  }
}  /* READSEQ */


/* ----------------------------------------------------------  */
/* ----------------- MAIN  PROCEDURE  -----------------------  */
main(argc, argv)
int argc;
Char *argv[];
{
  PASCAL_MAIN(argc, argv);
  INDFILE.f = NULL;
  strcpy(INDFILE.name, "INDFILE");
  SEQFILE.f = NULL;
  strcpy(SEQFILE.name, "SEQFILE");
  ANOFILE.f = NULL;
  strcpy(ANOFILE.name, "ANOFILE");
  DBFILE.f = NULL;
  strcpy(DBFILE.name, "DBFILE");
  /* Determine which database is being processed. */
  READOPTIONS(&ARGNUM, &DATABASE, &COMPRESSION, &TAXONOMY);

  /* Initialize keywords according to database type. */
  INITKEYS(&NAMEID, &DEFID, &ACID, &SEQID, &SOURCEID, &ORGID);

  /* Open files. */
  OPENFILES(&DBFILE, &ANOFILE, &SEQFILE, &INDFILE, &ANOARGNUM);
  P_sun_argv(ANOFILENAME, 132, (int)ANOARGNUM);
      /* save ANOFILENAME for rewrite */

  CANO = 0;
  CSEQ = 0;

  /* Advance to first entry and throw away header information
     by rewriting ANOFILE.*/
  POSITION = 1;
  CRUNCHFLAG = (Char)CRUNCHOFFSET;
  ADVANCE(&CLINE, NAMEID, &POSITION, &CANO);

  strcpy(ANOFILE.name, P_trimname(ANOFILENAME, 132));
  if (ANOFILE.f != NULL)
    ANOFILE.f = freopen(ANOFILE.name, "w", ANOFILE.f);
  else
    ANOFILE.f = fopen(ANOFILE.name, "w");
  if (ANOFILE.f == NULL)
    _EscIO2(FileNotFound, ANOFILE.name);
  SETUPBUF(ANOFILE.f, Char);

  WRITELINE(&ANOFILE, CLINE, CLINE.LEN);
  putc('\n', ANOFILE.f);
  CANO = 1;
  GETWORD(CLINE, &NAME, &POSITION);
  WRITEWORD(&INDFILE, NAME, 16L);
  FIRSTANO = CANO;
  /* Note on integer output to INDFILE: */
  /* ATT Pascal writes integers using a field width of 10 unless a field
     width is set.  If the field width is too small to print the integer,
     the minimum required field width is used.  Thus, I have to set a
     field width to insure that the minimum field width is used in
     all cases, if I want to minimize the size of the output file.*/

  while (!BUFEOF(DBFILE.f)) {
    switch (DATABASE) {

    case GB:
    case PIR:
    case VECTOR:
      READDEF(&CLINE, &CANO, &CSEQ);
      READAC(&CLINE, &ACCESSION, &CANO);
      if (TAXONOMY && DATABASE == GB)
	READSOURCE(&CLINE, &CANO, &ORGANISM);
      READSEQ(&CLINE, &CSEQ);
      putc(' ', INDFILE.f);
      WRITEWORD(&INDFILE, ACCESSION, 8L);
      fprintf(INDFILE.f, " %8ld %8ld", FIRSTANO, FIRSTSEQ);
      if (TAXONOMY && DATABASE == GB) {
	putc(' ', INDFILE.f);
	WRITELINE(&INDFILE, ORGANISM, ORGANISM.LEN);
      }
      putc('\n', INDFILE.f);
      break;
      /* GB,PIR,VECTOR */

    case EMBL:
      READAC(&CLINE, &ACCESSION, &CANO);
      READDEF(&CLINE, &CANO, &CSEQ);
      READSEQ(&CLINE, &CSEQ);
      putc(' ', INDFILE.f);
      WRITEWORD(&INDFILE, ACCESSION, 8L);
      fprintf(INDFILE.f, " %8ld %8ld\n", FIRSTANO, FIRSTSEQ);
      break;
      /* EMBL */

    }

    /* Advance to next LOCUS line. Read in NAME. */
    /* Write NAME to INDFILE and set FIRSTANO to current line number*/
    if (BUFEOF(DBFILE.f))
      break;
    POSITION = 1;
    ADVANCE(&CLINE, NAMEID, &POSITION, &CANO);
    GETWORD(CLINE, &NAME, &POSITION);
    WRITEWORD(&INDFILE, NAME, 16L);
    FIRSTANO = CANO;
  }  /* while not eof */


  /*!!!CLOSE(ANOFILE); CLOSE(SEQFILE); CLOSE(INDFILE)*/
  if (DBFILE.f != NULL)
    fclose(DBFILE.f);
  if (ANOFILE.f != NULL)
    fclose(ANOFILE.f);
  if (SEQFILE.f != NULL)
    fclose(SEQFILE.f);
  if (INDFILE.f != NULL)
    fclose(INDFILE.f);
  exit(EXIT_SUCCESS);
}  /* SPLITDB  */



/* End. */
