/* Output from p2c 1.21alpha-07.Dec.93, the Pascal-to-C translator */
/* From input file "multidigest.p" */


/**********************************************************/
/*                                                        */
/* MULTIDIGEST  Version  3/29/06  Standard Pascal         */
/*              Brian Fristensky                          */
/*              Dept. of Plant Science                    */
/*              University of Manitoba                    */
/*              Winnipeg, MB R3T 2N2  CANADA              */
/*                                                        */
/* Copyright (c) 1987,1990 by Brian Fristensky            */
/* !!! in comment indicates feature which may need change */
/******************************************************** */

/* REVISION HISTORY
29 Mar 2006 Changed name of program from DIGEST to MULTIDIGEST,
            to avoid conflict with Solaris 10 digest command.
27 Mar 2006 READSITES revised to read new BACHREST output files
            with parameter information. NOT backward compatible
            with old BACHREST files!!
12 Aug 2001  MAXENZ increased from 500 to 4000 due to growth in REBASE
            MAXSITES increased to 100000
*/


#include <p2c.h>


/*,INFILE,OUTFILE*/
/*!!! Some Pascals may require file parameters in program heading */

#define MAXSITES        100000L   /*Max. number of sites in all digests */

#define MAXENZ          4000   /*Max. number of enzymes read in      */
#define MAXFRAGS        500   /*Max. number of fragments in a digest*/
#define MAXPAT          23   /*Max. length of a restriction site   */
#define MAXWORD         23   /*Max. length of a WORD               */
#define MAXLINE         150

#define VERSION         "MULTIDIGEST       Version 29 Mar 2006"


/* BEGIN MODULE TYPE.WORD */
/*   <word>::= <non-blank char>[<non-blank char>] */

typedef struct WORD {
  long LEN;
  Char STR[MAXWORD];
} WORD;

/* END MODULE TYPE.WORD         VERSION= 'SUNMODS     Version  6/26/01'; */

/* BEGIN MODULE TYPE.LINE */
typedef Char CHARARRAY[MAXLINE];

typedef struct LINE {
  CHARARRAY STR;
  long LEN;
} LINE;

/* END MODULE TYPE.LINE         VERSION= 'SUNMODS     Version  6/26/01'; */

/* RESTRICTION ENZYME */

typedef struct ENZYME {
  WORD NAME, RSEQ;
  long FIRSTSITE;   /* Rest. sites are stored in  */
  /* SITES[FIRSTSITE..FIRSTSITE+NUMSITES]*/
  long NUMSITES;
} ENZYME;

/* LIST OF RESTRICTION SITES FOUND */

typedef struct FRAGMENT {
  long START, FINISH, SIZE, STARTENZ, FINISHENZ;
  struct FRAGMENT *PREV, *NEXT;
} FRAGMENT;

typedef struct FRAGSFOUND {
  long LNUM;
  FRAGMENT *HEAD, *TAIL;
} FRAGSFOUND;


Static _TEXT INFILE, OUTFILE;
Static LINE IFN, OFN, HLINE;
Static ENZYME LIST[MAXENZ + 1];   /*Enzyme sites found */
Static long SITES[MAXSITES];
Static long NUMENZ;   /* Number of sites found*/
Static WORD NAME;   /* Name of sequence searched */
Static FRAGSFOUND FOUND, PFOUND;   /* Frags. in current digest */
Static FRAGMENT *FREEFRAG;   /*Points to freelist of sites*/
Static long FREESITE;
Static boolean CIRCULAR;   /*=true if seq. is circular */
Static long I, CHOICE, SEQLEN;   /* length of sequence */


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


/* END MODULE INPLINE         VERSION= 'SUNMODS     Version  6/26/01'; */

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


/* END MODULE WRITELINE         VERSION= 'SUNMODS     Version  6/26/01'; */

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


/* END MODULE GETFILE         VERSION= 'SUNMODS     Version  6/26/01'; */

/**************************************************/
/*  WORD  I/O procedures                          */
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


/* END MODULE INPWORD         VERSION= 'SUNMODS     Version  6/26/01'; */

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


/* END MODULE READWORD         VERSION= 'SUNMODS     Version  6/26/01'; */

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


/* END MODULE WRITEWORD         VERSION= 'SUNMODS     Version  6/26/01'; */

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


/* END MODULE GETINTEGER         VERSION= 'SUNMODS     Version  6/26/01'; */

/**************************************************************/
/*  Read a list of restriction sites from INFILE.             */
/**************************************************************/
Static Void READSITES()
{
  WORD DUMMY, TOPOLOGY;
  long I, NUMREAD;
  Char CH;
  boolean FRAGSFOUND_;
  _TEXT TEMP;
  long FORLIM;
  ENZYME *WITH;

  printf("Reading restriction site data...\n");
  /* Read title box */
  fscanf(INFILE.f, "%*[^\n]");
  getc(INFILE.f);
  fscanf(INFILE.f, "%*[^\n]");
  getc(INFILE.f);
  fscanf(INFILE.f, "%*[^\n]");
  getc(INFILE.f);
  READWORD(&INFILE, &NAME);
  READWORD(&INFILE, &DUMMY);
  TEMP.f = stdout;
  *TEMP.name = '\0';
  /*!!!*/
  WRITEWORD(&TEMP, NAME, NAME.LEN);
  putchar(' ');
  TEMP.f = stdout;
  *TEMP.name = '\0';
  WRITEWORD(&TEMP, DUMMY, DUMMY.LEN);
  putchar(' ');

  READWORD(&INFILE, &TOPOLOGY);
  READWORD(&INFILE, &DUMMY);
  TEMP.f = stdout;
  *TEMP.name = '\0';
  /*!!!*/
  WRITEWORD(&TEMP, TOPOLOGY, TOPOLOGY.LEN);
  putchar(' ');
  TEMP.f = stdout;
  *TEMP.name = '\0';
  WRITEWORD(&TEMP, DUMMY, DUMMY.LEN);
  putchar(' ');
  if (TOPOLOGY.STR[0] == 'C')
    CIRCULAR = true;
  else
    CIRCULAR = false;
  fscanf(INFILE.f, "%ld%*[^\n]", &SEQLEN);
  getc(INFILE.f);
  /*!!!*/
  printf("%12ld\n", SEQLEN);

  /* Read parameter box */
  for (I = 1; I <= 13; I++) {
    fscanf(INFILE.f, "%*[^\n]");
    getc(INFILE.f);
  }

  NUMENZ = 0;
  FREESITE = 1;   /*first element in site array */

  while (!BUFEOF(INFILE.f)) {
    NUMENZ++;
    WITH = &LIST[NUMENZ];
    READWORD(&INFILE, &WITH->NAME);
    READWORD(&INFILE, &WITH->RSEQ);
    TEMP.f = stdout;
    *TEMP.name = '\0';
    /*!!!*/
    WRITEWORD(&TEMP, WITH->NAME, WITH->NAME.LEN);
    putchar(' ');
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITEWORD(&TEMP, WITH->RSEQ, WITH->RSEQ.LEN);
    printf(" \n");

    /* Old code, prior to reading REBASE style sites */
    /* Read cutting position.  Asymmetric enzymes will have a */
    /* second cutting position in parentheses.                */
    /*        read(INFILE,CUT); read(INFILE,CH);
              while INFILE^=' ' do read(INFILE,CH);
              if INFILE^='(' then begin
                while INFILE^<> ')' do read(INFILE,CH);
                read(INFILE,CH) end;                      */

    /* Read in restriction sites for enzyme */
    fscanf(INFILE.f, "%ld%*[^\n]", &WITH->NUMSITES);
    getc(INFILE.f);

    /* Check to see if fragments are shown. BACHREST will
       print a message in parentheses,
       rather than fragments, if there
       are too many fragments. If that case, we skip
       this enzyme. */
    while (P_peek(INFILE.f) == ' ') {
      CH = getc(INFILE.f);
      if (CH == '\n')
	CH = ' ';
    }
    if (P_peek(INFILE.f) == '(') {
      FRAGSFOUND_ = false;
      fscanf(INFILE.f, "%*[^\n]");
      getc(INFILE.f);
    } else
      FRAGSFOUND_ = true;

    if (WITH->NUMSITES > 0 && FRAGSFOUND_) {
      WITH->FIRSTSITE = FREESITE;
      I = WITH->FIRSTSITE;
      FORLIM = WITH->NUMSITES;
      for (NUMREAD = 1; NUMREAD <= FORLIM; NUMREAD++) {
	fscanf(INFILE.f, "%ld%*[^\n]", &SITES[I-1]);
	getc(INFILE.f);
	I++;
      }
      FREESITE += WITH->NUMSITES;
    } else
      NUMENZ--;
    /* Enzymes with no sites are ignored */

    if (!CIRCULAR) {
      fscanf(INFILE.f, "%*[^\n]");
      getc(INFILE.f);
    }
    while ((P_peek(INFILE.f) == ' ') & (!BUFEOF(INFILE.f))) {
      CH = getc(INFILE.f);
      if (CH == '\n')
	CH = ' ';
    }
  }
  /*!!!*/
  if (INFILE.f != NULL)
    fclose(INFILE.f);
  INFILE.f = NULL;
}  /* READSITES */


/* BEGIN MODULE LINKED */
/*********************************************************/
/*  Linked-list operations for restriction fragment list.*/
/*********************************************************/

/*Get a new fragment from freelist.*/
Static Void GETFRAG(NEWFRAG)
FRAGMENT **NEWFRAG;
{
  if (FREEFRAG == NULL)
    *NEWFRAG = (FRAGMENT *)Malloc(sizeof(FRAGMENT));
  else {
    *NEWFRAG = FREEFRAG;
    FREEFRAG = FREEFRAG->NEXT;
  }
}


/*Add a fragment after DEST*/
Static Void ADDFRAG(AFRAG, DEST)
FRAGMENT **AFRAG, **DEST;
{
  FRAGMENT *TEMP;

  TEMP = (*DEST)->NEXT;
  (*DEST)->NEXT = *AFRAG;
  (*AFRAG)->PREV = *DEST;
  (*AFRAG)->NEXT = TEMP;
  TEMP->PREV = *AFRAG;
}


/*Return a list to the top of freelist*/
Static Void RIDOF(HEAD, TAIL)
FRAGMENT **HEAD, **TAIL;
{
  FRAGMENT *TEMPHEAD, *TEMPTAIL;

  if ((*HEAD)->NEXT == *TAIL)
    return;
  TEMPHEAD = (*HEAD)->NEXT;
  TEMPTAIL = (*TAIL)->PREV;
  (*HEAD)->NEXT = *TAIL;
  (*TAIL)->PREV = *HEAD;
  TEMPHEAD->PREV = NULL;
  TEMPTAIL->NEXT = FREEFRAG;
  FREEFRAG = TEMPHEAD;
}


/* END MODULE LINKED         VERSION= 'SUNMODS     Version  6/26/01'; */

/************************************************/
/*  Re-initialize linked lists.                 */
/************************************************/
Static Void INITLIST(F)
FRAGSFOUND *F;
{
  F->HEAD->NEXT = F->TAIL;
  F->TAIL->PREV = F->HEAD;
  F->HEAD->START = 0;
  F->HEAD->SIZE = 0;
  F->TAIL->START = SEQLEN + 1;
  F->TAIL->SIZE = SEQLEN + 1;
}  /* INITLIST */


/************************************************/
/*  Perform one-time initializations.           */
/************************************************/
Static Void INITIALIZE()
{
  LIST[0].NAME.LEN = 0;
  /* Initialize the linked list FOUND */
  FREEFRAG = NULL;
  FOUND.HEAD = (FRAGMENT *)Malloc(sizeof(FRAGMENT));
  FOUND.TAIL = (FRAGMENT *)Malloc(sizeof(FRAGMENT));
  PFOUND.HEAD = (FRAGMENT *)Malloc(sizeof(FRAGMENT));
  PFOUND.TAIL = (FRAGMENT *)Malloc(sizeof(FRAGMENT));
  INITLIST(&FOUND);
  INITLIST(&PFOUND);
}


/* Local variables for WHICH: */
struct LOC_WHICH {
  _TEXT *OUTFILE;
  long NUMDISPLAYED;
  Char INCSYM[2];
  boolean INCLUDED[MAXENZ];
} ;

/* Local variables for DISPLAY: */
struct LOC_DISPLAY {
  struct LOC_WHICH *LINK;
} ;

/* Local variables for MERGE: */
struct LOC_MERGE {
  struct LOC_DISPLAY *LINK;
  long ANS;
  FRAGMENT *RIGHTMOST;
} ;

/*  Add a node to the linked list FOUND, telling the cutting site.*/
Local Void ADDSITE(POSN, LINK)
long POSN;
struct LOC_MERGE *LINK;
{
  FRAGMENT *NEWFRAG;

  GETFRAG(&NEWFRAG);
  NEWFRAG->START = POSN;
  NEWFRAG->STARTENZ = LINK->ANS;
  /* Find the rightmost site <= POSN */
  while (LINK->RIGHTMOST->NEXT->START <= POSN)
    LINK->RIGHTMOST = LINK->RIGHTMOST->NEXT;
  /* if RIGHTMOST site <> POSN, add a site */
  if (POSN <= LINK->RIGHTMOST->START)
    return;
  ADDFRAG(&NEWFRAG, &LINK->RIGHTMOST);
  FOUND.LNUM++;
  LINK->RIGHTMOST = NEWFRAG;
}  /* ADDSITE */

/*  Merge an enzyme with the list of sites found */
Local Void MERGE(ANS_, LINK)
long ANS_;
struct LOC_DISPLAY *LINK;
{
  struct LOC_MERGE V;
  long NUMADDED = 0;
  long S;
  ENZYME *WITH;

  V.LINK = LINK;
  V.ANS = ANS_;
  WITH = &LIST[V.ANS];
  WRITEWORD(LINK->LINK->OUTFILE, WITH->NAME, 10L);
  WRITEWORD(LINK->LINK->OUTFILE, WITH->RSEQ, (long)MAXPAT);
  fprintf(LINK->LINK->OUTFILE->f, "%4ld\n", WITH->NUMSITES);
  S = WITH->FIRSTSITE;
  V.RIGHTMOST = FOUND.HEAD;
  while (NUMADDED < WITH->NUMSITES) {
    ADDSITE(SITES[S-1], &V);
    NUMADDED++;
    S++;
  }
  LINK->LINK->INCLUDED[V.ANS-1] = true;
}  /* MERGE */

/* Display enzymes on screen and add enzymes to FOUND */
Local Void DISPLAY(FIRST, LAST, LINK)
long FIRST, LAST;
struct LOC_WHICH *LINK;
{
  struct LOC_DISPLAY V;
  long ANS, I, WIDTH;
  _TEXT TEMP;


  V.LINK = LINK;
  do {
    fprintf(stdout, "\f");
    WIDTH = 0;
    for (I = FIRST; I <= LAST; I++) {
      printf("%3ld)%c", I, LINK->INCSYM[LINK->INCLUDED[I-1]]);
      TEMP.f = stdout;
      *TEMP.name = '\0';
      WRITEWORD(&TEMP, LIST[I].NAME, 10L);
      WIDTH++;
      if (WIDTH == 5) {
	putchar('\n');
	WIDTH = 0;
      }
    }
    printf("\nType number of an enzyme to include in this digest\n");
    printf("or 0 to continue:\n");
    GETINTEGER(&ANS, 0L, LAST);
    scanf("%*[^\n]");
    getchar();
    if (ANS >= FIRST && ANS <= LAST)
      MERGE(ANS, &V);
    else if (ANS != 0)
      printf("%12ld OUT OF RANGE\n", ANS);
  } while (ANS != 0);
  LINK->NUMDISPLAYED = LAST;
}  /* DISPLAY */


/***************************************************************/
/*  Ask user for enzymes he wishes included in the digest      */
/***************************************************************/
Static Void WHICH(OUTFILE_)
_TEXT *OUTFILE_;
{
  struct LOC_WHICH V;
  long FIRST, LAST, I, FORLIM;
  _TEXT TEMP;


  V.OUTFILE = OUTFILE_;
  /* Initialize inclusion array.  INCLUDED[n] is true when LIST[n] has
     been included in the current digest.  INCSYM holds a '+' to indicate
     inclusion and a ' ' to indicate non-inclusion. */
  V.INCSYM[true] = '+';
  V.INCSYM[0] = ' ';
  FORLIM = NUMENZ;
  for (I = 0; I < FORLIM; I++)
    V.INCLUDED[I] = false;

  /* Print instructions */
  printf("The names of enzymes with known sites in ");
  TEMP.f = stdout;
  *TEMP.name = '\0';
  WRITEWORD(&TEMP, NAME, NAME.LEN + 1);
  printf("\nwill be displayed a screenful at a time.\n");
  printf("You will be asked to specify enzymes one at a time,\n");
  printf("to include in this digest.\n");
  printf("There are %12ld enzymes listed\n", NUMENZ);
  printf("Press RETURN to begin.\n");
  scanf("%*[^\n]");
  getchar();
  V.NUMDISPLAYED = 0;
  FOUND.LNUM = 0;
  while (V.NUMDISPLAYED < NUMENZ) {
    FIRST = V.NUMDISPLAYED + 1;
    if (NUMENZ - V.NUMDISPLAYED > 100)
      LAST = V.NUMDISPLAYED + 100;
    else
      LAST = NUMENZ;
    DISPLAY(FIRST, LAST, &V);
  }
}  /* WHICH */


/* Local variables for REPORT: */
struct LOC_REPORT {
  _TEXT *OUTFILE;
  FRAGMENT *ORDER[MAXFRAGS];
} ;

/*Calculate sizes and ends of fragments*/
Local Void CALCULATE(LINK)
struct LOC_REPORT *LINK;
{
  FRAGMENT *CURRENTFRAG;
  long SITE = 1;

  /* if LINEAR add a fragment to head of list*/
  if (!CIRCULAR && FOUND.HEAD->NEXT->START != 1) {
    GETFRAG(&CURRENTFRAG);
    ADDFRAG(&CURRENTFRAG, &FOUND.HEAD);
    CURRENTFRAG->START = 1;
    CURRENTFRAG->STARTENZ = 0;
    FOUND.LNUM++;
  } else
    CURRENTFRAG = FOUND.HEAD->NEXT;

  if (FOUND.LNUM <= 0)
    return;
  /*Calculate ends and size of each fragment and assign it */
  /* a place in the ORDER array, to be sorted later.    */
  while (CURRENTFRAG->NEXT != FOUND.TAIL) {
    CURRENTFRAG->FINISH = CURRENTFRAG->NEXT->START - 1;
    CURRENTFRAG->FINISHENZ = CURRENTFRAG->NEXT->STARTENZ;
    CURRENTFRAG->SIZE = CURRENTFRAG->FINISH - CURRENTFRAG->START + 1;
    LINK->ORDER[SITE-1] = CURRENTFRAG;
    SITE++;
    CURRENTFRAG = CURRENTFRAG->NEXT;
  }
  /*Last fragment in list is a special case*/
  if (CIRCULAR) {
    CURRENTFRAG->FINISH = FOUND.HEAD->NEXT->START - 1;
    CURRENTFRAG->FINISHENZ = FOUND.HEAD->NEXT->STARTENZ;
    CURRENTFRAG->SIZE = SEQLEN - CURRENTFRAG->START + CURRENTFRAG->FINISH + 1;
    if (CURRENTFRAG->FINISH == 0)
      CURRENTFRAG->FINISH = SEQLEN;
  } else {
    CURRENTFRAG->FINISH = SEQLEN;
    CURRENTFRAG->FINISHENZ = 0;
    CURRENTFRAG->SIZE = CURRENTFRAG->FINISH - CURRENTFRAG->START + 1;
  }
  LINK->ORDER[SITE-1] = CURRENTFRAG;
}

Local Void SWAP(FIRST, SECOND)
FRAGMENT **FIRST, **SECOND;
{
  FRAGMENT *TEMP;

  TEMP = *FIRST;
  *FIRST = *SECOND;
  *SECOND = TEMP;
}

/* BEGIN MODULE SORT */
/*  Invariant:  The array elements > TOP are sorted.*/
/*    TOP >= unsorted elements.                     */
Local Void BUBBLESORT(TOP, BOTTOM, LINK)
long TOP, BOTTOM;
struct LOC_REPORT *LINK;
{
  long SMALLEST, NEXT;

  while (TOP >= BOTTOM) {
    /*bubble smallest unsorted number to the top of sorted list*/
    SMALLEST = BOTTOM;
    NEXT = BOTTOM + 1;
    while (NEXT <= TOP) {
      if (LINK->ORDER[SMALLEST-1]->SIZE < LINK->ORDER[NEXT-1]->SIZE)
	SWAP(&LINK->ORDER[SMALLEST-1], &LINK->ORDER[NEXT-1]);
      SMALLEST = NEXT;
      NEXT = SMALLEST + 1;
    }
    TOP--;
  }
}  /* BUBBLESORT */

/* Local variables for PARTIAL: */
struct LOC_PARTIAL {
  struct LOC_REPORT *LINK;
  FRAGMENT *STARTFRAG, *ENDFRAG;
  long LENGTH;
} ;

/* Create a partial fragment */
Local Void MAKEFRAG(F, LINK)
FRAGMENT **F;
struct LOC_PARTIAL *LINK;
{
  FRAGMENT *WITH;

  GETFRAG(F);
  WITH = *F;
  WITH->STARTENZ = LINK->STARTFRAG->STARTENZ;
  WITH->START = LINK->STARTFRAG->START;
  WITH->SIZE = LINK->LENGTH;
  WITH->FINISH = LINK->ENDFRAG->FINISH;
  WITH->FINISHENZ = LINK->ENDFRAG->FINISHENZ;
}  /* MAKEFRAG */

/* END MODULE SORT         VERSION= 'SUNMODS     Version  6/26/01'; */

/* Calculate fragments in a partial digest, using complete digest
   as a guide. */
Local Void PARTIAL(FOUND, PFOUND, LINK)
FRAGSFOUND *FOUND, *PFOUND;
struct LOC_REPORT *LINK;
{
  struct LOC_PARTIAL V;
  FRAGMENT *RIGHTMOST, *PFRAG;
  long I, J, FORLIM;

  V.LINK = LINK;
  /* if CIRCULAR, circularize the complete digest linked list. */
  /* with FOUND */

  if (CIRCULAR) {
    FOUND->HEAD->NEXT->PREV = FOUND->TAIL->PREV;
    FOUND->TAIL->PREV->NEXT = FOUND->HEAD->NEXT;
  }

  /* Using the complete digest list as a guide, calculate the
     partial fragments and store in PFOUND.  PFOUND is sorted as it
     is built.  Invariant: PFOUND is sorted by SIZE, increasing
     from HEAD to TAIL. */
  V.STARTFRAG = FOUND->HEAD->NEXT;
  PFOUND->LNUM = 0;
  FORLIM = FOUND->LNUM;
  for (I = 1; I <= FORLIM; I++) {  /* for each starting point */
    V.ENDFRAG = V.STARTFRAG;
    V.LENGTH = V.STARTFRAG->SIZE;
    RIGHTMOST = PFOUND->HEAD;
    while (V.LENGTH <= SEQLEN) {  /* for each end point */
      /* Creat a parital fragment*/
      MAKEFRAG(&PFRAG, &V);
      /* Place the fragment into PFOUND list */
      while (V.LENGTH > RIGHTMOST->NEXT->SIZE)
	RIGHTMOST = RIGHTMOST->NEXT;
      ADDFRAG(&PFRAG, &RIGHTMOST);
      PFOUND->LNUM++;
      RIGHTMOST = PFRAG;
      V.ENDFRAG = V.ENDFRAG->NEXT;
      V.LENGTH += V.ENDFRAG->SIZE;
    }
    V.STARTFRAG = V.STARTFRAG->NEXT;
  }

  /* Get rid of the complete digest list and reassign FOUND to the
     parital list. */
  RIDOF(&FOUND->HEAD, &FOUND->TAIL);
  FOUND->HEAD->NEXT = PFOUND->HEAD->NEXT;
  FOUND->TAIL->PREV = PFOUND->TAIL->PREV;
  FOUND->HEAD->NEXT->PREV = FOUND->HEAD;
  FOUND->TAIL->PREV->NEXT = FOUND->TAIL;
  FOUND->LNUM = PFOUND->LNUM;

  /* Assign sorted values in FOUND to the ORDER array, in descending
     order of size. */
  RIGHTMOST = FOUND->TAIL->PREV;
  FORLIM = FOUND->LNUM;
  for (J = 0; J < FORLIM; J++) {
    LINK->ORDER[J] = RIGHTMOST;
    RIGHTMOST = RIGHTMOST->PREV;
  }
  /* Reinitialize PFOUND */
  PFOUND->HEAD->NEXT = PFOUND->TAIL;
  PFOUND->TAIL->PREV = PFOUND->HEAD;
}  /* PARTIAL*/

/*  Print the list of sites found.        */
Local Void PRINTLIST(LINK)
struct LOC_REPORT *LINK;
{
  long SITE, FORLIM;
  FRAGMENT *WITH;

  FORLIM = FOUND.LNUM;
  for (SITE = 0; SITE < FORLIM; SITE++) {
    WITH = LINK->ORDER[SITE];
    fprintf(LINK->OUTFILE->f, "%*ld", (int)(MAXPAT + 22), WITH->SIZE);
    fprintf(LINK->OUTFILE->f, "%8ld", WITH->START);
    WRITEWORD(LINK->OUTFILE, LIST[WITH->STARTENZ].NAME, 10L);
    fprintf(LINK->OUTFILE->f, "%6ld", WITH->FINISH);
    WRITEWORD(LINK->OUTFILE, LIST[WITH->FINISHENZ].NAME, 10L);
    putc('\n', LINK->OUTFILE->f);
  }
  RIDOF(&FOUND.HEAD, &FOUND.TAIL);
  putc('\n', LINK->OUTFILE->f);
}


/*********************************************/
/* Compile a report for output.              */
/*********************************************/
Static Void REPORT(OUTFILE_)
_TEXT *OUTFILE_;
{
  struct LOC_REPORT V;
  Char RESPONSE;

  V.OUTFILE = OUTFILE_;
  CALCULATE(&V);
  do {
    printf("Type C for complete digest, P for parital:\n");
    scanf("%c%*[^\n]", &RESPONSE);
    getchar();
    if (RESPONSE == '\n')
      RESPONSE = ' ';
  } while (RESPONSE != 'p' && RESPONSE != 'P' && RESPONSE != 'c' &&
	   RESPONSE != 'C');
  switch (RESPONSE) {

  case 'C':
  case 'c':
    BUBBLESORT(FOUND.LNUM, 1L, &V);
    break;

  case 'P':
  case 'p':
    PARTIAL(&FOUND, &PFOUND, &V);
    break;
  }
  PRINTLIST(&V);
}  /* REPORT */


/***************************************************************/
/* Prompt user for sites and call SEARCH.                      */
/***************************************************************/
Static Void PROMPT(OUTFILE)
_TEXT *OUTFILE;
{
  Char ANSWER;

  fprintf(OUTFILE->f, "%s\n", VERSION);
  WRITEWORD(OUTFILE, NAME, NAME.LEN);
  fprintf(OUTFILE->f, "  Topology: ");
  switch (CIRCULAR) {

  case true:
    fprintf(OUTFILE->f, " CIRCULAR");
    break;

  case false:
    fprintf(OUTFILE->f, " LINEAR");
    break;
  }
  fprintf(OUTFILE->f, " Length: %12ld bp\n", SEQLEN);
  fprintf(OUTFILE->f, "%*s\n", (int)(MAXPAT + 14), "# of");
  fprintf(OUTFILE->f, "%*s\n",
	  (int)(MAXPAT + 46), "Sites   Frags   Begin            End");
  do {
    WHICH(OUTFILE);
    REPORT(OUTFILE);
    do {
      printf("Type  D to generate a digest, Q to quit\n");
      scanf("%c%*[^\n]", &ANSWER);
      getchar();
      if (ANSWER == '\n')
	ANSWER = ' ';
    } while (ANSWER != 'q' && ANSWER != 'Q' && ANSWER != 'd' && ANSWER != 'D');
  } while (ANSWER != 'q' && ANSWER != 'Q');   /* PROMPT */
}


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
  INFILE.f = NULL;
  *INFILE.name = '\0';
  printf("%50s\n\n", VERSION);
  /* END MODULE STARTUP         VERSION= 'SUNMODS     Version  6/26/01'; */
  INITIALIZE();
  /* Initialize horizontal output line */

  /* Open restriction site file */
  for (I = 1; I <= MAXLINE; I++)
    HLINE.STR[I-1] = '_';
  HLINE.LEN = MAXLINE;
  printf("Enter restriction site filename:\n");
  GETFILE(&INFILE, 'I', &IFN);
  NAME.LEN = 0;
  READSITES();
  if (NAME.LEN == 0) {
    printf("Type name to appear on output:\n");
    INPWORD(&NAME);
    scanf("%*[^\n]");
    getchar();
  }
  INITLIST(&FOUND);
  INITLIST(&PFOUND);

  OFN.LEN = 0;   /* indicates that OUTFILE is not yet open */

  /* MAIN MENU */
  do {
    putchar('\n');
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, HLINE, 80L);
    printf("\nMULTIDIGEST%25s\n", "MAIN MENU");
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, HLINE, 80L);
    printf("\nInput file:        ");
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, IFN, IFN.LEN);
    printf("\nOutput file:       ");
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, OFN, OFN.LEN);
    putchar('\n');
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, HLINE, 80L);
    printf("\n%20c1) Read in a new input file\n", ' ');
    printf("%20c2) Open a new output file\n", ' ');
    printf("%20c3) Generate digests (output to screen)\n", ' ');
    printf("%20c4) Generate digests (output to file)\n", ' ');
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, HLINE, 80L);
    printf("\nType the number of your choice  (0 to quit program)\n");
    GETINTEGER(&CHOICE, 0L, 4L);
    scanf("%*[^\n]");
    getchar();

    switch (CHOICE) {

    case 0:
      /* blank case */
      break;

    case 1:
      printf("Enter restriction site filename:\n");
      GETFILE(&INFILE, 'I', &IFN);
      NAME.LEN = 0;
      READSITES();
      if (NAME.LEN == 0) {
	printf("Type name to appear on output:\n");
	INPWORD(&NAME);
	scanf("%*[^\n]");
	getchar();
      }
      INITLIST(&FOUND);
      INITLIST(&PFOUND);
      break;

    case 2:
      printf("Type output filename:\n");
      GETFILE(&OUTFILE, 'O', &OFN);
      break;

    /*!!!         CLOSE(OUTFILE); */
    case 3:
      TEMP.f = stdout;
      *TEMP.name = '\0';
      PROMPT(&TEMP);
      break;

    case 4:
      if (OFN.LEN > 0)
	PROMPT(&OUTFILE);
      else {
	printf(">>> NO OUTPUT FILE CURRENTLY OPEN\n");
	scanf("%*[^\n]");
	getchar();
      }
      break;
    }
  } while (CHOICE != 0);
  /*!!!*/
  if (OUTFILE.f != NULL)
    fclose(OUTFILE.f);
  OUTFILE.f = NULL;
  if (INFILE.f != NULL)
    fclose(INFILE.f);
  if (OUTFILE.f != NULL)
    fclose(OUTFILE.f);
  exit(EXIT_SUCCESS);
}  /* MULTIDIGEST */




/* End. */
