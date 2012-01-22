/* Output from p2c 1.21alpha-07.Dec.93, the Pascal-to-C translator */
/* From input file "numseq.p" */


/***********************************************************/
/*                                                         */
/*  NUMSEQ    VERSION   6/26/2001  Standard Pascal         */
/*            Brian Fristensky                             */
/*            Dept. of Plant Science                       */
/*            University of Manitoba                       */
/*            Winnipeg, MB R3T 2N2  CANADA                 */
/*                                                         */
/*  Copyright (c) l982,1986,1988,1990  by Brian Fristensky */
/*  !!! in comment indicates feature which may need change */
/***********************************************************/

/* Revision history
06/26/2001 - Changed parameter menu to allow sequence lengths
             up to 11 digits for compatibility with GenBank.
11/05/2000 - In previous versions, if START and FINISH imply
             circularity in a linear molecule, NUMSEQ would
             just print through the ends. For example, if
             START > FINISH on the input strand, NUMSEQ
             would incorrectly print to the end, and then
             resume printing at 1. Now, NUMSEQ would correctly
             "fall off the end" of a linear molecule.

*/


#include <p2c.h>


/*,INFILE,OUTFILE,GCFILE*/
/*!!!  Some Pascals require file parameters in program heading */

#define MAXSEQ          750000L

#define MAXLINE         150
#define MAXWORD         25

#define VERSION         "NUMSEQ            Version   6/26/2001"


typedef enum {
  T, C, A, G, R, Y, M, W, S, K, D, H, V, B, N
} NUCLEOTIDE;
typedef long SS;

typedef enum {
  ALA, CYS, ASP, GLU, PHE, GLY, HIS, ILE, LYS, LEU, MET, ASN, PRO, GLN, ARG,
  SER, THR, VALINE, TRP, TYR, ASX, GLX, TERM, UNKX
} AMINOACID;
/* Note: VAL is a reserved word in some Pascals */

typedef NUCLEOTIDE NA[15];
typedef NUCLEOTIDE SEQUENCE[MAXSEQ];

typedef Char LETTERS[3];

/* This table holds the amino acid assignments for the 64 codons.*/
typedef AMINOACID GENETICCODE[(long)G - (long)T + 1][(long)G - (long)T + 1]
		  [(long)G - (long)T + 1];

typedef NUCLEOTIDE RULE[3];


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


/* File variables */
Static _TEXT INFILE;   /* input sequence file */
Static _TEXT OUTFILE;   /* output sequence file */
Static _TEXT GCFILE;   /* genetic code file */
Static LINE IFN, OFN, GCFN;   /* input, output, and genetic code filenames*/

/* Sequence variables */
Static SEQUENCE SEQ;
Static long SEQLEN, CHOICE;
Static WORD SEQNAME;   /* Req'd by READSEQ but not used */
Static NA INUC, COMP;
Static SS NUCSET[15];
Static Char NUCHAR[15];

/* Genetic code variables */
Static GENETICCODE GC;
Static RULE RULELIST[64];
Static AMINOACID AALIST[64];
Static long NUMRULES;
Static LETTERS AAS[24];

/* GLOBAL PARAMETERS */
Static long START, FINISH, STARTNO, GROUP, GPL, SENSE, FRAMES, STRANDS, I;
Static Char NUCCASE, COORD, WHICH, FORM, KIND, NUMBERS, NUCS, PEPTIDES;
Static boolean CIRCULAR;
Static LINE HLINE;
Static Char ANSWER;


/* BEGIN MODULE INPLINE */
/*  Read a line from the console.    */
Static Void INPLINE(W_)
LINE *W_;
{
  long I;
  Char CH;
  Char BACKSPACE = '\b';

  /*!!!   if eoln(input) then readln; */
  /*!!!   get(input); */
  W_->LEN = 0;
  while (!P_eoln(stdin)) {
    CH = getchar();
    if (CH == '\n')
      CH = ' ';
    /*!!!       if not eoln then*/
    if (CH == BACKSPACE && W_->LEN > 0) {
      putchar(BACKSPACE);
      W_->LEN--;
    } else if (W_->LEN < MAXLINE) {
      W_->LEN++;
      W_->STR[W_->LEN - 1] = CH;
    }
  }
  scanf("%*[^\n]");
  getchar();
  for (I = W_->LEN; I < MAXLINE; I++)
    W_->STR[I] = ' ';
}  /* INPLINE */


/* END MODULE INPLINE         VERSION= 'SUNMODS     Version  6/26/01'; */

/* BEGIN MODULE WRITELINE */
/*  Write a line to a file using L char, left-justified.  */
Static Void WRITELINE(F, W_, L)
_TEXT *F;
LINE W_;
long L;
{
  long I;

  for (I = 1; I <= L; I++) {
    if (I <= W_.LEN)
      putc(W_.STR[I-1], F->f);
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

/* BEGIN MODULE INPWORD */
/* Read a WORD from the terminal. */
Static Void INPWORD(W_)
WORD *W_;
{
  long I;
  Char CH;
  Char BACKSPACE = '\b';

  W_->LEN = 0;
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
    if (CH == BACKSPACE && W_->LEN > 0) {
      putchar(BACKSPACE);
      W_->LEN--;
    } else if (W_->LEN < MAXWORD) {
      W_->LEN++;
      W_->STR[W_->LEN - 1] = CH;
    }
  }
  /*readln;*/
  for (I = W_->LEN; I < MAXWORD; I++)
    W_->STR[I] = ' ';
}  /* INPWORD */


/* END MODULE INPWORD         VERSION= 'SUNMODS     Version  6/26/01'; */

/* BEGIN MODULE READWORD */
/*  Read a word from a textfile           */
Static Void READWORD(F, W_)
_TEXT *F;
WORD *W_;
{
  long I;
  Char CH;

  W_->LEN = 0;
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
    if (W_->LEN < MAXWORD) {
      W_->LEN++;
      W_->STR[W_->LEN - 1] = getc(F->f);
      if (W_->STR[W_->LEN - 1] == '\n')
	W_->STR[W_->LEN - 1] = ' ';
    } else {
      CH = getc(F->f);
      if (CH == '\n')
	CH = ' ';
    }
  }
  for (I = W_->LEN; I < MAXWORD; I++)
    W_->STR[I] = ' ';
}  /* READWORD */


/* END MODULE READWORD         VERSION= 'SUNMODS     Version  6/26/01'; */

/* BEGIN MODULE WRITEWORD */
/*  Write a word to a file using L char, left-justified.  */
Static Void WRITEWORD(F, W_, L)
_TEXT *F;
WORD W_;
long L;
{
  long I;

  for (I = 1; I <= L; I++) {
    if (I <= W_.LEN)
      putc(W_.STR[I-1], F->f);
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

/* BEGIN MODULE NUC */
/*****************************************************************/
/*  Convert a character to the appropriate nucleotide symbol.    */
/*****************************************************************/
Static NUCLEOTIDE NUC(CH)
Char CH;
{
  NUCLEOTIDE Result;

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

  case 'T':
  case 't':
  case 'U':
  case 'u':
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
  }
  return Result;
}


/* END MODULE NUC         VERSION= 'SUNMODS     Version  6/26/01'; */

/* BEGIN MODULE AA */
/*****************************************************************/
/*  Convert a character to the appropriate amino acid symbol.    */
/*****************************************************************/
Static AMINOACID AA(CH)
Char CH;
{
  AMINOACID Result;

  switch (CH) {

  case 'A':
    Result = ALA;
    break;

  case 'C':
    Result = CYS;
    break;

  case 'D':
    Result = ASP;
    break;

  case 'E':
    Result = GLU;
    break;

  case 'F':
    Result = PHE;
    break;

  case 'G':
    Result = GLY;
    break;

  case 'H':
    Result = HIS;
    break;

  case 'I':
    Result = ILE;
    break;

  case 'K':
    Result = LYS;
    break;

  case 'L':
    Result = LEU;
    break;

  case 'M':
    Result = MET;
    break;

  case 'N':
    Result = ASN;
    break;

  case 'P':
    Result = PRO;
    break;

  case 'Q':
    Result = GLN;
    break;

  case 'R':
    Result = ARG;
    break;

  case 'S':
    Result = SER;
    break;

  case 'T':
    Result = THR;
    break;

  case 'V':
    Result = VALINE;
    break;

  case 'W':
    Result = TRP;
    break;

  case 'Y':
    Result = TYR;
    break;

  case 'B':
    Result = ASX;
    break;

  case 'Z':
    Result = GLX;
    break;

  case 'X':
    Result = UNKX;
    break;

  case '*':
    Result = TERM;
    break;
  }
  return Result;
}  /* AA */


/* END MODULE AA         VERSION= 'SUNMODS     Version  6/26/01'; */

/*******************************************************/
/* Initialize tables of nucleotides, amino acids etc.  */
/* This is done only once when the program begins.     */
/*******************************************************/
Static Void INITTABLES()
{
  /*   Initialize INUC and COMP, the arrays which */
  /*   hold the  values of the NUCLEOTIDEs.      */
  INUC[(long)T] = T;
  COMP[(long)T] = A;
  INUC[(long)C] = C;
  COMP[(long)C] = G;
  INUC[(long)A] = A;
  COMP[(long)A] = T;
  INUC[(long)G] = G;
  COMP[(long)G] = C;
  INUC[(long)R] = R;
  COMP[(long)R] = Y;
  INUC[(long)D] = D;
  COMP[(long)D] = H;
  INUC[(long)V] = V;
  COMP[(long)V] = B;
  INUC[(long)M] = M;
  COMP[(long)M] = K;
  INUC[(long)K] = K;
  COMP[(long)K] = M;
  INUC[(long)B] = B;
  COMP[(long)B] = V;
  INUC[(long)H] = H;
  COMP[(long)H] = D;
  INUC[(long)Y] = Y;
  COMP[(long)Y] = R;
  INUC[(long)W] = W;
  COMP[(long)W] = W;
  INUC[(long)S] = S;
  COMP[(long)S] = S;
  INUC[(long)N] = N;
  COMP[(long)N] = N;

  /* Initialize the nucleotide sets */
  NUCSET[(long)N] = (1L << ((long)N + 1)) - (1L << ((long)T));
  NUCSET[(long)R] = (1L << ((long)A)) | (1L << ((long)G)) | (1L << ((long)R));
  NUCSET[(long)Y] = (1L << ((long)C)) | (1L << ((long)T)) | (1L << ((long)Y));
  NUCSET[(long)M] = (1L << ((long)A)) | (1L << ((long)C)) | (1L << ((long)M));
  NUCSET[(long)W] = (1L << ((long)A)) | (1L << ((long)T)) | (1L << ((long)W));
  NUCSET[(long)S] = (1L << ((long)C)) | (1L << ((long)G)) | (1L << ((long)S));
  NUCSET[(long)K] = (1L << ((long)G)) | (1L << ((long)T)) | (1L << ((long)K));
  NUCSET[(long)D] = (1L << ((long)T)) | (1L << ((long)A)) | (1L << ((long)G)) |
		    (1L << ((long)W)) | (1L << ((long)K)) | (1L << ((long)R));
  NUCSET[(long)H] = (1L << ((long)T)) | (1L << ((long)C)) | (1L << ((long)A)) |
		    (1L << ((long)Y)) | (1L << ((long)W)) | (1L << ((long)M));
  NUCSET[(long)V] = (1L << ((long)C)) | (1L << ((long)A)) | (1L << ((long)G)) |
		    (1L << ((long)M)) | (1L << ((long)S)) | (1L << ((long)R));
  NUCSET[(long)B] = (1L << ((long)T)) | (1L << ((long)C)) | (1L << ((long)G)) |
		    (1L << ((long)Y)) | (1L << ((long)K)) | (1L << ((long)S));
  NUCSET[(long)T] = 1L << ((long)T);
  NUCSET[(long)C] = 1L << ((long)C);
  NUCSET[(long)A] = 1L << ((long)A);
  NUCSET[(long)G] = 1L << ((long)G);

  /* Initialize GC using universal genetic code */
  GC[0][0][0] = PHE;
  GC[0][(long)C - (long)T][0] = SER;
  GC[0][(long)A - (long)T][0] = TYR;
  GC[0][(long)G - (long)T][0] = CYS;
  GC[0][0][(long)C - (long)T] = PHE;
  GC[0][(long)C - (long)T][(long)C - (long)T] = SER;
  GC[0][(long)A - (long)T][(long)C - (long)T] = TYR;
  GC[0][(long)G - (long)T][(long)C - (long)T] = CYS;
  GC[0][0][(long)A - (long)T] = LEU;
  GC[0][(long)C - (long)T][(long)A - (long)T] = SER;
  GC[0][(long)A - (long)T][(long)A - (long)T] = TERM;
  GC[0][(long)G - (long)T][(long)A - (long)T] = TERM;
  GC[0][0][(long)G - (long)T] = LEU;
  GC[0][(long)C - (long)T][(long)G - (long)T] = SER;
  GC[0][(long)A - (long)T][(long)G - (long)T] = TERM;
  GC[0][(long)G - (long)T][(long)G - (long)T] = TRP;
  GC[(long)C - (long)T][0][0] = LEU;
  GC[(long)C - (long)T][(long)C - (long)T][0] = PRO;
  GC[(long)C - (long)T][(long)A - (long)T][0] = HIS;
  GC[(long)C - (long)T][(long)G - (long)T][0] = ARG;
  GC[(long)C - (long)T][0][(long)C - (long)T] = LEU;
  GC[(long)C - (long)T][(long)C - (long)T][(long)C - (long)T] = PRO;
  GC[(long)C - (long)T][(long)A - (long)T][(long)C - (long)T] = HIS;
  GC[(long)C - (long)T][(long)G - (long)T][(long)C - (long)T] = ARG;
  GC[(long)C - (long)T][0][(long)A - (long)T] = LEU;
  GC[(long)C - (long)T][(long)C - (long)T][(long)A - (long)T] = PRO;
  GC[(long)C - (long)T][(long)A - (long)T][(long)A - (long)T] = GLN;
  GC[(long)C - (long)T][(long)G - (long)T][(long)A - (long)T] = ARG;
  GC[(long)C - (long)T][0][(long)G - (long)T] = LEU;
  GC[(long)C - (long)T][(long)C - (long)T][(long)G - (long)T] = PRO;
  GC[(long)C - (long)T][(long)A - (long)T][(long)G - (long)T] = GLN;
  GC[(long)C - (long)T][(long)G - (long)T][(long)G - (long)T] = ARG;
  GC[(long)A - (long)T][0][0] = ILE;
  GC[(long)A - (long)T][(long)C - (long)T][0] = THR;
  GC[(long)A - (long)T][(long)A - (long)T][0] = ASN;
  GC[(long)A - (long)T][(long)G - (long)T][0] = SER;
  GC[(long)A - (long)T][0][(long)C - (long)T] = ILE;
  GC[(long)A - (long)T][(long)C - (long)T][(long)C - (long)T] = THR;
  GC[(long)A - (long)T][(long)A - (long)T][(long)C - (long)T] = ASN;
  GC[(long)A - (long)T][(long)G - (long)T][(long)C - (long)T] = SER;
  GC[(long)A - (long)T][0][(long)A - (long)T] = ILE;
  GC[(long)A - (long)T][(long)C - (long)T][(long)A - (long)T] = THR;
  GC[(long)A - (long)T][(long)A - (long)T][(long)A - (long)T] = LYS;
  GC[(long)A - (long)T][(long)G - (long)T][(long)A - (long)T] = ARG;
  GC[(long)A - (long)T][0][(long)G - (long)T] = MET;
  GC[(long)A - (long)T][(long)C - (long)T][(long)G - (long)T] = THR;
  GC[(long)A - (long)T][(long)A - (long)T][(long)G - (long)T] = LYS;
  GC[(long)A - (long)T][(long)G - (long)T][(long)G - (long)T] = ARG;
  GC[(long)G - (long)T][0][0] = VALINE;
  GC[(long)G - (long)T][(long)C - (long)T][0] = ALA;
  GC[(long)G - (long)T][(long)A - (long)T][0] = ASP;
  GC[(long)G - (long)T][(long)G - (long)T][0] = GLY;
  GC[(long)G - (long)T][0][(long)C - (long)T] = VALINE;
  GC[(long)G - (long)T][(long)C - (long)T][(long)C - (long)T] = ALA;
  GC[(long)G - (long)T][(long)A - (long)T][(long)C - (long)T] = ASP;
  GC[(long)G - (long)T][(long)G - (long)T][(long)C - (long)T] = GLY;
  GC[(long)G - (long)T][0][(long)A - (long)T] = VALINE;
  GC[(long)G - (long)T][(long)C - (long)T][(long)A - (long)T] = ALA;
  GC[(long)G - (long)T][(long)A - (long)T][(long)A - (long)T] = GLU;
  GC[(long)G - (long)T][(long)G - (long)T][(long)A - (long)T] = GLY;
  GC[(long)G - (long)T][0][(long)G - (long)T] = VALINE;
  GC[(long)G - (long)T][(long)C - (long)T][(long)G - (long)T] = ALA;
  GC[(long)G - (long)T][(long)A - (long)T][(long)G - (long)T] = GLU;
  GC[(long)G - (long)T][(long)G - (long)T][(long)G - (long)T] = GLY;
}  /* INITTABLES */


/* **************************************************** */
/* Initialize parameters for a new sequence.            */
/* **************************************************** */
Static Void INITPARAM()
{
  /* Set default values for global parameters */
  START = 1;
  FINISH = SEQLEN;
  STARTNO = START;
  SENSE = 1;
  GROUP = 10;
  GPL = 7;
  NUCCASE = 'U';
  FRAMES = 1;
  FORM = 'L';
  COORD = 'S';
  WHICH = 'I';
  STRANDS = 1;
  KIND = 'D';
  NUMBERS = 'Y';
  NUCS = 'Y';
  PEPTIDES = 'N';
}  /* INITPARAM */


/* ************************************************************ */
/* Initialize nucleotide and amino acid output strings          */
/* according to parameters.                                     */
/* ************************************************************ */
Static Void INITSTRING()
{
  /* Initialize nucleotide output characters */
  switch (NUCCASE) {

  case 'U':
    if (KIND == 'D')
      NUCHAR[(long)T] = 'T';
    else
      NUCHAR[(long)T] = 'U';
    NUCHAR[(long)C] = 'C';
    NUCHAR[(long)A] = 'A';
    NUCHAR[(long)G] = 'G';
    NUCHAR[(long)R] = 'R';
    NUCHAR[(long)D] = 'D';
    NUCHAR[(long)V] = 'V';
    NUCHAR[(long)M] = 'M';
    NUCHAR[(long)K] = 'K';
    NUCHAR[(long)B] = 'B';
    NUCHAR[(long)H] = 'H';
    NUCHAR[(long)Y] = 'Y';
    NUCHAR[(long)W] = 'W';
    NUCHAR[(long)S] = 'S';
    NUCHAR[(long)N] = 'N';
    break;

  case 'L':
    if (KIND == 'D')
      NUCHAR[(long)T] = 't';
    else
      NUCHAR[(long)T] = 'u';
    NUCHAR[(long)C] = 'c';
    NUCHAR[(long)A] = 'a';
    NUCHAR[(long)G] = 'g';
    NUCHAR[(long)R] = 'r';
    NUCHAR[(long)D] = 'd';
    NUCHAR[(long)V] = 'v';
    NUCHAR[(long)M] = 'm';
    NUCHAR[(long)K] = 'k';
    NUCHAR[(long)B] = 'b';
    NUCHAR[(long)H] = 'h';
    NUCHAR[(long)Y] = 'y';
    NUCHAR[(long)W] = 'w';
    NUCHAR[(long)S] = 's';
    NUCHAR[(long)N] = 'n';
    break;
  }/* case NUCCASE */

  /* Initialize amino acid output strings. */
  switch (FORM) {

  case 'S':
    memcpy(AAS[(long)PHE], "F  ", sizeof(LETTERS));
    memcpy(AAS[(long)LEU], "L  ", sizeof(LETTERS));
    memcpy(AAS[(long)ILE], "I  ", sizeof(LETTERS));
    memcpy(AAS[(long)MET], "M  ", sizeof(LETTERS));
    memcpy(AAS[(long)VALINE], "V  ", sizeof(LETTERS));
    memcpy(AAS[(long)SER], "S  ", sizeof(LETTERS));
    memcpy(AAS[(long)PRO], "P  ", sizeof(LETTERS));
    memcpy(AAS[(long)THR], "T  ", sizeof(LETTERS));
    memcpy(AAS[(long)ALA], "A  ", sizeof(LETTERS));
    memcpy(AAS[(long)TYR], "Y  ", sizeof(LETTERS));
    memcpy(AAS[(long)HIS], "H  ", sizeof(LETTERS));
    memcpy(AAS[(long)GLN], "Q  ", sizeof(LETTERS));
    memcpy(AAS[(long)ASN], "N  ", sizeof(LETTERS));
    memcpy(AAS[(long)LYS], "K  ", sizeof(LETTERS));
    memcpy(AAS[(long)ASP], "D  ", sizeof(LETTERS));
    memcpy(AAS[(long)GLU], "E  ", sizeof(LETTERS));
    memcpy(AAS[(long)CYS], "C  ", sizeof(LETTERS));
    memcpy(AAS[(long)TRP], "W  ", sizeof(LETTERS));
    memcpy(AAS[(long)ARG], "R  ", sizeof(LETTERS));
    memcpy(AAS[(long)GLY], "G  ", sizeof(LETTERS));
    memcpy(AAS[(long)ASX], "B  ", sizeof(LETTERS));
    memcpy(AAS[(long)GLX], "Z  ", sizeof(LETTERS));
    memcpy(AAS[(long)TERM], "*  ", sizeof(LETTERS));
    memcpy(AAS[(long)UNKX], "X  ", sizeof(LETTERS));
    break;

  case 'L':
    memcpy(AAS[(long)PHE], "Phe", sizeof(LETTERS));
    memcpy(AAS[(long)LEU], "Leu", sizeof(LETTERS));
    memcpy(AAS[(long)ILE], "Ile", sizeof(LETTERS));
    memcpy(AAS[(long)MET], "MET", sizeof(LETTERS));
    memcpy(AAS[(long)VALINE], "Val", sizeof(LETTERS));
    memcpy(AAS[(long)SER], "Ser", sizeof(LETTERS));
    memcpy(AAS[(long)PRO], "Pro", sizeof(LETTERS));
    memcpy(AAS[(long)THR], "Thr", sizeof(LETTERS));
    memcpy(AAS[(long)ALA], "Ala", sizeof(LETTERS));
    memcpy(AAS[(long)TYR], "Tyr", sizeof(LETTERS));
    memcpy(AAS[(long)HIS], "His", sizeof(LETTERS));
    memcpy(AAS[(long)GLN], "Gln", sizeof(LETTERS));
    memcpy(AAS[(long)ASN], "Asn", sizeof(LETTERS));
    memcpy(AAS[(long)LYS], "Lys", sizeof(LETTERS));
    memcpy(AAS[(long)ASP], "Asp", sizeof(LETTERS));
    memcpy(AAS[(long)GLU], "Glu", sizeof(LETTERS));
    memcpy(AAS[(long)CYS], "Cys", sizeof(LETTERS));
    memcpy(AAS[(long)TRP], "Trp", sizeof(LETTERS));
    memcpy(AAS[(long)ARG], "Arg", sizeof(LETTERS));
    memcpy(AAS[(long)GLY], "Gly", sizeof(LETTERS));
    memcpy(AAS[(long)ASX], "Asx", sizeof(LETTERS));
    memcpy(AAS[(long)GLX], "Glx", sizeof(LETTERS));
    memcpy(AAS[(long)TERM], "   ", sizeof(LETTERS));
    memcpy(AAS[(long)UNKX], "---", sizeof(LETTERS));
    break;
  }/* case FORM */
}  /*INITSTRING*/


/*************************************************************/
/*  Read in an alternative genetic code from a datafile.     */
/*  The datafile takes the form of a two-dimensional table   */
/*  of amino acids as found in Watson, J.D. (1976), Molecular*/
/*  Biology of the Gene pg.356, Table 13-7. The file         */
/*  should contain a single CAPITAL letter amino acid code   */
/*  for each position corresponding to a given codon.  This  */
/*  program will read in the first 64 characters that could  */
/*  be amino acid or stop ( * ) symbols. Blanks or lowercase */
/*  letters may be included for clarity.                     */
/*************************************************************/
Static Void READCODE(F, GCFN, GC)
_TEXT *F;
LINE *GCFN;
AMINOACID (*GC)[(long)G - (long)T + 1][(long)G - (long)T + 1];
{
  NUCLEOTIDE N1, N2, N3;
  Char AACHAR;

  printf("Type name of genetic code file:\n");
  GETFILE(F, 'I', GCFN);
  for (N1 = T; (long)N1 <= (long)G; N1 = (NUCLEOTIDE)((long)N1 + 1)) {
    for (N3 = T; (long)N3 <= (long)G; N3 = (NUCLEOTIDE)((long)N3 + 1))
    {   /* This is correct. Think about it. */
      for (N2 = T; (long)N2 <= (long)G; N2 = (NUCLEOTIDE)((long)N2 + 1)) {
	AACHAR = ' ';
	while (!BUFEOF(F->f) && AACHAR != '*' && AACHAR != 'Y' &&
	       AACHAR != 'W' && AACHAR != 'V' && AACHAR != 'T' &&
	       AACHAR != 'S' && AACHAR != 'R' && AACHAR != 'Q' &&
	       AACHAR != 'P' && AACHAR != 'N' && AACHAR != 'M' &&
	       AACHAR != 'L' && AACHAR != 'K' && AACHAR != 'I' &&
	       AACHAR != 'H' && AACHAR != 'G' && AACHAR != 'F' &&
	       AACHAR != 'E' && AACHAR != 'D' && AACHAR != 'C' &&
	       AACHAR != 'A') {
	  AACHAR = getc(F->f);
	  if (AACHAR == '\n')
	    AACHAR = ' ';
	}
	GC[(long)N1 - (long)T][(long)N2 - (long)T]
	  [(long)N3 - (long)T] = AA(AACHAR);
      }
    }
  }
  /*!!!*/
  if (F->f != NULL)
    fclose(F->f);
  F->f = NULL;
}  /* READCODE */


/* Given a set, return the nucleotide symbol corresponding to the minimal
   set which contains all members of that set */
Local NUCLEOTIDE NUCSYM(NS)
SS NS;
{
  NUCLEOTIDE N1 = T;

  while ((NS & (~NUCSET[(long)N1])) != 0)
    N1 = (NUCLEOTIDE)((long)N1 + 1);
  return N1;
}  /* NUCSYM */

/* Add a rule to the end of RULELIST */
Local Void ADDRULE(N1, N2, N3, AMA)
NUCLEOTIDE N1, N2, N3;
AMINOACID AMA;
{
  NUMRULES++;
  RULELIST[NUMRULES-1][0] = N1;
  RULELIST[NUMRULES-1][1] = N2;
  RULELIST[NUMRULES-1][2] = N3;
  AALIST[NUMRULES-1] = AMA;
}  /* ADDRULE */

/* Delete a rule from the RULELIST. Move the other rules down the list.*/
Local Void DELRULE(I)
long I;
{
  long J;

  J = I + 1;
  while (J <= NUMRULES) {
    memcpy(RULELIST[I-1], RULELIST[J-1], sizeof(RULE));
    AALIST[I-1] = AALIST[J-1];
    I++;
    J++;
  }
  NUMRULES--;
}  /* DELRULE */


/*************************************************************/
/* Using the 64 unambiguous codons stored in GC, make the    */
/* the optimal list of ambiguity rules, assigning codons to  */
/* amino acids. All possible ambiguous codons are taken into */
/* account.                                                  */
/*************************************************************/
Static Void MAKERULES(GC)
AMINOACID (*GC)[(long)G - (long)T + 1][(long)G - (long)T + 1];
{
  NUCLEOTIDE N1, N2, N3;   /* current nucleotides */
  AMINOACID AMA;   /* current amino acid */
  long RA, RB, FIRST, LASTOLDRULE;   /* All are indices of rules */
  SS S1, S2, S3;

  /* For each amino acid, make the minimal set of rules such that all
     possible codons for that amino acid are included in some rule */
  NUMRULES = 0;
  for (AMA = ALA; (long)AMA <= (long)TERM; AMA = (AMINOACID)((long)AMA + 1)) {
    FIRST = NUMRULES + 1;

    /* STEP 1: Add to RULELIST all codons for a given amino acid */
    for (N2 = T; (long)N2 <= (long)G; N2 = (NUCLEOTIDE)((long)N2 + 1))
    {   /* This is correct, too. */
      for (N1 = T; (long)N1 <= (long)G; N1 = (NUCLEOTIDE)((long)N1 + 1)) {
	for (N3 = T; (long)N3 <= (long)G; N3 = (NUCLEOTIDE)((long)N3 + 1)) {
	  if (GC[(long)N1 - (long)T][(long)N2 - (long)T]
	      [(long)N3 - (long)T] == AMA)
	    ADDRULE(N1, N2, N3, AMA);
	}
      }
    }

    /* STEP2: Generalization of the rules.  If any two positions in rule A
       are both equal to the corresponding positions in rule B, then update
       the remaining position of rule A to include both members at that
       position. Delete rule B. */
    RA = FIRST;
    while (RA < NUMRULES) {
      RB = RA + 1;
      while (RB <= NUMRULES) {
	if (RULELIST[RA-1][0] == RULELIST[RB-1][0] &&
	    RULELIST[RA-1][1] == RULELIST[RB-1][1]) {
	  RULELIST[RA-1][2] = NUCSYM(NUCSET[(long)RULELIST[RA-1][2]] |
				     NUCSET[(long)RULELIST[RB-1][2]]);
	  DELRULE(RB);
	  continue;
	}
	if (RULELIST[RA-1][0] == RULELIST[RB-1][0] &&
	    RULELIST[RA-1][2] == RULELIST[RB-1][2]) {
	  RULELIST[RA-1][1] = NUCSYM(NUCSET[(long)RULELIST[RA-1][1]] |
				     NUCSET[(long)RULELIST[RB-1][1]]);
	  DELRULE(RB);
	} else if (RULELIST[RA-1][1] == RULELIST[RB-1][1] &&
		   RULELIST[RA-1][2] == RULELIST[RB-1][2]) {
	  RULELIST[RA-1][0] = NUCSYM(NUCSET[(long)RULELIST[RA-1][0]] |
				     NUCSET[(long)RULELIST[RB-1][0]]);
	  DELRULE(RB);
	} else
	  RB++;
      }
      RA++;
    }

    /* STEP 3: Generalization of the rules.  If any two positions in rule A
       both intersect with the corresponding positions in rule B, then
       create a new rule comprised of the intersection at those two
       positions and the sum of the remaining positions of rule A and
       rule B. */
    RA = FIRST;
    LASTOLDRULE = NUMRULES;   /* The last rule made in STEP 2 */
    while (RA < LASTOLDRULE) {
      RB = RA + 1;
      while (RB <= LASTOLDRULE) {
	S1 = NUCSET[(long)RULELIST[RA-1][0]] & NUCSET[(long)RULELIST[RB-1][0]];
	S2 = NUCSET[(long)RULELIST[RA-1][1]] & NUCSET[(long)RULELIST[RB-1][1]];
	S3 = NUCSET[(long)RULELIST[RA-1][2]] & NUCSET[(long)RULELIST[RB-1][2]];
	if (S1 != 0 && S2 != 0)
	  ADDRULE(NUCSYM(S1), NUCSYM(S2),
	    NUCSYM(NUCSET[(long)RULELIST[RA-1][2]] | NUCSET[(long)RULELIST[RB-1]
							    [2]]), AMA);
	if (S1 != 0 && S3 != 0)
	  ADDRULE(NUCSYM(S1),
		  NUCSYM(NUCSET[(long)RULELIST[RA-1][1]] |
			 NUCSET[(long)RULELIST[RB-1][1]]), NUCSYM(S3), AMA);
	if (S2 != 0 && S3 != 0)
	  ADDRULE(NUCSYM(NUCSET[(long)RULELIST[RA-1][0]] |
			 NUCSET[(long)RULELIST[RB-1][0]]), NUCSYM(S2),
		  NUCSYM(S3), AMA);
	RB++;
      }
      RA++;
    }

    /* STEP 4: Delete any rule which is entirely contained within
       another rule */
    RA = FIRST;
    while (RA < NUMRULES) {
      RB = RA + 1;
      while (RB <= NUMRULES) {
	if (((1L << ((long)RULELIST[RA-1][0])) & NUCSET[(long)RULELIST[RB-1]
							[0]]) != 0 &&
	    ((1L << ((long)RULELIST[RA-1][1])) & NUCSET[(long)RULELIST[RB-1]
							[1]]) != 0 &&
	    ((1L << ((long)RULELIST[RA-1][2])) & NUCSET[(long)RULELIST[RB-1]
							[2]]) != 0)
	  DELRULE(RA);
	else
	  RB++;
      }
      RA++;
    }

  }  /* for AMA */

  /* STEP 5: Add a rule for unknown codons */
  NUMRULES++;
  RULELIST[NUMRULES-1][0] = N;
  RULELIST[NUMRULES-1][1] = N;
  RULELIST[NUMRULES-1][2] = N;
  AALIST[NUMRULES-1] = UNKX;

  /*!!!*/
  /* Print rule list */
  /*    FORM:='L';INITSTRING;
        writeln('The following rules implement the genetic code:');
        for RA:= 1 to NUMRULES do begin
          for RB:= 1 to 3 do write(NUCHAR[RULELIST[RA][RB]]);
          write(' ':7);
          writeln(AAS[AALIST[RA]]:10)
          end;
        writeln('Press RETURN to continue:');
        readln */
}  /* MAKERULES */


/* BEGIN MODULE SAMEWORD */
/* Compare two WORDS for equality */
Static boolean SAMEWORD(W1, W2)
WORD *W1, *W2;
{
  long I;
  boolean T_;

  if (W1->LEN == W2->LEN) {
    T_ = true;
    I = 1;
    while (I <= W1->LEN && T_) {
      if (W1->STR[I-1] == W2->STR[I-1])
	I++;
      else
	T_ = false;
    }
    return T_;
  } else
    return false;
}  /* SAMEWORD */


/* END MODULE SAMEWORD         VERSION= 'SUNMODS     Version  6/26/01'; */

/* BEGIN MODULE READSEQ */
/*  ******************************************* */
/*  Read a DNA or RNA sequence from SFILE       */
/*   and store it in S.                         */
/*  ******************************************* */
Static Void READSEQ(SFILE, S_, SEQLEN, NAME, CIRCULAR)
_TEXT *SFILE;
NUCLEOTIDE *S_;
long *SEQLEN;
WORD *NAME;
boolean *CIRCULAR;
{
  Char CH, FILETYPE, LASTNUM;
  boolean TOOLONG = false, BADFILE = false;
  WORD ID, LOCUS, ORI, CIRC;

  /* Prompt user for sequence file type */
  printf("The following file formats can be read:\n");
  printf("  F:free format   B:BIONET    G:GENBANK\n");
  do {
    printf("Type letter of format (F|B|G)\n");
    scanf("%c%*[^\n]", &CH);
    getchar();
    if (CH == '\n')
      CH = ' ';
  } while (CH != 'g' && CH != 'G' && CH != 'b' && CH != 'B' && CH != 'f' &&
	   CH != 'F');
  switch (CH) {

  case 'F':
  case 'f':
    FILETYPE = 'F';
    break;

  case 'B':
  case 'b':
    FILETYPE = 'B';
    break;

  case 'G':
  case 'g':
    FILETYPE = 'G';
    break;
  }
  printf("Reading input file...\n");
  NAME->LEN = 0;

  /* For BIONET or GENBANK, read in sequence name and topology */
  /* Advance to beginning of sequence */
  if (FILETYPE == 'B') {
    /* First non-comment line is name. Name may be blank.*/
    while (P_peek(SFILE->f) == ';') {
      fscanf(SFILE->f, "%*[^\n]");
      getc(SFILE->f);
    }
    while ((P_peek(SFILE->f) == ' ') & (!P_eoln(SFILE->f)))
      getc(SFILE->f);
    if (!P_eoln(SFILE->f))
      READWORD(SFILE, NAME);
    else
      NAME->LEN = 0;
    if (!BUFEOF(SFILE->f)) {
      fscanf(SFILE->f, "%*[^\n]");
      getc(SFILE->f);
    }
  }  /* BIONET */

  else if (FILETYPE == 'G') {
    /* Initialize identifiers */
    LOCUS.STR[0] = 'L';
    LOCUS.STR[1] = 'O';
    LOCUS.STR[2] = 'C';
    LOCUS.STR[3] = 'U';
    LOCUS.STR[4] = 'S';
    LOCUS.LEN = 5;
    ORI.STR[0] = 'O';
    ORI.STR[1] = 'R';
    ORI.STR[2] = 'I';
    ORI.STR[3] = 'G';
    ORI.STR[4] = 'I';
    ORI.STR[5] = 'N';
    ORI.LEN = 6;
    CIRC.STR[0] = 'c';
    CIRC.STR[1] = 'i';
    CIRC.STR[2] = 'r';
    CIRC.STR[3] = 'c';
    CIRC.STR[4] = 'u';
    CIRC.STR[5] = 'l';
    CIRC.STR[6] = 'a';
    CIRC.STR[7] = 'r';
    CIRC.LEN = 8;
    /* Advance to LOCUS line. Read in NAME and topology. */
    while (!((P_peek(SFILE->f) == 'L') | BUFEOF(SFILE->f))) {
      fscanf(SFILE->f, "%*[^\n]");
      getc(SFILE->f);
    }
    if (!BUFEOF(SFILE->f)) {
      READWORD(SFILE, &ID);
      if (SAMEWORD(&ID, &LOCUS)) {
	if (!BUFEOF(SFILE->f))
	  READWORD(SFILE, NAME);
	if (!BUFEOF(SFILE->f))   /* skip seq. length */
	  READWORD(SFILE, &ID);
	if (!BUFEOF(SFILE->f))   /* skip. 'bp' */
	  READWORD(SFILE, &ID);
	/* After 'bp', there is an optional field telling
	   the type of molecule (ss-RNA, ds-DNA etc.)
	   Since this field is optional, we must test the
	   next two tokens to see if they are 'circular' */
	*CIRCULAR = false;
	if (!BUFEOF(SFILE->f))
	  READWORD(SFILE, &ID);
	if (SAMEWORD(&ID, &CIRC))
	  *CIRCULAR = true;
	else {
	  if (!BUFEOF(SFILE->f))
	    READWORD(SFILE, &ID);
	  if (SAMEWORD(&ID, &CIRC))
	    *CIRCULAR = true;
	}
      }  /* SAMEWORD(ID,LOCUS) */
      else
	BADFILE = true;
    }

    /* Advance to ORIGIN line. Sequence begins on next line */
    if (!BUFEOF(SFILE->f)) {
      do {
	fscanf(SFILE->f, "%*[^\n]");
	getc(SFILE->f);
	if (P_peek(SFILE->f) == 'O')
	  READWORD(SFILE, &ID);
      } while (!(SAMEWORD(&ID, &ORI) | BUFEOF(SFILE->f)));
      if (SAMEWORD(&ID, &ORI)) {
	fscanf(SFILE->f, "%*[^\n]");
	getc(SFILE->f);
      }
    }
    if (BUFEOF(SFILE->f))
      BADFILE = true;
  }  /* GENBANK */

  /* Read in sequence */
  *SEQLEN = 0;
  if (BADFILE) {
    printf(">>> ERROR: Not a GENBANK file. No sequence read.\n");
    /*!!!   CLOSE(SFILE);*/
    return;
  }
  while (!BUFEOF(SFILE->f)) {
    while (!P_eoln(SFILE->f)) {
      CH = getc(SFILE->f);
      if (CH == '\n')
	CH = ' ';
      if (CH == 'b' || CH == 'B' || CH == 'v' || CH == 'V' || CH == 'h' ||
	  CH == 'H' || CH == 'd' || CH == 'D' || CH == 'k' || CH == 'K' ||
	  CH == 's' || CH == 'S' || CH == 'w' || CH == 'W' || CH == 'm' ||
	  CH == 'M' || CH == 'y' || CH == 'Y' || CH == 'r' || CH == 'R' ||
	  CH == 'n' || CH == 'N' || CH == 'u' || CH == 'U' || CH == 't' ||
	  CH == 'T' || CH == 'g' || CH == 'G' || CH == 'c' || CH == 'C' ||
	  CH == 'a' || CH == 'A') {
/* p2c: numseq.p, line 805: 
 * Note: Line breaker spent 0.0 seconds, 5000 tries on line 1323 [251] */
	if (*SEQLEN < MAXSEQ - 2) {
	  /*write(CH);*/
	  (*SEQLEN)++;
	  S_[*SEQLEN - 1] = NUC(CH);
	} else
	  TOOLONG = true;
	continue;
      }
      if (CH == ';') {   /*begin    comment in input file */
	fscanf(SFILE->f, "%*[^\n]");
	getc(SFILE->f);
      }
      /*;writeln end*/
      else if (CH == '2' || CH == '1')
	LASTNUM = CH;
    }
    fscanf(SFILE->f, "%*[^\n]");
    getc(SFILE->f);   /*writeln*/
  }
  if (TOOLONG)
    printf(">>> WARNING! Sequence length exceeds MAXSEQ-2. Seq. truncated.\n");

  if (FILETYPE == 'F') {
    do {
      printf("Is sequence circular or linear? (Type C or L)\n");
      scanf("%c%*[^\n]", &CH);
      getchar();
      if (CH == '\n')
	CH = ' ';
    } while (CH != 'l' && CH != 'L' && CH != 'c' && CH != 'C');
    switch (CH) {

    case 'C':
    case 'c':
      *CIRCULAR = true;
      break;

    case 'L':
    case 'l':
      *CIRCULAR = false;
      break;
    }
    return;
  }
  if (FILETYPE != 'B')
    return;
  if (LASTNUM == '1')
    *CIRCULAR = false;
  else
    *CIRCULAR = true;
}  /* READSEQ */


typedef Char LETTERS_[10];
typedef long CHSET[4];


/*  Read an integer parameter from the console and check */
/*    that it is in range.                               */
Local Void GETNUMBER(P, PNAME, LOW, HIGH)
long *P;
Char *PNAME;
long LOW, HIGH;
{
  printf("Type new value for %.10s  (CURRENT VALUE: %12ld)\n", PNAME, *P);
  GETINTEGER(P, LOW, HIGH);
  scanf("%*[^\n]");
  getchar();
}  /* GETNUMBER */

Local Char TOUPPER(CH)
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

/* BEGIN MODULE GETCHAR */
/* Read a character from the console and check */
/*  for correct response.                      */
Local Void GETCHAR(CH, PNAME, ALLOWED)
Char *CH;
Char *PNAME;
long *ALLOWED;
{
  /* BEGIN MODULE TOUPPER */
  /* Change a character from lower to uppercase */
  /* END MODULE TOUPPER         VERSION= 'SUNMODS     Version  6/26/01'; */
  putchar('\n');
  do {
    printf("Type new value for %.10s  (CURRENT VALUE: %c)\n", PNAME, *CH);
    scanf("%c%*[^\n]", CH);
    getchar();
    if (*CH == '\n')
      *CH = ' ';
    *CH = TOUPPER(*CH);
    if (!P_inset(*CH, ALLOWED))
      printf("Inappropriate response: %c\n", *CH);
  } while (!P_inset(*CH, ALLOWED));   /* GETCHAR */
}

/* END MODULE GETCHAR         VERSION= 'SUNMODS     Version  6/26/01'; */

/* Display  parameters on screen */
Local Void DISPLAY()
{
  _TEXT TEMP;

  printf("Name: ");
  TEMP.f = stdout;
  *TEMP.name = '\0';
  WRITEWORD(&TEMP, SEQNAME, 20L);
  switch (CIRCULAR) {

  case true:
    printf("Topology:   CIRCULAR");
    break;

  case false:
    printf("Topology:     LINEAR");
    break;
  }
  printf("%14s%11ld nt\n", "Length: ", SEQLEN);
  TEMP.f = stdout;
  *TEMP.name = '\0';
  WRITELINE(&TEMP, HLINE, 80L);
  printf("\n%12cParameter   Description/Response                 Value\n",
	 ' ');
  TEMP.f = stdout;
  *TEMP.name = '\0';
  WRITELINE(&TEMP, HLINE, 80L);
  printf("\n%12c 1)START    first nucleotide printed%16ld\n", ' ', START);
  printf("%12c 2)FINISH   last  nucleotide printed%16ld\n", ' ', FINISH);
  printf("%12c 3)NUCCASE   U:(A,G,C,T...), l:(a,g,c,t...)%9c\n", ' ', NUCCASE);
  printf("%12c 4)STARTNO  number of starting nucleotide%11ld\n", ' ', STARTNO);
  printf("%12c 5)GROUP    number every GROUP nucleotides%10ld\n", ' ', GROUP);
  printf("%12c 6)GPL      number of GROUPs printed per line%7ld\n", ' ', GPL);
  printf("%12c 7)WHICH    I: input strand  O: opposite strand%5c\n",
	 ' ', WHICH);
  printf("%12c 8)STRANDS  1: one  strand,  2:both strands%9ld\n",
	 ' ', STRANDS);
  printf("%12c 9)KIND     R:RNA            D:DNA%18c\n", ' ', KIND);
  printf("%12c10)NUMBERS  Number  the sequence    (Y or N)%8c\n",
	 ' ', NUMBERS);
  printf("%12c11)NUCS     Print nucleotide seq.   (Y or N)%8c\n", ' ', NUCS);
  printf("%12c12)PEPTIDES Print amino acid seq.   (Y or N)%8c\n",
	 ' ', PEPTIDES);
  printf("%12c13)FRAMES   1 for this frame, 3 for 3 frames%8ld\n",
	 ' ', FRAMES);
  printf("%12c14)FORM     L:3 letter amino acid, S: 1 letter%6c\n", ' ', FORM);
  TEMP.f = stdout;
  *TEMP.name = '\0';
  WRITELINE(&TEMP, HLINE, 80L);
  putchar('\n');
}  /* DISPLAY */


/* END MODULE READSEQ         VERSION= 'SUNMODS     Version  6/26/01'; */

/* **************************************************** */
/* Prompt user for parameters used by program.          */
/* **************************************************** */
Static Void PARAMETERS()
{
  long RESPONSE;
  long SET[4];
  long SET1[4];
  long SET2[4];
  long SET3[4];
  long SET4[4];

  /* Prompt user for new parameter values */
  do {
    fprintf(stdout, "\f");
    DISPLAY();
    if (WHICH == 'O')
      printf("Be sure START and FINISH values are appropriate for WHICH=O\n");
    printf("Type number of parameter you wish to change (0 to continue)\n");
    GETINTEGER(&RESPONSE, 0L, 14L);
    scanf("%*[^\n]");
    getchar();
    if ((unsigned long)RESPONSE < 32 && ((1L << RESPONSE) & 0x7ffe) != 0) {
      switch (RESPONSE) {

      case 1:
	GETNUMBER(&START, "START     ", 1L, SEQLEN);
	STARTNO = START;
	break;

      case 2:
	GETNUMBER(&FINISH, "FINISH    ", 1L, SEQLEN);
	break;

      case 3:
	P_addset(P_expset(SET, 0L), 'U');
	GETCHAR(&NUCCASE, "NUCCASE   ", P_addset(SET, 'L'));
	break;

      case 4:
	GETNUMBER(&STARTNO, "STARTNO   ", -MAXSEQ, MAXSEQ);
	break;

      case 5:
	GETNUMBER(&GROUP, "GROUP     ", 3L, 160L);
	break;

      case 6:
	GETNUMBER(&GPL, "GPL       ", 1L, 160L / GROUP);
	break;

      case 7:
	P_addset(P_expset(SET1, 0L), 'I');
	GETCHAR(&WHICH, "WHICH     ", P_addset(SET1, 'O'));
	break;

      case 8:
	GETNUMBER(&STRANDS, "STRANDS   ", 1L, 2L);
	break;

      case 9:
	P_addset(P_expset(SET2, 0L), 'R');
	GETCHAR(&KIND, "KIND      ", P_addset(SET2, 'D'));
	break;

      case 10:
	P_addset(P_expset(SET3, 0L), 'Y');
	GETCHAR(&NUMBERS, "NUMBERS   ", P_addset(SET3, 'N'));
	break;

      case 11:
	P_addset(P_expset(SET3, 0L), 'Y');
	GETCHAR(&NUCS, "NUCS      ", P_addset(SET3, 'N'));
	break;

      case 12:
	P_addset(P_expset(SET3, 0L), 'Y');
	GETCHAR(&PEPTIDES, "PEPTIDES  ", P_addset(SET3, 'N'));
	break;

      case 13:
	GETNUMBER(&FRAMES, "FRAMES    ", 1L, 3L);
	break;

      case 14:
	P_addset(P_expset(SET4, 0L), 'S');
	GETCHAR(&FORM, "FORM      ", P_addset(SET4, 'L'));
	break;
      }
    }
    if (PEPTIDES == 'Y' && GROUP % 3 > 0) {
      printf(">>>> GROUP must be divisible by 3! Press <RETURN>\n");
      scanf("%*[^\n]");
      getchar();
      RESPONSE = -1;
    }
/* p2c: numseq.p, line 918:
 * Note: Using % for possibly-negative arguments [317] */
  } while (RESPONSE != 0);
  if (STARTNO == START)
    COORD = 'S';
  else
    COORD = 'U';
}  /* PARAMETERS */


/* Local variables for PRINTSEQ: */
struct LOC_PRINTSEQ {
  _TEXT *F;
  long GROUPWIDTH, DELTA, SKIP, NUCSKIP, FRAMESHIFT, THISLINE, NUMBER;
} ;

/* Skip the indicated number of spaces */
Local Void ADVANCE(NUM, LINK)
long NUM;
struct LOC_PRINTSEQ *LINK;
{
  long I;

  for (I = 1; I <= NUM; I++)
    putc(' ', LINK->F->f);
}  /* ADVANCE */

/* Write a line of numbers */
Local Void WRITENUM(LINK)
struct LOC_PRINTSEQ *LINK;
{
  long WIDTH = 0;
  long NEXTNUM;

  while (WIDTH < LINK->THISLINE) {
    /* 0 position not permitted in output. */
    NEXTNUM = LINK->NUMBER + LINK->DELTA;
    switch (COORD) {

    case 'S':
      if (WHICH == 'I' && NEXTNUM > SEQLEN)
	LINK->NUMBER += LINK->DELTA - SEQLEN;
      else if (WHICH == 'O' && NEXTNUM < 1)
	LINK->NUMBER = SEQLEN + NEXTNUM;
      else
	LINK->NUMBER = NEXTNUM;
      break;

    case 'U':
      if (NEXTNUM == 0)
	LINK->NUMBER = SENSE;
      else if ((double)LINK->NUMBER / NEXTNUM < 0)
	LINK->NUMBER = NEXTNUM + SENSE;
      else
	LINK->NUMBER = NEXTNUM;
      break;
    }
    fprintf(LINK->F->f, "%*ld", (int)LINK->GROUPWIDTH, LINK->NUMBER);
    ADVANCE(LINK->SKIP, LINK);
    WIDTH++;
  }
  putc('\n', LINK->F->f);
}  /* WRITENUM */

/* Compute the next nucleotide position */
Local long NEXT(POS, LINK)
long POS;
struct LOC_PRINTSEQ *LINK;
{
  POS += SENSE;
  if (CIRCULAR) {
    if (POS > SEQLEN)
      return 1;
    else if (POS < 1)
      return SEQLEN;
    else
      return POS;
  } else if (POS > SEQLEN || POS < 1)
    return 0;
  else
    return POS;
}  /* NEXT */

/* Write a line of bases */
Local Void WRITENUC(NUC, POS, LINK)
NUCLEOTIDE *NUC;
long *POS;
struct LOC_PRINTSEQ *LINK;
{
  long WIDTH = 0, I = 0;

  while (WIDTH < LINK->THISLINE) {
    *POS = NEXT(*POS, LINK);
    if (*POS > 0)
      putc(NUCHAR[(long)NUC[(long)SEQ[*POS - 1]]], LINK->F->f);
    WIDTH++;
    I++;
    if (I == LINK->NUCSKIP) {
      ADVANCE(LINK->SKIP, LINK);
      I = 0;
    }
  }
  putc('\n', LINK->F->f);
}  /* WRITENUC */

/* Given a codon, return an amino acid */
Local AMINOACID AAMATCH(N1, N2, N3)
NUCLEOTIDE N1, N2, N3;
{
  long I = 1;   /* rule currently pointed to */
  boolean NOTFOUND = true;

  while (NOTFOUND) {
    if (((1L << ((long)N1)) & NUCSET[(long)RULELIST[I-1][0]]) == 0) {
      I++;
      continue;
    }
    if (((1L << ((long)N2)) & NUCSET[(long)RULELIST[I-1][1]]) != 0) {
      if (((1L << ((long)N3)) & NUCSET[(long)RULELIST[I-1][2]]) != 0)
	NOTFOUND = false;
      else
	I++;
    } else
      I++;
  }
  return (AALIST[I-1]);
}  /* AAMATCH */

/* Write a line of amino acids */
Local Void WRITEAA(NUC, POS, NUMCODONS, LINK)
NUCLEOTIDE *NUC;
long POS, NUMCODONS;
struct LOC_PRINTSEQ *LINK;
{
  long WIDTH = 0;
  long P2, P3;

  POS = NEXT(POS, LINK) + LINK->FRAMESHIFT * SENSE;
  ADVANCE(LINK->FRAMESHIFT, LINK);
  while (WIDTH < NUMCODONS) {
    P2 = NEXT(POS, LINK);
    P3 = NEXT(P2, LINK);
    if (P2 > 0 && P3 > 0)
      fprintf(LINK->F->f, "%.3s",
	      AAS[(long)AAMATCH(NUC[(long)SEQ[POS-1]], NUC[(long)SEQ[P2-1]],
				NUC[(long)SEQ[P3-1]])]);
    ADVANCE(LINK->SKIP, LINK);
    WIDTH++;
    POS = NEXT(P3, LINK);
  }
  putc('\n', LINK->F->f);
}  /* WRITEAA */


/*************************************************************/
/*  Write the sequence to OUTFILE in the specified format.   */
/*************************************************************/
Static Void PRINTSEQ(F_)
_TEXT *F_;
{
  struct LOC_PRINTSEQ V;
  long I, POS, POSITION, GROUPSLEFT, NUCPERLINE, NUCLEFT, NUMCODONS, S2, S3;
  boolean LASTLINE = false;
  LINE NAME;
  long FORLIM;

  V.F = F_;
  printf("Type title to appear on output (<RETURN> for blank):\n");
  INPLINE(&NAME);
  WRITELINE(V.F, NAME, NAME.LEN);
  if (NAME.LEN > 0)
    putc('\n', V.F->f);
  /* Set parameters dependant on the strand */
  switch (WHICH) {

  case 'I':
    SENSE = 1;
    V.DELTA = GROUP;
    if (START <= FINISH)
      NUCLEFT = FINISH - START + 1;
    else if (CIRCULAR)
      NUCLEFT = SEQLEN - START + FINISH + 1;
    else
      NUCLEFT = SEQLEN - START + 1;
    break;

  case 'O':
    SENSE = -1;
    if (COORD == 'S')
      V.DELTA = -GROUP;
    else
      V.DELTA = GROUP;
    if (START >= FINISH)
      NUCLEFT = START - FINISH + 1;
    else if (CIRCULAR)
      NUCLEFT = START + SEQLEN - FINISH + 1;
    else
      NUCLEFT = START;
    break;
  }
  NUCPERLINE = GROUP * GPL;
  /* Set parameters dependant upon amino acid seq. */
  if (PEPTIDES == 'Y') {
    if (FRAMES == 3) {
      V.SKIP = 0;
      V.NUCSKIP = NUCPERLINE + 1;
    } else {
      V.SKIP = 1;
      V.NUCSKIP = 3;
    }
    V.GROUPWIDTH = GROUP + V.SKIP * (GROUP / 3) - V.SKIP;
    S2 = SENSE * 2;
    S3 = SENSE * 3;
  } else {
    V.SKIP = 1;
    V.NUCSKIP = GROUP;
    V.GROUPWIDTH = GROUP;
  }

  POSITION = START - SENSE;
  if (COORD == 'S')
    V.NUMBER = START - SENSE;
  else
    V.NUMBER = STARTNO - 1;
  GROUPSLEFT = NUCLEFT / GROUP;

  do {
    /* Write numbers */
    if (GROUPSLEFT < GPL)
      V.THISLINE = GROUPSLEFT;
    else
      V.THISLINE = GPL;
    if (NUMBERS == 'Y')
      WRITENUM(&V);
    GROUPSLEFT -= V.THISLINE;

    /* Write bases */
    if (NUCLEFT < NUCPERLINE)
      V.THISLINE = NUCLEFT;
    else
      V.THISLINE = NUCPERLINE;
    POS = POSITION;
    if (NUCS == 'Y') {
      switch (WHICH) {

      case 'I':
	WRITENUC(INUC, &POS, &V);
	break;

      case 'O':
	WRITENUC(COMP, &POS, &V);
	break;
      }
      if (STRANDS == 2) {
	POS = POSITION;
	switch (WHICH) {

	case 'I':
	  WRITENUC(COMP, &POS, &V);
	  break;

	case 'O':
	  WRITENUC(INUC, &POS, &V);
	  break;
	}
      }
    } else {
      FORLIM = V.THISLINE;
      for (I = 1; I <= FORLIM; I++)
	POS = NEXT(POS, &V);
    }
    NUCLEFT -= V.THISLINE;

    /* Write amino acids */
    if (PEPTIDES == 'Y') {
      V.FRAMESHIFT = 0;
      while (V.FRAMESHIFT < FRAMES) {
	NUMCODONS = V.THISLINE / 3;
	switch (WHICH) {

	case 'I':
	  WRITEAA(INUC, POSITION, NUMCODONS, &V);
	  break;

	case 'O':
	  WRITEAA(COMP, POSITION, NUMCODONS, &V);
	  break;
	}
	V.FRAMESHIFT++;
      }
    }
    POSITION = POS;

    /* Write a blank line */
    if (NUMBERS == 'Y')
      putc('\n', V.F->f);
    if (NUCLEFT == 0)
      LASTLINE = true;
  } while (!LASTLINE);   /* PRINTSEQ */
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
  GCFILE.f = NULL;
  *GCFILE.name = '\0';
  OUTFILE.f = NULL;
  *OUTFILE.name = '\0';
  INFILE.f = NULL;
  *INFILE.name = '\0';
  printf("%50s\n\n", VERSION);
  /* END MODULE STARTUP         VERSION= 'SUNMODS     Version  6/26/01'; */

  /* Set up initial values for nucleotide tables and genetic code */
  INITTABLES();
  MAKERULES(GC);

  /* Initialize horizontal output line */

  /* Read in initial sequence and set up parameters. */
  for (I = 1; I <= MAXLINE; I++)
    HLINE.STR[I-1] = '_';
  HLINE.LEN = MAXLINE;
  printf("Enter sequence filename:\n");
  GETFILE(&INFILE, 'I', &IFN);
  READSEQ(&INFILE, SEQ, &SEQLEN, &SEQNAME, &CIRCULAR);
  INITPARAM();

  /* Open initial output file */
  printf("Type output filename:\n");
  GETFILE(&OUTFILE, 'O', &OFN);

  /* MAIN MENU */
  GCFN.LEN = 0;
  do {
    putchar('\n');
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, HLINE, 80L);
    printf("\nNUMSEQ%30s\n", "MAIN MENU");
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
    printf("\nGenetic code file: ");
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, GCFN, GCFN.LEN);
    putchar('\n');
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, HLINE, 80L);
    printf("\n%20c1) Read in a new sequence\n", ' ');
    printf("%20c2) Open a new output file\n", ' ');
    printf("%20c3) Read in an alternative genetic code\n", ' ');
    printf("%20c4) Change parameters\n", ' ');
    printf("%20c5) Write output to screen\n", ' ');
    printf("%20c6) Write output to file\n", ' ');
    TEMP.f = stdout;
    *TEMP.name = '\0';
    WRITELINE(&TEMP, HLINE, 80L);
    printf("\nType the number of your choice  (0 to quit program)\n");
    GETINTEGER(&CHOICE, 0L, 6L);
    scanf("%*[^\n]");
    getchar();

    switch (CHOICE) {

    case 0:
      /* blank case */
      break;

    case 1:
      printf("Enter sequence filename:\n");
      GETFILE(&INFILE, 'I', &IFN);
      READSEQ(&INFILE, SEQ, &SEQLEN, &SEQNAME, &CIRCULAR);
      START = 1;
      FINISH = SEQLEN;
      STARTNO = START;
      break;

    case 2:
      if (OUTFILE.f != NULL)
	fclose(OUTFILE.f);
      OUTFILE.f = NULL;
      printf("Type output filename:\n");
      GETFILE(&OUTFILE, 'O', &OFN);
      break;

    /*!!!*/
    case 3:
      READCODE(&GCFILE, &GCFN, GC);
      MAKERULES(GC);
      break;

    case 4:
      PARAMETERS();
      break;

    case 5:
      INITSTRING();
      TEMP.f = stdout;
      *TEMP.name = '\0';
      PRINTSEQ(&TEMP);
      printf("Press RETURN to continue");
      scanf("%*[^\n]");
      getchar();
      break;

    case 6:
      INITSTRING();
      PRINTSEQ(&OUTFILE);
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
  if (GCFILE.f != NULL)
    fclose(GCFILE.f);
  exit(EXIT_SUCCESS);
}  /* NUMSEQ */



/* End. */
