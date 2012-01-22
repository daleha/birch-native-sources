/* Output from p2c 1.21alpha-07.Dec.93, the Pascal-to-C translator */
/* From input file "ribosome.p" */


/***********************************************************/
/*                                                         */
/*  RIBOSOME  VERSION   8/12/93  Standard Pascal           */
/*            Brian Fristensky                             */
/*            Dept. of Plant Science                       */
/*            University of Manitoba                       */
/*            Winnipeg, MB R3T 2N2   CANADA                */
/*                                                         */
/*  Copyright (c) l989, 1990 by Brian Fristensky           */
/*  !!! in comment indicates feature which may need change */
/***********************************************************/

#include <p2c.h>


/*,GCFILE*/
/*!!!  Some Pascals require file parameters in program heading */

#define MAXSEQ          32700
#define MAXWORD         25
#define STARTARGNUM     1


typedef enum {
  T, C, A, G, R, Y, M, W, S, K, D, H, V, B, N
} NUCLEOTIDE;
typedef long SS;

typedef enum {
  ALA, CYS, ASP, GLU, PHE, GLY, HIS, ILE, LYS, LEU, MET, ASN, PRO, GLN, ARG,
  SER, THR, VALINE, TRP, TYR, ASX, GLX, TERM, UNKX
} AMINOACID;
/* Note: VAL is a reserved word in some Pascals */

typedef NUCLEOTIDE SEQUENCE[MAXSEQ];

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

/* END MODULE TYPE.WORD         VERSION= 'SUNMODS     Version  9/12/91'; */


/* File variables */
Static _TEXT GCFILE;   /* genetic code file */
Static long ARGNUM;
Static Char ARGUMENT[132];

/* Sequence variables */
Static SEQUENCE SEQ;
Static long SEQLEN;
Static WORD SEQNAME;
Static SS NUCSET[15];
Static Char NUCHAR[15];

/* Genetic code variables */
Static GENETICCODE GC;
Static RULE RULELIST[64];
Static AMINOACID AALIST[64];
Static long NUMRULES;
Static Char AAS[24];

/* GLOBAL PARAMETERS */
Static long START, FINISH;


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


/* END MODULE READWORD         VERSION= 'SUNMODS     Version  9/12/91'; */

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


/* END MODULE WRITEWORD         VERSION= 'SUNMODS     Version  9/12/91'; */

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


/* END MODULE NUC         VERSION= 'SUNMODS     Version  9/12/91'; */

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


/* END MODULE AA         VERSION= 'SUNMODS     Version  9/12/91'; */

/*******************************************************/
/* Initialize tables of nucleotides, amino acids etc.  */
/* This is done only once when the program begins.     */
/*******************************************************/
Static Void INITTABLES()
{
  /*   Initialize NUCHAR the array which */
  /*   hold the  values of the NUCLEOTIDEs.      */
  NUCHAR[(long)T] = 'T';
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


/* ************************************************************ */
/* Initialize amino acid output strings according to parameters.*/
/* ************************************************************ */
Static Void INITSTRING()
{
  AAS[(long)PHE] = 'F';
  AAS[(long)LEU] = 'L';
  AAS[(long)ILE] = 'I';
  AAS[(long)MET] = 'M';
  AAS[(long)VALINE] = 'V';
  AAS[(long)SER] = 'S';
  AAS[(long)PRO] = 'P';
  AAS[(long)THR] = 'T';
  AAS[(long)ALA] = 'A';
  AAS[(long)TYR] = 'Y';
  AAS[(long)HIS] = 'H';
  AAS[(long)GLN] = 'Q';
  AAS[(long)ASN] = 'N';
  AAS[(long)LYS] = 'K';
  AAS[(long)ASP] = 'D';
  AAS[(long)GLU] = 'E';
  AAS[(long)CYS] = 'C';
  AAS[(long)TRP] = 'W';
  AAS[(long)ARG] = 'R';
  AAS[(long)GLY] = 'G';
  AAS[(long)ASX] = 'B';
  AAS[(long)GLX] = 'Z';
  AAS[(long)TERM] = '*';
  AAS[(long)UNKX] = 'X';
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
Static Void READCODE(F, GC)
_TEXT *F;
AMINOACID (*GC)[(long)G - (long)T + 1][(long)G - (long)T + 1];
{
  NUCLEOTIDE N1, N2, N3;
  Char AACHAR;

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
/* p2c: ribosome.p, line 227: Note:
 * Line breaker spent 1.2+0.23 seconds, 2661 tries on line 563 [251] */
	  AACHAR = getc(F->f);
	  if (AACHAR == '\n')
	    AACHAR = ' ';
	}
	GC[(long)N1 - (long)T][(long)N2 - (long)T]
	  [(long)N3 - (long)T] = AA(AACHAR);
      }
    }
  }
  /*!!! CLOSE(F) */
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
/* p2c: ribosome.p, line 378: Note:
 * Line breaker spent 1.1+0.73 seconds, 2060 tries on line 724 [251] */
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
  /*    INITSTRING;
        writeln('The following rules implement the genetic code:');
        for RA:= 1 to NUMRULES do begin
          for RB:= 1 to 3 do write(NUCHAR[RULELIST[RA][RB]]);
          write(' ':7);
          writeln(AAS[AALIST[RA]]:10)
          end
        writeln;
        readln */
}  /* MAKERULES */


/*************************************************************/
/*  Read a nucleic acid sequence from SEQFILE.               */
/*************************************************************/
Static Void READSEQUENCE(SEQFILE, NAME, SEQ, START, FINISH)
_TEXT *SEQFILE;
WORD *NAME;
NUCLEOTIDE *SEQ;
long *START, *FINISH;
{
  Char CH;
  boolean TOOLONG = false;

  SEQLEN = 0;
  /* Read in the sequence name, if present. Name may be preceeded by
     comment lines. */
  while (P_peek(SEQFILE->f) == ';') {
    fscanf(SEQFILE->f, "%*[^\n]");
    getc(SEQFILE->f);
  }
  if (P_peek(SEQFILE->f) == '>') {
    READWORD(SEQFILE, NAME);
    fscanf(SEQFILE->f, "%*[^\n]");
    getc(SEQFILE->f);
  }

  /* Read in the sequence */
  while ((!BUFEOF(SEQFILE->f)) & (P_peek(SEQFILE->f) != '>')) {
    if (P_peek(SEQFILE->f) != ';') {   /* comment line, ignore */
      while (!P_eoln(SEQFILE->f)) {
	CH = getc(SEQFILE->f);
	if (CH == '\n')
	  CH = ' ';
	if (CH == 'b' || CH == 'B' || CH == 'v' || CH == 'V' || CH == 'h' ||
	    CH == 'H' || CH == 'd' || CH == 'D' || CH == 'k' || CH == 'K' ||
	    CH == 's' || CH == 'S' || CH == 'w' || CH == 'W' || CH == 'm' ||
	    CH == 'M' || CH == 'y' || CH == 'Y' || CH == 'r' || CH == 'R' ||
	    CH == 'n' || CH == 'N' || CH == 'u' || CH == 'U' || CH == 't' ||
	    CH == 'T' || CH == 'g' || CH == 'G' || CH == 'c' || CH == 'C' ||
	    CH == 'a' || CH == 'A') {
/* p2c: ribosome.p, line 415: Note:
 * Line breaker spent 1.8+0.38 seconds, 5000 tries on line 795 [251] */
	  if (SEQLEN < MAXSEQ - 2) {
	    SEQLEN++;
	    SEQ[SEQLEN-1] = NUC(CH);
	  } else
	    TOOLONG = true;
	}
      }  /* eoln */
    }
    fscanf(SEQFILE->f, "%*[^\n]");
    getc(SEQFILE->f);
  }  /* eof */
  if (TOOLONG) {
    printf(">>> WARNING! Sequence length exceeds MAXSEQ-2. Seq. truncated.\n");
  }  /* TOOLONG */
  *START = 1;
  *FINISH = SEQLEN;
}  /* READSEQUENCE */


#define CODONSPERLINE   50


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


/*************************************************************/
/*  Translate sequence and write to output.                  */
/*************************************************************/
Static Void TRANSLATE(F, NAME, SEQ, START, FINISH)
_TEXT *F;
WORD *NAME;
NUCLEOTIDE *SEQ;
long *START, *FINISH;
{
  long P1, P2, P3, CODONSLEFT, EXTRANUCS, THISLINE, TOTALDISTANCE,
       TOTALCODONS, WIDTH;


  if (NAME->LEN > 0) {
    WRITEWORD(F, *NAME, NAME->LEN);
    putc('\n', F->f);
  }
  TOTALDISTANCE = *FINISH - *START + 1;
  TOTALCODONS = TOTALDISTANCE / 3;

  /* If sequence is truncated within last codon, add one or two N's to
     complete the codon. Sequence ambiguity may allow another aa to
     be assigned. */
  EXTRANUCS = TOTALDISTANCE % 3;
/* p2c: ribosome.p, line 453:
 * Note: Using % for possibly-negative arguments [317] */
  if (EXTRANUCS > 0) {
    SEQ[SEQLEN] = N;
    if (EXTRANUCS == 2)
      SEQ[SEQLEN+1] = N;
    TOTALCODONS++;
  }  /* EXTRANUCS > 0 */

  /* Translate the sequence */
  P1 = *START;
  CODONSLEFT = TOTALCODONS;
  while (CODONSLEFT > 0) {
    if (CODONSLEFT >= CODONSPERLINE)
      THISLINE = CODONSPERLINE;
    else
      THISLINE = CODONSLEFT;
    for (WIDTH = 1; WIDTH <= THISLINE; WIDTH++) {
      P2 = P1 + 1;
      P3 = P2 + 1;
      fputc(AAS[(long)AAMATCH(SEQ[P1-1], SEQ[P2-1], SEQ[P3-1])], F->f);
      P1 = P3 + 1;
    }  /* for WIDTH */
    putc('\n', F->f);
    CODONSLEFT -= THISLINE;
  }  /* CODONSLEFT>0 */
}  /* TRANSLATE */

#undef CODONSPERLINE


/* ----------------------------------------------------------  */
/* ----------------- MAIN  PROCEDURE  -----------------------  */
main(argc, argv)
int argc;
Char *argv[];
{
  _TEXT TEMP;

  PASCAL_MAIN(argc, argv);
  GCFILE.f = NULL;
  *GCFILE.name = '\0';
  /* Set up initial values for nucleotide tables and genetic code */
  INITTABLES();
  INITSTRING();

  /* Read in a genetic code file, if specified by command line option -g*/
  if (P_argc > 1) {
    ARGNUM = 2;
    P_sun_argv(ARGUMENT, 132, (int)ARGNUM);
    if (ARGUMENT[0] == '-' && ARGUMENT[1] == 'g') {
      ARGNUM = 3;
      FILEARGS(&GCFILE, 'I', &ARGNUM);
      READCODE(&GCFILE, GC);
    }
  }

  /* Derive ambiguity rules based on genetic code */
  MAKERULES(GC);

  /* Read each sequence and translate while writing to output.*/
  while (!P_eof(stdin)) {
    TEMP.f = stdin;
    *TEMP.name = '\0';
    READSEQUENCE(&TEMP, &SEQNAME, SEQ, &START, &FINISH);
    TEMP.f = stdout;
    *TEMP.name = '\0';
    TRANSLATE(&TEMP, &SEQNAME, SEQ, &START, &FINISH);
    /* Advance to first line beginning with '>' */
    while ((P_peek(stdin) != '>') & (!P_eof(stdin))) {
      scanf("%*[^\n]");
      getchar();
    }
  }
  if (GCFILE.f != NULL)
    fclose(GCFILE.f);
  exit(EXIT_SUCCESS);
}  /* RIBOSOME */



/* End. */
