/* Output from p2c 1.21alpha-07.Dec.93, the Pascal-to-C translator */
/* From input file "prot2nuc.p" */


/* ********************************************************  */
/*                                                           */
/*   PROT2NUC  Version  8/10/94, Standard Pascal             */
/*             Brian Fristensky                              */
/*             Dept. of Plant Science                        */
/*             University of Manitoba                        */
/*             Winnipeg, MB R3T 2N2 CANADA                   */
/*                                                           */
/* SYNOPSIS                                                  */
/*    prot2nuc [-l<n> -g<n>]                                 */
/*                                                           */
/* DESCRIPTION                                               */
/*    Reverse translates protein to nucleic acid sequence.   */
/*    File is Pearson (.wrp) format file. see Lipman et al   */
/*                                                           */
/*    -l<n>  integer, length of output line in codons        */
/*           (default = 20)                                  */
/*                                                           */
/*    -g<n>  integer; number ever n codons (default = 5)     */
/*                                                           */
/* Copyright (c) 1994 by Brian Fristensky.                   */
/* !!! in comment indicates feature which may need change.   */
/*  *******************************************************  */


#include <p2c.h>


/*!!! Some Pascals require file parameters in program heading */

#define VERSION         "PROT2NUC       Version  8/10/94"

#define MAXSEQ          9000
#define MAXWORD         25
#define MAXLINE         70   /* only used for ARGUMENT */
/* BEGIN MODULE STARTARGNUM */
#define STARTARGNUM     1
    /* SUN Pascal: ARG(1) is 1st command line argument*/


/*STARTARGNUM=2;       ATT Pascal: ARG(2) is 1st command line argument*/
/* END MODULE STARTARGNUM         VERSION= 'SUNMODS     Version  8/ 9/94'; */


typedef enum {
  GLY, ALA, VALINE, LEU, ILE, MET, PHE, PRO, SER, THR, CYS, ASN, GLN, TYR,
  TRP, ASP, GLU, HIS, LYS, ARG, ASX, GLX, TERM, UNKX
} AMINOACID;
typedef AMINOACID PROTEIN[MAXSEQ];
typedef Char TRIPLET[3];
typedef TRIPLET DEGENERATE[24];

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


/* Command line argument variables. */
Static CHARARRAY ARGUMENT;   /*     command line argument */
Static long ARGNUM;   /* # of     "     "      "   */
Static long POS;   /* posn. in "     "      "   */
Static boolean UNREADOPTIONS;

/* Amino acid sequence variables */
Static WORD NAME;
Static PROTEIN PR;
Static long SEQLEN;

/* Variables for printing */
Static long CODONSPERLINE, GROUP;
Static Char AACHAR[24];
Static DEGENERATE CODON1, CODON2;


#define C2STR           "   "


/* **************************************************** */
/*              INITIALIZATION  PROCEDURES              */
/* **************************************************** */
Static Void INITIALIZE()
{
  AMINOACID AA1;

  /*   Initialize AACHAR array which holds the            */
  /*   character values of the AMINOACIDs.                */
  AACHAR[(long)GLY] = 'G';
  AACHAR[(long)ALA] = 'A';
  AACHAR[(long)VALINE] = 'V';
  AACHAR[(long)LEU] = 'L';
  AACHAR[(long)ILE] = 'I';
  AACHAR[(long)MET] = 'M';
  AACHAR[(long)PHE] = 'F';
  AACHAR[(long)PRO] = 'P';
  AACHAR[(long)SER] = 'S';
  AACHAR[(long)THR] = 'T';
  AACHAR[(long)CYS] = 'C';
  AACHAR[(long)ASN] = 'N';
  AACHAR[(long)GLN] = 'Q';
  AACHAR[(long)TYR] = 'Y';
  AACHAR[(long)TRP] = 'W';
  AACHAR[(long)ASP] = 'D';
  AACHAR[(long)GLU] = 'E';
  AACHAR[(long)HIS] = 'H';
  AACHAR[(long)LYS] = 'K';
  AACHAR[(long)ARG] = 'R';
  AACHAR[(long)ASX] = 'B';
  AACHAR[(long)GLX] = 'Z';
  AACHAR[(long)TERM] = '*';
  AACHAR[(long)UNKX] = 'X';

  /* Initialize CODON1 & CODON2. All aminoacids have an equivalent
     TRIPLET in the CODON1 array. For amino acids that require two
     degenerate codons (eg. Leu = TTR, CTN), a triplet is needed
     for CODON2. For amino acids that can be represented in a single
     degenerate codon (eg. Cys = TGY), the TRIPLET in CODON2 is
     the value of the constant C2STR. */
  memcpy(CODON1[(long)PHE], "TTy", sizeof(TRIPLET));
  memcpy(CODON1[(long)LEU], "CTn", sizeof(TRIPLET));
  memcpy(CODON1[(long)SER], "TCn", sizeof(TRIPLET));
  memcpy(CODON1[(long)TYR], "TAy", sizeof(TRIPLET));
  memcpy(CODON1[(long)TERM], "TAr", sizeof(TRIPLET));
  memcpy(CODON1[(long)CYS], "TGy", sizeof(TRIPLET));
  memcpy(CODON1[(long)TRP], "TGG", sizeof(TRIPLET));
  memcpy(CODON1[(long)PRO], "CCn", sizeof(TRIPLET));
  memcpy(CODON1[(long)HIS], "CAy", sizeof(TRIPLET));
  memcpy(CODON1[(long)GLN], "CAr", sizeof(TRIPLET));
  memcpy(CODON1[(long)ARG], "CGn", sizeof(TRIPLET));
  memcpy(CODON1[(long)ILE], "ATh", sizeof(TRIPLET));
  memcpy(CODON1[(long)MET], "ATG", sizeof(TRIPLET));
  memcpy(CODON1[(long)THR], "ACn", sizeof(TRIPLET));
  memcpy(CODON1[(long)ASN], "AAy", sizeof(TRIPLET));
  memcpy(CODON1[(long)LYS], "AAr", sizeof(TRIPLET));
  memcpy(CODON1[(long)VALINE], "GTn", sizeof(TRIPLET));
  memcpy(CODON1[(long)ALA], "GCn", sizeof(TRIPLET));
  memcpy(CODON1[(long)ASP], "GAy", sizeof(TRIPLET));
  memcpy(CODON1[(long)GLU], "GAr", sizeof(TRIPLET));
  memcpy(CODON1[(long)GLY], "GGn", sizeof(TRIPLET));
  memcpy(CODON1[(long)ASX], "rAy", sizeof(TRIPLET));
  memcpy(CODON1[(long)GLX], "sAr", sizeof(TRIPLET));
  memcpy(CODON1[(long)UNKX], "nnn", sizeof(TRIPLET));

  for (AA1 = GLY; (long)AA1 <= (long)UNKX; AA1 = (AMINOACID)((long)AA1 + 1))
    memcpy(CODON2[(long)AA1], C2STR, sizeof(TRIPLET));
  memcpy(CODON2[(long)LEU], "TTr", sizeof(TRIPLET));
  memcpy(CODON2[(long)TERM], "TGA", sizeof(TRIPLET));
  memcpy(CODON2[(long)ARG], "AGr", sizeof(TRIPLET));
  memcpy(CODON2[(long)SER], "AGy", sizeof(TRIPLET));
}  /* INITIALIZE */

#undef C2STR


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


/* END MODULE AA         VERSION= 'SUNMODS     Version  8/ 9/94'; */

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

/* BEGIN MODULE NUMBER */
/* Extract an integer from a CHARARRAY (see TYPE.LINE), starting
   at the current position. */
Static long NUMBER(TARGET, POS)
Char *TARGET;
long *POS;
{
  long N = 0, ORDZERO = '0';

  /* evaluate characteristic */
  while (isdigit(TARGET[*POS - 1])) {
    N = N * 10 + TARGET[*POS - 1] - ORDZERO;
    (*POS)++;
  }
  return N;
}  /* NUMBER */


/* END MODULE NUMBER         VERSION= 'SUNMODS     Version  8/ 9/94'; */

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

/*  ******************************************* */
/*  Read a protein    sequence from SEQFILE     */
/*   and store it in PR.                        */
/*  ******************************************* */
Static Void READPRO(SEQFILE, PR, SEQLEN, NAME)
_TEXT *SEQFILE;
AMINOACID *PR;
long *SEQLEN;
WORD *NAME;
{
  Char CH;
  boolean TOOLONG = false;

  *SEQLEN = 0;
  /* Read in the sequence name, if present. Name may be preceeded by
     comment lines. */
  while (P_peek(SEQFILE->f) == ';') {
    fscanf(SEQFILE->f, "%*[^\n]");
    getc(SEQFILE->f);
  }
  if (P_peek(SEQFILE->f) == '>') {
    CH = getc(SEQFILE->f);
    if (CH == '\n')
      CH = ' ';
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
	if (CH == 'X' || CH == '*' || CH == 'Z' || CH == 'B' || CH == 'R' ||
	    CH == 'K' || CH == 'H' || CH == 'E' || CH == 'D' || CH == 'W' ||
	    CH == 'Y' || CH == 'Q' || CH == 'N' || CH == 'C' || CH == 'T' ||
	    CH == 'S' || CH == 'P' || CH == 'F' || CH == 'M' || CH == 'I' ||
	    CH == 'L' || CH == 'V' || CH == 'A' || CH == 'G') {
	  if (*SEQLEN < MAXSEQ - 2) {
	    (*SEQLEN)++;
	    PR[*SEQLEN - 1] = AA(CH);
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
}  /* READPRO */


/*  ******************************************* */
/*  Print a header listing the amino acid       */
/*  and nucleic acid one-letter codes.          */
/*  ******************************************* */
Static Void HEADER()
{
  Char STR1[256];

  printf("%s\n\n", VERSION);
  printf("     IUPAC-IUP AMINO ACID SYMBOLS\n");
  printf("     [J. Biol. Chem. 243, 3557-3559 (1968)]\n\n");
  printf("          Phe         F          Leu         L          Ile         I\n");
  printf("          Met         M          Val         V          Ser         S\n");
  printf("          Pro         P          Thr         T          Ala         A\n");
  printf("          Tyr         Y          His         H          Gln         Q\n");
  printf("          Asn         N          Lys         K          Asp         D\n");
  printf("          Glu         E          Cys         C          Trp         W\n");
  printf("          Arg         R          Gly         G          STOP        *\n");
  printf(
    "          Asx         B          Glx         Z          UNKNOWN     X\n\n\n");
  printf("     IUPAC-IUB SYMBOLS FOR NUCLEOTIDE NOMENCLATURE\n");
  printf("     [Cornish-Bowden (1985) Nucl. Acids Res. 13: 3021-3030.]\n\n");
  printf("     Symbol         Meaning              | Symbol         Meaning\n");
  printf(
    "     ------------------------------------+---------------------------------\n");
  printf("     G              Guanine              | k              G or T\n");
  printf("     A              Adenine              | s              G or C\n");
  printf("     C              Cytosine             | w              A or T\n");
  printf("     T              Thymine              | h              A or C or T\n");
  printf("     U              Uracil               | b              G or T or C\n");
  printf("     r              Purine (A or G)      | v              G or C or A\n");
  printf("     y              Pyrimidine (C or T)  | d              G or T or A\n");
  printf(
    "     m              A or C               | n              G or A or T or C\n\n");
}  /* HEADER */


/* Local variables for PRINTSEQ: */
struct LOC_PRINTSEQ {
  _TEXT *F;
  PROTEIN PR;
  long SEQLEN, GROUPWIDTH;
} ;

/* Write a line of numbers */
Local Void WRITENUM(THISLINE, NUMBER, LINK)
long THISLINE, *NUMBER;
struct LOC_PRINTSEQ *LINK;
{
  long WIDTH = 0;
  long NEXTNUM;

  while (WIDTH < THISLINE) {
    NEXTNUM = *NUMBER + GROUP;
    if (NEXTNUM <= LINK->SEQLEN) {
      *NUMBER = NEXTNUM;
      fprintf(LINK->F->f, "%*ld", (int)LINK->GROUPWIDTH, *NUMBER);
    }
    WIDTH += GROUP;
  }
  putc('\n', LINK->F->f);
}  /* WRITENUM */

/* Write a line of amino acids */
Local Void WRITEAA(THISLINE, POS, LINK)
long THISLINE, *POS;
struct LOC_PRINTSEQ *LINK;
{
  long WIDTH = 0;

  while (WIDTH < THISLINE) {
    (*POS)++;
    fprintf(LINK->F->f, "%c  ", AACHAR[(long)LINK->PR[*POS - 1]]);
    WIDTH++;
  }
  putc('\n', LINK->F->f);
}  /* WRITEAA */

/* Write a line of degenerate triplets */
Local Void WRITENUC(THISLINE, NUCSEQ, POS, LINK)
long THISLINE;
TRIPLET *NUCSEQ;
long *POS;
struct LOC_PRINTSEQ *LINK;
{
  long WIDTH = 0;

  while (WIDTH < THISLINE) {
    (*POS)++;
    fprintf(LINK->F->f, "%.3s", NUCSEQ[(long)LINK->PR[*POS - 1]]);
    WIDTH++;
  }
  putc('\n', LINK->F->f);
}  /* WRITENUC */


/*************************************************************/
/*  Write the sequence to OUTFILE in the specified format.   */
/*************************************************************/
Static Void PRINTSEQ(F_, PR_, SEQLEN_, NAME)
_TEXT *F_;
AMINOACID *PR_;
long SEQLEN_;
WORD NAME;
{
  struct LOC_PRINTSEQ V;
  long POS;
  long POSITION = 0;
  long CODONSLEFT, THISLINE;
  long NUMBER = 0;
  boolean LASTLINE = false;

  V.F = F_;
  memcpy(V.PR, PR_, sizeof(PROTEIN));
  V.SEQLEN = SEQLEN_;
  WRITEWORD(V.F, NAME, NAME.LEN);
  if (NAME.LEN > 0)
    putc('\n', V.F->f);
  CODONSLEFT = V.SEQLEN;
  V.GROUPWIDTH = GROUP * 3;   /* width in characters, rather than codons */

  do {
    if (CODONSLEFT < CODONSPERLINE)
      THISLINE = CODONSLEFT;
    else
      THISLINE = CODONSPERLINE;

    /* Write numbers */
    WRITENUM(THISLINE, &NUMBER, &V);

    /* Write amino acids */
    POS = POSITION;
    WRITEAA(THISLINE, &POS, &V);

    /* Write bases */
    POS = POSITION;
    WRITENUC(THISLINE, CODON1, &POS, &V);
    POS = POSITION;
    WRITENUC(THISLINE, CODON2, &POS, &V);
    POSITION = POS;

    /* Write a blank line */
    putc('\n', V.F->f);
    CODONSLEFT -= THISLINE;
    if (CODONSLEFT == 0)
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

  PASCAL_MAIN(argc, argv);
  INITIALIZE();
  CODONSPERLINE = 25;
  GROUP = 5;

  /* Read options from the command line. */
  ARGNUM = STARTARGNUM;
  UNREADOPTIONS = true;
  while (UNREADOPTIONS) {
    if (ARGNUM >= P_argc) {
      UNREADOPTIONS = false;
      break;
    }  /* if ARGNUM <= argc */
    P_sun_argv(ARGUMENT, sizeof(CHARARRAY), (int)ARGNUM);
    if (ARGUMENT[0] == '-') {
      if (ARGUMENT[1] == 'g' || ARGUMENT[1] == 'l') {
	POS = 3;
	switch (ARGUMENT[1]) {

	case 'l':
	  CODONSPERLINE = NUMBER(ARGUMENT, &POS);
	  break;

	case 'g':
	  GROUP = NUMBER(ARGUMENT, &POS);
	  break;
	}
      }  /* l,g */
    }  /* '-' */
    else
      UNREADOPTIONS = false;
    ARGNUM++;
  }

  /* GROUP must divide evenly into CODONSPERLINE */
  if (CODONSPERLINE % GROUP != 0) {
/* p2c: prot2nuc.p, line 373:
 * Note: Using % for possibly-negative arguments [317] */
    GROUP = 5;
    CODONSPERLINE = 25;
  }

  TEMP.f = stdin;
  *TEMP.name = '\0';
  /* Read the protein and print the reverse translation. */
  READPRO(&TEMP, PR, &SEQLEN, &NAME);
  HEADER();
  TEMP.f = stdout;
  *TEMP.name = '\0';
  PRINTSEQ(&TEMP, PR, SEQLEN, NAME);
  exit(EXIT_SUCCESS);
}  /* PROT2NUC */




/* End. */
