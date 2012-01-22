/* Output from p2c 1.21alpha-07.Dec.93, the Pascal-to-C translator */
/* From input file "prostat.p" */


/*$DEBUG- */
/* ********************************************************  */
/*                                                           */
/*    PROSTAT  Version  8/30/93  Standard Pascal             */
/*             Brian Fristensky                              */
/*             Dept. of Plant Science                        */
/*             University of Manitoba                        */
/*             Winnipeg, MB  Canada  R3T 2N2                 */
/*                                                           */
/* Copyright (c) 1984 - 1993 by Brian Fristensky.            */
/* !!! in comment indicates feature which may need change.   */
/*  *******************************************************  */


#include <p2c.h>


/*,SFILE,OUTFILE*/
/*!!! Some Pascals require file parameters in program heading */

#define MAXSEQ          10000
#define MAXRANGE        0
#define MAXLINE         150
#define MAXWORD         25

#define VERSION         "PROSTAT      Version  8/30/93"


typedef enum {
  GLY, ALA, VALINE, LEU, ILE, MET, PHE, PRO, SER, THR, CYS, ASN, GLN, TYR,
  TRP, ASP, GLU, HIS, LYS, ARG, ASX, GLX, TERM, UNKX
} AMINOACID;
typedef AMINOACID PROTEIN[MAXSEQ];
typedef long AAI[24];
typedef double AAR[24];

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


Static _TEXT SFILE, OUTFILE;
Static LINE SFN, OFN;

/* Variables associated with the test protein */
Static PROTEIN SEQ;
Static long SEQLEN;
Static WORD NAME;
Static double MW;
Static AAI AACOMP;   /* number of each amino acid in the protein */
Static AAR AAWT;   /* molecular weight of each amino acid */


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

/* **************************************************** */
/*              INITIALIZATION  PROCEDURES              */
/* **************************************************** */
Static Void INITIALIZE()
{
  AMINOACID AA1;

  for (AA1 = GLY; (long)AA1 <= (long)UNKX; AA1 = (AMINOACID)((long)AA1 + 1))
    AACOMP[(long)AA1] = 0;

  /* Molecular weights of amino acids, water subtracted out */
  AAWT[(long)GLY] = 57.0;
  AAWT[(long)ALA] = 71.0;
  AAWT[(long)VALINE] = 99.0;
  AAWT[(long)LEU] = 113.0;
  AAWT[(long)ILE] = 113.0;
  AAWT[(long)MET] = 131.0;
  AAWT[(long)PHE] = 147.0;
  AAWT[(long)PRO] = 97.0;
  AAWT[(long)SER] = 87.0;
  AAWT[(long)THR] = 101.0;
  AAWT[(long)CYS] = 103.0;
  AAWT[(long)ASN] = 114.0;
  AAWT[(long)GLN] = 128.0;
  AAWT[(long)TYR] = 163.0;
  AAWT[(long)TRP] = 186.0;
  AAWT[(long)ASP] = 115.0;
  AAWT[(long)GLU] = 129.0;
  AAWT[(long)HIS] = 137.0;
  AAWT[(long)LYS] = 128.0;
  AAWT[(long)ARG] = 156.0;
  AAWT[(long)ASX] = 114.5;
  AAWT[(long)GLX] = 128.5;
  AAWT[(long)TERM] = 0.0;
  AAWT[(long)UNKX] = 119.0;
}  /* INITIALIZE */


/* BEGIN MODULE READPRO */
/*  ******************************************* */
/*  Read a protein    sequence from SFILE       */
/*   and store it in PR.                        */
/*  ******************************************* */
Static Void READPRO(SFILE, PR, LEN, NAME)
_TEXT *SFILE;
AMINOACID *PR;
long *LEN;
WORD *NAME;
{
  Char CH, ANSWER, FILETYPE;
  long J;
  boolean TOOLONG = false;
  long FORLIM;

  /* Prompt user for sequence file type */
  printf("The following file formats can be read:\n");
  printf("  F:free format    N:NBRF\n");
  do {
    printf("Type letter of format (F|N)\n");
    scanf("%c%*[^\n]", &ANSWER);
    getchar();
    if (ANSWER == '\n')
      ANSWER = ' ';
  } while (ANSWER != 'n' && ANSWER != 'N' && ANSWER != 'f' && ANSWER != 'F');
  switch (ANSWER) {

  case 'F':
  case 'f':
    FILETYPE = 'F';
    break;

  case 'N':
  case 'n':
    FILETYPE = 'N';
    break;
  }
  printf("Reading input file...\n");

  /* NBRF: read sequence name and title line */
  /* 'Official' NBRF files have a name line followed by a title line.
      In the name line, four data characters preceed the name itself.
      These are deleted by the program.  Short NBRF files have a title
      following the name, on the same line. The name is not preceeded
      by data characters. */
  NAME->LEN = 0;
  if (FILETYPE == 'N') {
    if (P_peek(SFILE->f) == '>') {
      CH = getc(SFILE->f);
      if (CH == '\n')
	CH = ' ';
    }
    READWORD(SFILE, NAME);
    if (NAME->STR[2] == ';') {
      FORLIM = NAME->LEN;
      /* with NAME */
      for (J = 4; J <= FORLIM; J++)
	NAME->STR[J-4] = NAME->STR[J-1];
      NAME->LEN -= 3;
      fscanf(SFILE->f, "%*[^\n]");
      getc(SFILE->f);
    }
    fscanf(SFILE->f, "%*[^\n]");
    getc(SFILE->f);
  }  /* NBRF */

  /* Read in the sequence */
  J = 0;
  while (!BUFEOF(SFILE->f)) {
    while (!P_eoln(SFILE->f)) {
      CH = getc(SFILE->f);
      if (CH == '\n')
	CH = ' ';
      if (CH != 'X' && CH != '*' && CH != 'Z' && CH != 'B' && CH != 'R' &&
	  CH != 'K' && CH != 'H' && CH != 'E' && CH != 'D' && CH != 'W' &&
	  CH != 'Y' && CH != 'Q' && CH != 'N' && CH != 'C' && CH != 'T' &&
	  CH != 'S' && CH != 'P' && CH != 'F' && CH != 'M' && CH != 'I' &&
	  CH != 'L' && CH != 'V' && CH != 'A' && CH != 'G') {
	if (CH == ';') {   /* comment in input file */
	  fscanf(SFILE->f, "%*[^\n]");
	  getc(SFILE->f);
	}
	continue;
      }
      if (J >= MAXSEQ - MAXRANGE) {
	TOOLONG = true;
	continue;
      }
      J++;
      switch (CH) {

      case 'G':
	PR[J-1] = GLY;
	break;

      case 'A':
	PR[J-1] = ALA;
	break;

      case 'V':
	PR[J-1] = VALINE;
	break;

      case 'L':
	PR[J-1] = LEU;
	break;

      case 'I':
	PR[J-1] = ILE;
	break;

      case 'M':
	PR[J-1] = MET;
	break;

      case 'F':
	PR[J-1] = PHE;
	break;

      case 'P':
	PR[J-1] = PRO;
	break;

      case 'S':
	PR[J-1] = SER;
	break;

      case 'T':
	PR[J-1] = THR;
	break;

      case 'C':
	PR[J-1] = CYS;
	break;

      case 'N':
	PR[J-1] = ASN;
	break;

      case 'Q':
	PR[J-1] = GLN;
	break;

      case 'Y':
	PR[J-1] = TYR;
	break;

      case 'W':
	PR[J-1] = TRP;
	break;

      case 'D':
	PR[J-1] = ASP;
	break;

      case 'E':
	PR[J-1] = GLU;
	break;

      case 'H':
	PR[J-1] = HIS;
	break;

      case 'K':
	PR[J-1] = LYS;
	break;

      case 'R':
	PR[J-1] = ARG;
	break;

      case 'X':
	PR[J-1] = UNKX;
	break;

      case 'B':
	PR[J-1] = ASX;
	break;

      case 'Z':
	PR[J-1] = GLX;
	break;

      case '*':
	PR[J-1] = TERM;
	if (FILETYPE == 'N')
	  goto _L86;   /*ignore rest of file*/
	break;
      }
    }
    fscanf(SFILE->f, "%*[^\n]");
    getc(SFILE->f);
  }

  /* !!! */
_L86:
  /*CLOSE(SFILE); */
  /* branch destination for end of NBRF seq. */
  if (TOOLONG)
    printf(">>> WARNING! Sequence exceeds MAXSEQ-MAXRANGE. Sequence truncated.\n");
  *LEN = J;
}  /* READPRO */


#define WATER           18


/* END MODULE READPRO         VERSION= 'SUNMODS     Version  9/12/91'; */

/***********************************************************************/
/* Calculate amino acid composition and molecular weight of a protein. */
/***********************************************************************/
Static Void PROPARAM(P, SEQLEN, AACOMP, MW)
AMINOACID *P;
long SEQLEN;
long *AACOMP;
double *MW;
{
  long I;
  AMINOACID AA;

  /* Tabulate number of each amino acid in the protein */
  for (I = 0; I < SEQLEN; I++)
    AACOMP[(long)P[I]]++;

  /* Calculate molecular weight.*/
  for (AA = GLY; (long)AA <= (long)UNKX; AA = (AMINOACID)((long)AA + 1))
    *MW += AACOMP[(long)AA] * AAWT[(long)AA];
  *MW += WATER;   /* accounts for hydration of N & C termini */
}  /* MW */

#undef WATER


/**********************************************/
/* Print a report of the findings.            */
/**********************************************/
Static Void REPORT()
{
  printf("Type output filename:\n");
  GETFILE(&OUTFILE, 'O', &OFN);
  fprintf(OUTFILE.f, "%s\n", VERSION);
  WRITEWORD(&OUTFILE, NAME, (long)MAXWORD);
  fprintf(OUTFILE.f, "%10ld aa\n", SEQLEN);
  fprintf(OUTFILE.f, "Molecular weight:%10.1f\n\n", MW);
  fprintf(OUTFILE.f, "Nonpolar side chains:\n");
  fprintf(OUTFILE.f, "Gly(G)   %12ld (%5.3f)\n",
	  AACOMP[(long)GLY], (double)AACOMP[(long)GLY] / SEQLEN);
  fprintf(OUTFILE.f, "Ala(A)   %12ld (%5.3f)\n",
	  AACOMP[(long)ALA], (double)AACOMP[(long)ALA] / SEQLEN);
  fprintf(OUTFILE.f, "Val(V)   %12ld (%5.3f)\n",
	  AACOMP[(long)VALINE], (double)AACOMP[(long)VALINE] / SEQLEN);
  fprintf(OUTFILE.f, "Leu(L)   %12ld (%5.3f)\n",
	  AACOMP[(long)LEU], (double)AACOMP[(long)LEU] / SEQLEN);
  fprintf(OUTFILE.f, "Ile(I)   %12ld (%5.3f)\n",
	  AACOMP[(long)ILE], (double)AACOMP[(long)ILE] / SEQLEN);
  fprintf(OUTFILE.f, "Met(M)   %12ld (%5.3f)\n",
	  AACOMP[(long)MET], (double)AACOMP[(long)MET] / SEQLEN);
  fprintf(OUTFILE.f, "Phe(F)   %12ld (%5.3f)\n",
	  AACOMP[(long)PHE], (double)AACOMP[(long)PHE] / SEQLEN);
  fprintf(OUTFILE.f, "Pro(P)   %12ld (%5.3f)\n",
	  AACOMP[(long)PRO], (double)AACOMP[(long)PRO] / SEQLEN);
  fprintf(OUTFILE.f, "Neutral polar side chains:\n");
  fprintf(OUTFILE.f, "Ser(S)   %12ld (%5.3f)\n",
	  AACOMP[(long)SER], (double)AACOMP[(long)SER] / SEQLEN);
  fprintf(OUTFILE.f, "Thr(T)   %12ld (%5.3f)\n",
	  AACOMP[(long)THR], (double)AACOMP[(long)THR] / SEQLEN);
  fprintf(OUTFILE.f, "Cys(C)   %12ld (%5.3f)\n",
	  AACOMP[(long)CYS], (double)AACOMP[(long)CYS] / SEQLEN);
  fprintf(OUTFILE.f, "Asn(N)   %12ld (%5.3f)\n",
	  AACOMP[(long)ASN], (double)AACOMP[(long)ASN] / SEQLEN);
  fprintf(OUTFILE.f, "Gln(Q)   %12ld (%5.3f)\n",
	  AACOMP[(long)GLN], (double)AACOMP[(long)GLN] / SEQLEN);
  fprintf(OUTFILE.f, "Tyr(Y)   %12ld (%5.3f)\n",
	  AACOMP[(long)TYR], (double)AACOMP[(long)TYR] / SEQLEN);
  fprintf(OUTFILE.f, "Trp(W)   %12ld (%5.3f)\n",
	  AACOMP[(long)TRP], (double)AACOMP[(long)TRP] / SEQLEN);
  fprintf(OUTFILE.f, "Charged polar side chains:\n");
  fprintf(OUTFILE.f, "Asp(D)   %12ld (%5.3f)\n",
	  AACOMP[(long)ASP], (double)AACOMP[(long)ASP] / SEQLEN);
  fprintf(OUTFILE.f, "Glu(E)   %12ld (%5.3f)\n",
	  AACOMP[(long)GLU], (double)AACOMP[(long)GLU] / SEQLEN);
  fprintf(OUTFILE.f, "His(H)   %12ld (%5.3f)\n",
	  AACOMP[(long)HIS], (double)AACOMP[(long)HIS] / SEQLEN);
  fprintf(OUTFILE.f, "Lys(K)   %12ld (%5.3f)\n",
	  AACOMP[(long)LYS], (double)AACOMP[(long)LYS] / SEQLEN);
  fprintf(OUTFILE.f, "Arg(R)   %12ld (%5.3f)\n",
	  AACOMP[(long)ARG], (double)AACOMP[(long)ARG] / SEQLEN);
  fprintf(OUTFILE.f, "Other:\n");
  fprintf(OUTFILE.f, "Asx(B)   %12ld (%5.3f)\n",
	  AACOMP[(long)ASX], (double)AACOMP[(long)ASX] / SEQLEN);
  fprintf(OUTFILE.f, "Glx(Z)   %12ld (%5.3f)\n",
	  AACOMP[(long)GLX], (double)AACOMP[(long)GLX] / SEQLEN);
  fprintf(OUTFILE.f, "Unk(X)   %12ld (%5.3f)\n",
	  AACOMP[(long)UNKX], (double)AACOMP[(long)UNKX] / SEQLEN);
}  /*REPORT*/


/* ----------------------------------------------------------  */
/* ----------------- MAIN  PROCEDURE  -----------------------  */
main(argc, argv)
int argc;
Char *argv[];
{
  Char STR2[256];

  /* BEGIN MODULE STARTUP */
  /* Peform operations which must be done at beginning of main
     procedure. */
  /*!!!   TERMIN(input);    Open input for interactive use */
  /*!!!   TERMOUT(output);   "   output "      "        "  */
  PASCAL_MAIN(argc, argv);
  OUTFILE.f = NULL;
  *OUTFILE.name = '\0';
  SFILE.f = NULL;
  *SFILE.name = '\0';
  printf("%50s\n\n", VERSION);
  /* END MODULE STARTUP         VERSION= 'SUNMODS     Version  9/12/91'; */
  printf("%s\n\n", VERSION);

  printf("Enter sequence filename:\n");
  GETFILE(&SFILE, 'I', &SFN);
  NAME.LEN = 0;
  READPRO(&SFILE, SEQ, &SEQLEN, &NAME);
  if (NAME.LEN == 0) {
    printf("Type name for protein:\n");
    INPWORD(&NAME);
    scanf("%*[^\n]");
    getchar();
  }
  INITIALIZE();
  if (SEQ[SEQLEN-1] == TERM)
    SEQLEN--;
  PROPARAM(SEQ, SEQLEN, AACOMP, &MW);
  REPORT();
  /*!!!*/
  if (OUTFILE.f != NULL)
    fclose(OUTFILE.f);
  OUTFILE.f = NULL;
  if (SFILE.f != NULL)
    fclose(SFILE.f);
  exit(EXIT_SUCCESS);
}  /* PROSTAT */



/* End. */
