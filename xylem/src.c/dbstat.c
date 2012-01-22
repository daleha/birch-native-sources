/* Output from p2c 1.21alpha-07.Dec.93, the Pascal-to-C translator */
/* From input file "dbstat.p" */


/* ********************************************************  */
/*                                                           */
/*    DBSTAT   Version  2/ 1/90, Standard Pascal             */
/*             Brian Fristensky                              */
/*             Dept. of Plant Science                        */
/*             University of Manitoba                        */
/*             Winnipeg, MB R3T 2N2 CANADA                   */
/*                                                           */
/* Copyright (c) 1990 by Brian Fristensky.                   */
/* !!! in comment indicates feature which may need change.   */
/*  *******************************************************  */


#include <p2c.h>


/*!!! Some Pascals require file parameters in program heading */

#define VERSION         "DBSTAT       Version  2/ 1/90"


typedef enum {
  G, A, V, L, I, M, F, P, S, T, C, N, Q, Y, W, D, E, H, K, R, B, Z, X
} AMINOACID;
typedef long AAR[23];


/* Variables associated with the test protein */
Static long DBSIZE;
Static AAR AACOMP;   /* number of each amino acid in the protein */


/* **************************************************** */
/*              INITIALIZATION  PROCEDURES              */
/* **************************************************** */
Static Void INITIALIZE()
{
  AMINOACID AA1;

  for (AA1 = G; (long)AA1 <= (long)X; AA1 = (AMINOACID)((long)AA1 + 1))
    AACOMP[(long)AA1] = 0;
}  /* INITIALIZE */


/*  ******************************************* */
/*  Read a protein library and tabulate number  */
/*  of amino acids and composition of library.  */
/*  ******************************************* */
Static Void READLIB(SFILE, DBSIZE, AACOMP)
_TEXT *SFILE;
long *DBSIZE;
long *AACOMP;
{
  Char CH;

  /* Read in the sequence */
  *DBSIZE = 0;
  while (!BUFEOF(SFILE->f)) {
    while (!P_eoln(SFILE->f)) {
      CH = getc(SFILE->f);
      if (CH == '\n')
	CH = ' ';
      if (CH != 'X' && CH != 'Z' && CH != 'B' && CH != 'R' && CH != 'K' &&
	  CH != 'H' && CH != 'E' && CH != 'D' && CH != 'W' && CH != 'Y' &&
	  CH != 'Q' && CH != 'N' && CH != 'C' && CH != 'T' && CH != 'S' &&
	  CH != 'P' && CH != 'F' && CH != 'M' && CH != 'I' && CH != 'L' &&
	  CH != 'V' && CH != 'A' && CH != 'G') {
	if (CH == '>' || CH == ';') {   /* comment line */
	  fscanf(SFILE->f, "%*[^\n]");
	  getc(SFILE->f);
	}
	continue;
      }
      (*DBSIZE)++;
      switch (CH) {

      case 'G':
	AACOMP[(long)G]++;
	break;

      case 'A':
	AACOMP[(long)A]++;
	break;

      case 'V':
	AACOMP[(long)V]++;
	break;

      case 'L':
	AACOMP[(long)L]++;
	break;

      case 'I':
	AACOMP[(long)I]++;
	break;

      case 'M':
	AACOMP[(long)M]++;
	break;

      case 'F':
	AACOMP[(long)F]++;
	break;

      case 'P':
	AACOMP[(long)P]++;
	break;

      case 'S':
	AACOMP[(long)S]++;
	break;

      case 'T':
	AACOMP[(long)T]++;
	break;

      case 'C':
	AACOMP[(long)C]++;
	break;

      case 'N':
	AACOMP[(long)N]++;
	break;

      case 'Q':
	AACOMP[(long)Q]++;
	break;

      case 'Y':
	AACOMP[(long)Y]++;
	break;

      case 'W':
	AACOMP[(long)W]++;
	break;

      case 'D':
	AACOMP[(long)D]++;
	break;

      case 'E':
	AACOMP[(long)E]++;
	break;

      case 'H':
	AACOMP[(long)H]++;
	break;

      case 'K':
	AACOMP[(long)K]++;
	break;

      case 'R':
	AACOMP[(long)R]++;
	break;

      case 'B':
	AACOMP[(long)B]++;
	break;

      case 'Z':
	AACOMP[(long)Z]++;
	break;

      case 'X':
	AACOMP[(long)X]++;
	break;
      }
    }
    if (BUFEOF(SFILE->f))
      break;
    fscanf(SFILE->f, "%*[^\n]");
    getc(SFILE->f);
  }
}  /* READLIB */


/**********************************************/
/* Print a report of the findings.            */
/**********************************************/
Static Void REPORT(OUTFILE)
_TEXT *OUTFILE;
{
  fprintf(OUTFILE->f, "%s\n\n\n", VERSION);
  fprintf(OUTFILE->f, "Nonpolar side chains:\n");
  fprintf(OUTFILE->f, "Gly(G)%12ld (%5.3f)\n",
	  AACOMP[(long)G], (double)AACOMP[(long)G] / DBSIZE);
  fprintf(OUTFILE->f, "Ala(A)%12ld (%5.3f)\n",
	  AACOMP[(long)A], (double)AACOMP[(long)A] / DBSIZE);
  fprintf(OUTFILE->f, "Val(V)%12ld (%5.3f)\n",
	  AACOMP[(long)V], (double)AACOMP[(long)V] / DBSIZE);
  fprintf(OUTFILE->f, "Leu(L)%12ld (%5.3f)\n",
	  AACOMP[(long)L], (double)AACOMP[(long)L] / DBSIZE);
  fprintf(OUTFILE->f, "Ile(I)%12ld (%5.3f)\n",
	  AACOMP[(long)I], (double)AACOMP[(long)I] / DBSIZE);
  fprintf(OUTFILE->f, "Met(M)%12ld (%5.3f)\n",
	  AACOMP[(long)M], (double)AACOMP[(long)M] / DBSIZE);
  fprintf(OUTFILE->f, "Phe(F)%12ld (%5.3f)\n",
	  AACOMP[(long)F], (double)AACOMP[(long)F] / DBSIZE);
  fprintf(OUTFILE->f, "Pro(P)%12ld (%5.3f)\n",
	  AACOMP[(long)P], (double)AACOMP[(long)P] / DBSIZE);
  fprintf(OUTFILE->f, "Neutral polar side chains:\n");
  fprintf(OUTFILE->f, "Ser(S)%12ld (%5.3f)\n",
	  AACOMP[(long)S], (double)AACOMP[(long)S] / DBSIZE);
  fprintf(OUTFILE->f, "Thr(T)%12ld (%5.3f)\n",
	  AACOMP[(long)T], (double)AACOMP[(long)T] / DBSIZE);
  fprintf(OUTFILE->f, "Cys(C)%12ld (%5.3f)\n",
	  AACOMP[(long)C], (double)AACOMP[(long)C] / DBSIZE);
  fprintf(OUTFILE->f, "Asn(N)%12ld (%5.3f)\n",
	  AACOMP[(long)N], (double)AACOMP[(long)N] / DBSIZE);
  fprintf(OUTFILE->f, "Glu(Q)%12ld (%5.3f)\n",
	  AACOMP[(long)Q], (double)AACOMP[(long)Q] / DBSIZE);
  fprintf(OUTFILE->f, "Tyr(Y)%12ld (%5.3f)\n",
	  AACOMP[(long)Y], (double)AACOMP[(long)Y] / DBSIZE);
  fprintf(OUTFILE->f, "Trp(W)%12ld (%5.3f)\n",
	  AACOMP[(long)W], (double)AACOMP[(long)W] / DBSIZE);
  fprintf(OUTFILE->f, "Charged polar side chains:\n");
  fprintf(OUTFILE->f, "Asp(D)%12ld (%5.3f)\n",
	  AACOMP[(long)D], (double)AACOMP[(long)D] / DBSIZE);
  fprintf(OUTFILE->f, "Glu(E)%12ld (%5.3f)\n",
	  AACOMP[(long)E], (double)AACOMP[(long)E] / DBSIZE);
  fprintf(OUTFILE->f, "His(H)%12ld (%5.3f)\n",
	  AACOMP[(long)H], (double)AACOMP[(long)H] / DBSIZE);
  fprintf(OUTFILE->f, "Lys(K)%12ld (%5.3f)\n",
	  AACOMP[(long)K], (double)AACOMP[(long)K] / DBSIZE);
  fprintf(OUTFILE->f, "Arg(R)%12ld (%5.3f)\n",
	  AACOMP[(long)R], (double)AACOMP[(long)R] / DBSIZE);
  fprintf(OUTFILE->f, "Other:\n");
  fprintf(OUTFILE->f, "Asx(B)%12ld (%5.3f)\n",
	  AACOMP[(long)B], (double)AACOMP[(long)B] / DBSIZE);
  fprintf(OUTFILE->f, "Arg(Z)%12ld (%5.3f)\n",
	  AACOMP[(long)Z], (double)AACOMP[(long)Z] / DBSIZE);
  fprintf(OUTFILE->f, "Unk(X)%12ld (%5.3f)\n",
	  AACOMP[(long)X], (double)AACOMP[(long)X] / DBSIZE);
}  /*REPORT*/


/* ----------------------------------------------------------  */
/* ----------------- MAIN  PROCEDURE  -----------------------  */
main(argc, argv)
int argc;
Char *argv[];
{
  _TEXT TEMP;

  PASCAL_MAIN(argc, argv);
  INITIALIZE();
  TEMP.f = stdin;
  *TEMP.name = '\0';
  READLIB(&TEMP, &DBSIZE, AACOMP);

  TEMP.f = stdout;
  *TEMP.name = '\0';
  REPORT(&TEMP);
  exit(EXIT_SUCCESS);
}  /* DBSTAT */




/* End. */
