/*	pam.c	19-June-86
	copyright (c) 1987 William R. Pearson
	read in the alphabet and pam matrix data
	designed for universal matcher

	This version reads BLAST format (square) PAM files
*/

/*  $Id: apam.c 842 2011-10-17 16:05:50Z wrp $ */
/* $Revision: 842 $  */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "defs.h"
#include "param.h"

#define XTERNAL
#include "uascii.h"
#include "upam.h"
#undef XTERNAL

extern void alloc_pam (int d1, int d2, struct pstruct *ppst);

/* pam_opts -- modify PAM matrix (pamoff, -MS) if -MS or +off is part
   of PAM matrix name, e.g. MD20-MS or BL50+2 
*/

void
pam_opts(char *smstr, struct pstruct *ppst) {
  char *bp;

  ppst->pam_ms = 0;
  ppst->pamoff = 0;

  if ((bp=strchr(smstr,'-'))!=NULL) {
    if (!strncmp(bp+1,"MS",2) || !strncmp(bp+1,"ms",2)) {
      ppst->pam_ms = 1;
    }
    else {
      ppst->pamoff=atoi(bp+1);
    }
    *bp = '\0';
   }
  else if ((bp=strchr(smstr,'+'))!=NULL) {
    ppst->pamoff= -atoi(bp+1);
    *bp = '\0';
  }
}

/* modified 13-Oct-2005 to accomodate assymetrical matrices */

int
initpam (char *mfname, struct pstruct *ppst)
{
   char    line[512], *lp;
   int     i, j, iaa, pval;
   int *hsq, nsq;
   int *sascii;
   char *sq;
   int ess_tmp, max_val, min_val;
   int have_es = 0;
   FILE   *fmat;

   pam_opts(mfname, ppst);

   if ((fmat = fopen (mfname, "r")) == NULL)
   {
      printf ("***WARNING*** cannot open scoring matrix file %s\n", mfname);
      fprintf (stderr,"***WARNING*** cannot open scoring matrix file %s\n", mfname);
      return 0;
   }

/* 
   the size of the alphabet is determined in advance 
*/
   hsq = ppst->hsq;
   sq = ppst->sq;

   ppst->nt_align = (ppst->dnaseq == SEQT_DNA || ppst->dnaseq == SEQT_RNA);

/* 
     look for alphabet line, skipping the comments
     alphabet ends up in line[]
*/
   while (fgets (line, sizeof(line), fmat) != NULL && line[0]=='#');

   /* decide whether this is a protein or DNA matrix */
   if (ppst->nt_align) sascii = &nascii[0];
   else sascii = &aascii[0];

/*
  re-initialize sascii[] for matrix alphabet
*/

   /* save ',' value used by FASTS/FASTM/FASTF */
   ess_tmp = sascii[','];

/*	clear out sascii	*/
   for (i = 0; i <= AAMASK; i++) sascii[i] = NA;

/*	set end of line stop	*/
   sascii[0] = sascii['\r'] = sascii['\n'] = EL;

   sascii[','] = ess_tmp;

/* read the alphabet - determine alphabet nsq */
   sq[0] = '\0';
   for (i = 0, nsq = 1; line[i]; i++) {
     if (line[i] == '*') have_es = 1;
     if (line[i] > ' ') sq[nsq++] = toupper (line[i]);
   }
   sq[nsq]='\0';
   nsq--;

   if (nsq < 20) {
     ppst->nt_align = 1;
     if (ppst->dnaseq <= SEQT_PROT) {
       sascii = &nascii[0];
       ppst->dnaseq = SEQT_DNA;
     }
     memcpy(qascii, sascii, sizeof(qascii));
   }

/* set end of sequence stop */
   fprintf(stderr,"sq[%d]: %s\n",nsq,sq+1);

/* initialize sascii */
   for (iaa = 1; iaa <= nsq; iaa++) {
     sascii[sq[iaa]] = iaa;
   }
   if (ppst->dnaseq==SEQT_DNA) {
     sascii['U'] = sascii['T'];
     sascii['u'] = sascii['t'];
   }
   else if (ppst->dnaseq==SEQT_RNA) {
     sascii['T'] = sascii['U'];
     sascii['t'] = sascii['u'];
   }

/* finished with sascii[] */

/* setup hnt (ambiguous nt hash) values */
   hsq[0] = 0;
   for (iaa = 1; iaa <= nsq; iaa++) {
     hsq[iaa]=iaa;
   }
   if (ppst->nt_align || nsq < 20) {	/* DNA ambiguitities */
     hsq[sascii['R']]=hsq[sascii['M']]=hsq[sascii['W']]=hsq[sascii['A']];
     hsq[sascii['D']]=hsq[sascii['H']]=hsq[sascii['V']]=hsq[sascii['A']];
     hsq[sascii['N']]=hsq[sascii['X']]=hsq[sascii['A']];
     hsq[sascii['Y']]=hsq[sascii['S']]=hsq[sascii['B']]=hsq[sascii['C']];
     hsq[sascii['K']]=hsq[sascii['G']];
   }
   else 	/* protein ambiguities */
     if (ppst->dnaseq == SEQT_UNK || ppst->dnaseq == SEQT_PROT ||
	 (ppst->nsq >= 20 && ppst->nsq <= 24)) {
     hsq[sascii['B']] = hsq[sascii['N']];
     hsq[sascii['Z']] = hsq[sascii['E']];
     hsq[sascii['X']] = hsq[sascii['A']];
   }
   /* here if non-DNA, non-protein sequence */
   else ppst->dnaseq = SEQT_OTHER;

/*
  check for 2D pam  - if not found, allocate it
*/

   if (!ppst->have_pam2) {
     alloc_pam (MAXSQ, MAXSQ, ppst);
     ppst->have_pam2 = 1;
   }

/*
  read the scoring matrix values
*/

   max_val = -1;
   min_val =  1;
   for (j=0; j <= nsq; j++) {
     ppst->pam2[0][0][j] = ppst->pam2[0][j][0] = -BIGNUM;
   }
   for (iaa = 1; iaa <= nsq; iaa++) {	/* read pam value line */
     if (fgets(line,sizeof(line),fmat)==NULL) {
       fprintf (stderr," error reading pam line: %s\n",line);
       exit (1);
     }
     /*     fprintf(stderr,"%d/%d %s",iaa,nsq,line); */
     strtok(line," \t\n");		/* skip the letter (residue) */

     for (j = 1; j <= nsq; j++) {	/* iaa limits to triangle */
       lp=strtok(NULL," \t\n");		/* get the number string */
       pval=ppst->pam2[0][iaa][j]=atoi(lp);	/* convert to integer */
       if (pval > max_val) max_val = pval;
       if (pval < min_val) min_val = pval;
     }
   }

   /* this should only be done for proteins */
   if (!ppst->nt_align && have_es==0) {
     sascii['*']=nsq;
     nsq++;
     sq[nsq]='*';
     sq[nsq+1]='\0';
     ppst->pam2[0][nsq][0] = ppst->pam2[0][0][nsq] = -BIGNUM;
     for (j=1; j<=nsq; j++) ppst->pam2[0][nsq][j]= -1;
     ppst->pam2[0][nsq][nsq]= max_val/2;
   }

   ppst->sqx[0]='\0';	/* initialize sqx[] */
   for (i=1; i<= nsq; i++) {
     ppst->sqx[i] = sq[i];
     ppst->sqx[i+nsq] = tolower(sq[i]);
   }
   /*
     if (sascii[aa[i]] < NA && sq[i] >= 'A' && sq[i] <= 'Z')
        sascii[aa[i] - 'A' + 'a'] = sascii[aa[i]]+nsq;
   */

   ppst->nsq = nsq;	/* save new nsq */
   ppst->nsqx = nsq*2;	/* save new nsqx */

   ppst->pam_h = max_val;
   ppst->pam_l = min_val;

   strncpy (ppst->pamfile, mfname, MAX_FN);
   ppst->pamfile[MAX_FN-1]='\0';
   strncpy (ppst->pam_name, ppst->pamfile, MAX_FN);

   if (ppst->pam_ms) {
     strncat(ppst->pam_name,"-MS",MAX_FN-strlen(ppst->pam_name)-1);
   }
   ppst->pam_name[MAX_FN-1]='\0';
   fclose (fmat);
   return 1;
}

/* make a DNA scoring from +match/-mismatch values */

void mk_n_pam(int *arr,int siz, int mat, int mis)
{
  int i, j, k;
  /* current default match/mismatch values */
  int max_mat = +5;
  int min_mis = -4;
  float f_val, f_scale;
  
  f_scale = (float)(mat - mis)/(float)(max_mat - min_mis);
  
  k = 0;
  for (i = 0; i<nnt-1; i++)
    for (j = 0; j <= i; j++ ) {
      if (arr[k] == max_mat) arr[k] = mat;
      else if (arr[k] == min_mis) arr[k] = mis;
      else if (arr[k] != -1) { 
	f_val = (arr[k] - min_mis)*f_scale + 0.5;
	arr[k] = f_val + mis;
      }
      k++;
    }
}

int
standard_pam(char *smstr, struct pstruct *ppst, int del_set, int gap_set) {

  struct std_pam_str *std_pam_p;

  pam_opts(smstr, ppst);

  for (std_pam_p = std_pams; std_pam_p->abbrev[0]; std_pam_p++ ) {
    if (strcmp(smstr,std_pam_p->abbrev)==0) {
      pam = std_pam_p->pam;
      strncpy(ppst->pam_name,std_pam_p->name,MAX_FN);
      ppst->pam_name[MAX_FN-1]='\0';
      if (ppst->pam_ms) {
	strncat(ppst->pam_name,"-MS",MAX_FN-strlen(ppst->pam_name)-1);
      }
      ppst->pam_name[MAX_FN-1]='\0';
#ifdef OLD_FASTA_GAP
      if (!del_set) ppst->gdelval = std_pam_p->gdel+std_pam_p->ggap;
#else
      if (!del_set) ppst->gdelval = std_pam_p->gdel;
#endif
      if (!gap_set) ppst->ggapval = std_pam_p->ggap;
      ppst->pamscale = std_pam_p->scale;
      return 1;
    }
  }
  return 0;
}

/* scan through the sorted (by decreasing entropy) std_pams[] array
   for a scoring matrix with enough entropy
*/
int
min_pam_bits(int n0_eff, double bit_thresh, struct pstruct *ppst, int del_set, int gap_set) {
  struct std_pam_str *std_pam_p;
  int curr_pam_idx = 0;

  pam_opts(ppst->pamfile, ppst);

  /* get the index for the current (standard) pam file */
  for (curr_pam_idx = 0; std_pams[curr_pam_idx].abbrev[0]; curr_pam_idx++) {
    if (strcmp(ppst->pamfile,std_pams[curr_pam_idx].abbrev)==0) break;
  }

  for ( ; curr_pam_idx > 0 ; curr_pam_idx-- ) {
    if (n0_eff * std_pams[curr_pam_idx].entropy >= bit_thresh) goto new_pam;
  }
  return 0;

 new_pam:
  std_pam_p = &std_pams[curr_pam_idx];

  pam = std_pam_p->pam;
  strncpy(ppst->pam_name,std_pam_p->name,MAX_FN);
  if (ppst->pam_ms) {
    strncat(ppst->pam_name,"-MS",MAX_FN-strlen(ppst->pamfile)-1);
  }
  ppst->pam_name[MAX_FN-1]='\0';
#ifdef OLD_FASTA_GAP
  if (!del_set) ppst->gdelval = std_pam_p->gdel+std_pam_p->ggap;
#else
  if (!del_set) ppst->gdelval = std_pam_p->gdel;
#endif
  if (!gap_set) ppst->ggapval = std_pam_p->ggap;
  ppst->pamscale = std_pam_p->scale;
  return 1;
}

/* ESS must match uascii.h */
#define ESS 49

/* build_xascii is only used for SEQT_UNK - it replaces the default
   input mapping (aascii[])  with a mapping that preserves any letter
   in either the aax[], ntx[], or othx[] alphabets, or in
   save_str[]. othx[] was added to support letters that are mapped,
   but are not (yet) in aax[], e.g. 'OoUu'.  Because build_xascii
   makes a qascii[] that is all ascii characters with values > '@',
   these values must be replaced using either aascii[] or nascii[] and
   the initial query sequence re-coded.
 */
void
build_xascii(int *qascii, char *save_str) {
  int i, max_save;
  int comma_val, term_val;
  int save_arr[MAX_SSTR];

  comma_val = qascii[','];
  term_val = qascii['*'];

  /* preserve special characters */
  for (i=0; i < MAX_SSTR && save_str[i]; i++ ) {
    save_arr[i] = qascii[save_str[i]];
  }
  max_save = i;

  for (i=1; i<128; i++) {
    qascii[i]=NA;
  }

  /* range of values in aax, ntx is from 1..naax,nntx - 
     do not zero-out qascii[0] - 9 Oct 2002 */

  for (i=1; i<naax; i++) {
    qascii[aax[i]]=aax[i];
  }

  for (i=1; i<nntx; i++) {
    qascii[ntx[i]]=ntx[i];
  }

  /* put it letters that are not in other alphabets because they are
     mapped */
  for (i=1; i<nothx; i++) {
    qascii[othx[i]]=othx[i];
  }

  qascii['\n']=qascii['\r']=qascii[0] = EL;

  qascii[','] = comma_val;
  qascii['*'] = term_val;

  for (i=0; i < max_save; i++) {
    qascii[save_str[i]]=save_arr[i];
  }
}

/* init_ascii()

   checks for lower case letters in *sq array;
   if not present, map lowercase to upper

*/
void
init_ascii(int is_ext, int *xascii, int p_nsq, int is_dna) {

  int isq, have_lc;
  char *sq, term_char;
  int nsq;
  
  if (is_dna==SEQT_UNK) return;

  term_char = xascii['*'];

  if (is_dna==SEQT_DNA || is_dna == SEQT_RNA) {
    if (is_ext) { 
      sq = &ntx[0];
      nsq = nntx;
    }
    else {sq = &nt[0]; nsq = nnt;}
  }
  else {
    if (is_ext) { sq = &aax[0]; nsq = naax; }
    else {sq = &aa[0]; nsq = naa;}
  }

  /* initialize xascii from sq[], checking for lower-case letters */
  /* this code guarantees that all characters in sq are represented in
     xascii[], but it does not guarantee that everything else in
     xascii[], particularly xascii[O,U,J], have appropriate values
 */
  have_lc = 0;
  for (isq = 1; isq <= nsq; isq++) {
     xascii[sq[isq]] = isq;
     if (sq[isq] >= 'a' && sq[isq] <= 'z') have_lc = 1;
  }
  
  /* no lower case letters in alphabet, map lower case to upper */
  if (have_lc != 1) { 
    for (isq = 1; isq <= nsq; isq++) {
      if (sq[isq] >= 'A' && sq[isq] <= 'Z') xascii[sq[isq]-'A'+'a'] = isq;
    }
    if (is_dna==1) xascii['u'] = xascii['t'];
  }

  xascii['*']=term_char;
}

void
validate_novel_aa(int *xascii, int p_nsq, int dnaseq) {
  int isq, err_p_nsq_limit;
  /* these checks need to be done after xascii[] has been
     re-initialized */

  if (dnaseq != SEQT_DNA && dnaseq!=SEQT_RNA) {
    if (xascii['O'] > p_nsq || xascii['o'] > p_nsq) { xascii['O'] = xascii['K']; xascii['o'] = xascii['k'];}
    if (xascii['U'] > p_nsq || xascii['u'] > p_nsq) { xascii['U'] = xascii['C']; xascii['u'] = xascii['c'];}
    if (xascii['J'] > p_nsq || xascii['j'] > p_nsq) { xascii['J'] = xascii['L']; xascii['j'] = xascii['l'];}
  }

  /* one final check for characters out of range (> nsq)*/
  err_p_nsq_limit = 0;
  for (isq = 'A'; isq <= 'Z'; isq++) {
    if (xascii[isq] < NA && xascii[isq] > p_nsq) {
      fprintf(stderr, " *** ERROR *** xascii['%c']:%d > %d\n",isq, xascii[isq], p_nsq);
      err_p_nsq_limit = 1;
    }
  }
  if (err_p_nsq_limit) {exit(1);}
}


void 
print_pam(struct pstruct *ppst) {
  int i, nsq, ip;
  char *sq;

  fprintf(stderr," ext_sq_set: %d\n",ppst->ext_sq_set);

  nsq = ppst->nsq;
  ip = 0;
  sq = ppst->sq;

  fprintf(stderr," sq[%d]: %s\n",nsq, sq);

  if (ppst->ext_sq_set) {
    nsq = ppst->nsqx;
    ip = 1;
    sq = ppst->sqx;
    fprintf(stderr," sq[%d]: %s\n",nsq, sq);
  }

  for (i=1; i<=nsq; i++) {
    fprintf(stderr," %c:%c - %3d\n",sq[i], sq[i], ppst->pam2[ip][i][i]);
  }
}
