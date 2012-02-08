/* copyright (c) 1996, 1997, 1998, 1999 William R. Pearson and the
   U. of Virginia */

/*  $Id: compacc2.c 872 2011-10-31 18:27:04Z wrp $ */
/* $Revision: 872 $  */

/* Concurrent read version */

#include <stdio.h>
#include <stdlib.h>
#if defined(UNIX)
#include <unistd.h>
#endif
#if defined(UNIX) || defined(WIN32)
#include <sys/types.h>
#endif

#include <limits.h>
#include <float.h>

#include <string.h>
#include <time.h>
#include <math.h>

#include "defs.h"
#include "param.h"
#include "structs.h"

#include "mm_file.h"
#include "best_stats.h"

#define XTERNAL
#include "uascii.h"
#include "upam.h"
#undef XTERNAL

extern void abort ();

#include "drop_func.h"	/* get init_work() */

void revcomp(unsigned char *seq, int n, int *c_nt);
extern void qshuffle(unsigned char *aa0, int n0, int nm0, void *);
#ifdef DEBUG
unsigned long adler32(unsigned long, const unsigned char *, unsigned int);
#endif

void
print_sum(FILE *fd, struct db_str *qtt, struct db_str *ntt, int in_mem, long mem_use);
int
check_seq_range(unsigned char *aa1b, int n1, int nsq, char *str);
/* print timing information */
extern void ptime (FILE *, long);

/* this function consolidates code in comp_lib4.c for non-threaded, and in
   work_thr2.c (threads) and work_comp2.c (worker nodes)
*/

void
init_aa0(unsigned char **aa0, int n0, int nm0,
	 unsigned char **aa0s, unsigned char **aa1s, 
	 int qframe, int qshuffle_flg, int max_tot,
	 struct pstruct *ppst, void **f_str, void **qf_str,
	 void *my_rand_state) {
  int id;

  /* note that aa[5,4,3,2] are never used, but are provided so that frame
     can range from 0 .. 5; likewise for f_str[5..2] */

  aa0[5] = aa0[4] = aa0[3] = aa0[2] = aa0[1] = aa0[0];

  /* zero out for SSE2/ALTIVEC -- make sure this is ALWAYS done */
  for (id=0; id < SEQ_PAD; id++) aa0[0][n0+id] = '\0';

  init_work (aa0[0], n0, ppst, &f_str[0]);
  f_str[5] = f_str[4] = f_str[3] = f_str[2] = f_str[1] = f_str[0];

  if (qframe == 2) {
    if ((aa0[1]=(unsigned char *)calloc((size_t)n0+2+SEQ_PAD,sizeof(unsigned char)))==NULL) {
      fprintf(stderr," cannot allocate aa01[%d]\n", n0);
    }
    *aa0[1]='\0';
    aa0[1]++;
    memcpy(aa0[1],aa0[0],n0+1);
    /* for ALTIVEC/SSE2, must pad with 16 NULL's */
    for (id=0; id<SEQ_PAD; id++) {aa0[1][n0+id]=0;}
    revcomp(aa0[1],n0,ppst->c_nt);
    init_work (aa0[1], n0, ppst, &f_str[1]);
  }

  if (qshuffle_flg) {
    if ((*aa0s=(unsigned char *)calloc(n0+2+SEQ_PAD,sizeof(char)))==NULL) {
      fprintf(stderr,"cannot allocate aa0s[%d]\n",n0+2);
      exit(1);
    }
    **aa0s='\0';
    (*aa0s)++;
    memcpy(*aa0s,aa0[0],n0);
    qshuffle(*aa0s,n0,nm0, my_rand_state);
    /* for SSE2/ALTIVEC, must pad with 16 NULL's */
    for (id=0; id<SEQ_PAD; id++) {(*aa0s)[n0+id]=0;}
    init_work (*aa0s, n0, ppst, qf_str);
  }

  /* always allocate shuffle space */
  if((*aa1s=calloc(max_tot+1,sizeof(char))) == NULL) {
    fprintf(stderr,"unable to allocate shuffled library sequence [%d]\n", max_tot);
    exit(1);
  }
  else {
    **aa1s=0;
    (*aa1s)++;
  }
}

/* because it is used to pre-allocate space, maxn has various
   constraints.  For "simple" comparisons, it is simply the length of
   the longest library sequence.  But for translated comparisons, it
   must be 3 or 6X the length of the query sequence. 

   In addition, however, it can be reduced to make certain that
   sequences are read in smaller chunks.  And, maxn affect how large
   overlaps must be when sequences are read in chunks.
*/

int
reset_maxn(struct mngmsg *m_msp, int over_len, int maxn) {

  /* reduce maxn if requested */
  if (m_msp->ldb_info.maxn > 0 && m_msp->ldb_info.maxn < maxn) maxn = m_msp->ldb_info.maxn;

  if (m_msp->qdnaseq==m_msp->ldb_info.ldnaseq || m_msp->qdnaseq==SEQT_DNA ||
      m_msp->qdnaseq == SEQT_RNA) {/* !TFAST - either FASTA or FASTX */

    if (m_msp->n0 > m_msp->max_tot - over_len) {
      fprintf(stderr," query sequence is too long %d > %d - %d %s\n",
	      m_msp->n0,
	      m_msp->max_tot, over_len,
	      m_msp->sqnam);
      exit(1);
    }

    m_msp->ldb_info.l_overlap = over_len;
    m_msp->ldb_info.maxt3 = maxn-m_msp->ldb_info.l_overlap;
  }
  else {	/* is TFAST */
    if (m_msp->n0 > MAXTST) {
      fprintf(stderr," query sequence is too long %d %s\n",m_msp->n0,m_msp->sqnam);
      exit(1);
    }

    if (m_msp->n0*3 > maxn ) {	/* n0*3 for the three frames - this
				   will only happen if maxn has been
				   set low manually */

      if (m_msp->n0*4+2 < m_msp->max_tot) { /* m_msg0*3 + m_msg0 */
	fprintf(stderr,
		" query sequence too long for library segment: %d - resetting to %d\n",
	      maxn,m_msp->n0*3);
	maxn = m_msp->ldb_info.maxn = m_msp->n0*3;
      }
      else {
	fprintf(stderr," query sequence too long for translated search: %d * 4 > %d %s\n",
	      m_msp->n0,maxn, m_msp->sqnam);
	exit(1);
      }
    }

    /* set up some constants for overlaps */
    m_msp->ldb_info.l_overlap = 3*over_len;
    m_msp->ldb_info.maxt3 = maxn-m_msp->ldb_info.l_overlap-3;
    m_msp->ldb_info.maxt3 -= m_msp->ldb_info.maxt3%3;
    m_msp->ldb_info.maxt3++;

    maxn = maxn - 3; maxn -= maxn%3; maxn++;
  }
  return maxn;
}


int
scanseq(unsigned char *seq, int n, char *str) {
  int tot,i;
  char aaray[128];		/* this must be set > nsq */
	
  for (i=0; i<128; i++)  aaray[i]=0;
  for (i=0; (size_t)i < strlen(str); i++) aaray[qascii[str[i]]]=1;
  for (i=tot=0; i<n; i++) tot += aaray[seq[i]];
  return tot;
}

/* subs_env takes a string, possibly with ${ENV}, and looks up all the
   potential environment variables and substitutes them into the
   string */

void subs_env(char *dest, char *src, int dest_size) {
  char *last_src, *bp, *bp1;

  last_src = src;

  if ((bp = strchr(src,'$'))==NULL) {
    strncpy(dest, src, dest_size);
    dest[dest_size-1] = '\0';
  }
  else {
    *dest = '\0';
    while (strlen(dest) < dest_size-1 && bp != NULL ) {
      /* copy stuff before ${*/
      *bp = '\0';
      strncpy(dest, last_src, dest_size);
      *bp = '$';

      /* copy ENV */
      if (*(bp+1) != '{') {
	strncat(dest, "$", dest_size - strlen(dest) -1);
	dest[dest_size-1] = '\0';
	bp += 1;
      }
      else {	/* have  ${ENV} - put it in */
	if ((bp1 = strchr(bp+2,'}'))==NULL) {
	  fprintf(stderr, "Unterminated ENV: %s\n",src);
	  break;
	}
	else {
	  *bp1 = '\0';
	  if (getenv(bp+2)!=NULL) {
	    strncat(dest, getenv(bp+2), dest_size - strlen(dest) - 1);
	    dest[dest_size-1] = '\0';
	    *bp1 = '}';
	  }
	  bp = bp1+1;	/* bump bp even if getenv == NULL */
	}
      }
      last_src = bp;

      /* now get the next ${ENV} if present */
      bp = strchr(last_src,'$');
    }
    /* now copy the last stuff */
    strncat(dest, last_src, dest_size - strlen(dest) - 1);
    dest[dest_size-1]='\0';
  }
}


void
selectbest(struct beststr **bptr, int k, int n)	/* k is rank in array */
{
  int v, i, j, l, r;
  struct beststr *tmptr;

  l=0; r=n-1;

  while ( r > l ) {
    v = bptr[r]->rst.score[0];
    i = l-1;
    j = r;
    do {
      while (bptr[++i]->rst.score[0] > v) ;
      while (bptr[--j]->rst.score[0] < v) ;
      tmptr = bptr[i]; bptr[i]=bptr[j]; bptr[j]=tmptr;
    } while (j > i);
    bptr[j]=bptr[i]; bptr[i]=bptr[r]; bptr[r]=tmptr;
    if (i>=k) r = i-1;
    if (i<=k) l = i+1;
  }
}

void
selectbestz(struct beststr **bptr, int k, int n)	/* k is rank in array */
{
  int i, j, l, r;
  struct beststr *tmptr;
  double v;

  l=0; r=n-1;

  while ( r > l ) {
    v = bptr[r]->zscore;
    i = l-1;
    j = r;
    do {
      while (bptr[++i]->zscore > v) ;
      while (bptr[--j]->zscore < v) ;
      tmptr = bptr[i]; bptr[i]=bptr[j]; bptr[j]=tmptr;
    } while (j > i);
    bptr[j]=bptr[i]; bptr[i]=bptr[r]; bptr[r]=tmptr;
    if (i>=k) r = i-1;
    if (i<=k) l = i+1;
  }
}

/* improved shellsort with high-performance increments */
/*
shellsort(itemType a[], int l, int r)
{ int i, j, k, h; itemType v;
 int incs[16] = { 1391376, 463792, 198768, 86961, 33936,
		  13776, 4592, 1968, 861, 336, 
		  112, 48, 21, 7, 3, 1 };
 for ( k = 0; k < 16; k++)
   for (h = incs[k], i = l+h; i <= r; i++) { 
       v = a[i]; j = i;
       while (j > h && a[j-h] > v) {
         a[j] = a[j-h]; j -= h;
       }
       a[j] = v; 
     } 
}
*/

/* ?improved? version of sortbestz using optimal increments and fewer
   exchanges */
void sortbestz(struct beststr **bptr, int nbest)
{
  int gap, i, j, k;
  struct beststr *tmp;
  double v;
  int incs[14] = { 198768, 86961, 33936,
		   13776, 4592, 1968, 861, 336, 
		   112, 48, 21, 7, 3, 1 };

  for ( k = 0; k < 14; k++) {
    gap = incs[k];
    for (i=gap; i < nbest; i++) {
      tmp = bptr[i];
      j = i;
      v = bptr[i]->zscore;
      while ( j >= gap && bptr[j-gap]->zscore < v) {
	bptr[j] = bptr[j - gap];
	j -= gap;
      }
      bptr[j] = tmp;
    }
  }
}

/* sort based on sequence index */
void sortbesti(struct beststr **bptr, int nbest)
{
  int gap, i, j, k;
  struct beststr *tmp;
  double v;
  int incs[12] = { 33936, 13776, 4592, 1968, 861, 336, 
		   112, 48, 21, 7, 3, 1 };

  for ( k = 0; k < 12; k++) {
    gap = incs[k];
    for (i=gap; i < nbest; i++) {
      tmp = bptr[i];
      j = i;
      v = bptr[i]->seq->index;
      while ( j >= gap && bptr[j-gap]->seq->index < v) {
	bptr[j] = bptr[j - gap];
	j -= gap;
      }
      bptr[j] = tmp;
    }
  }
}

void
sortbeste(struct beststr **bptr, int nbest)
{
  int gap, i, j, k;
  struct beststr *tmp;
  double v;
  int incs[14] = { 198768, 86961, 33936,
		   13776, 4592, 1968, 861, 336, 
		   112, 48, 21, 7, 3, 1 };

  for ( k = 0; k < 14; k++) {
    gap = incs[k]; 
    for (i=gap; i < nbest; i++) {
      j = i;
      tmp = bptr[i];
      v = tmp->rst.escore;
      while ( j >= gap && bptr[j-gap]->rst.escore > v) {
	bptr[j] = bptr[j - gap];
	j -= gap;
      }
      bptr[j] = tmp;
    }
  }

  /* sometimes there are many high scores with E()==0.0, sort
     those by z() score */

  j = 0;
  while (j < nbest && bptr[j]->rst.escore <= 2.0*DBL_MIN ) {j++;}
  if (j > 1) sortbestz(bptr,j);
}

extern char *prog_func;
extern char *verstr, *iprompt0, *refstr, *mp_verstr;
extern long tstart, tscan, tprev, tdone;	/* Timing */
#ifdef COMP_MLIB
extern long ttscan, ttdisp;
#endif
extern time_t tdstart, tddone;

/* ****************************************************************
   print command line arguments (argv_line)
   possibly HTML header
   !BLAST 
     please cite
     version
   BLAST
     Reference version
**************************************************************** */
void
print_header1(FILE *fd, const char *argv_line,
	      const struct mngmsg *m_msp, const struct pstruct *ppst) {

#ifdef PGM_DOC
  if (!(m_msp->markx & (MX_M8OUT+MX_MBLAST2))) fprintf(fd, "#%s\n",argv_line);
#endif

  if (m_msp->markx & MX_M11OUT) {
    fprintf(fd, "#:lav\n\nd {\n   \"%s\"\n}\n",argv_line+1);
  }

  if (m_msp->markx & MX_HTML) {
#ifdef HTML_HEAD    
    fprintf(fd,"<html>\n<head>\n<title>%s Results</title>\n</head>\n<body>\n",prog_func);
#endif
    fprintf(fd,"<pre>\n");
  }

  if (m_msp->std_output) {
    fprintf(fd,"%s\n",iprompt0);
    if (refstr != NULL && refstr[0] != '\0') {
      fprintf(fd," version %s%s\nPlease cite:\n %s\n",verstr,mp_verstr,refstr);
    }
    else {
      fprintf(fd," version %s%s\n",verstr,mp_verstr);
    }
  }

  if (m_msp->markx & MX_MBLAST2) {
    if (refstr != NULL && refstr[0] != '\0') {
      fprintf(fd,"%s %s%s\n\nReference: %s\n\n",prog_func, verstr,mp_verstr,refstr);
    }
    else {
      fprintf(fd,"%s %s%s\n\n",prog_func, verstr,mp_verstr);
    }
  }

  fflush(fd);
}

/* ****************************************************************
   MX_HTML: <pre>
   Query:
     1>>>accession description # aa
   Annotation:
   Library:
**************************************************************** */
void
print_header2(FILE *fd, int qlib, char *info_qlabel, unsigned char **aa0,
	      const struct mngmsg *m_msp, const struct pstruct *ppst,
	      const char * info_lib_range_p) {
  int j;
  char tmp_str[MAX_STR];
  double db_tt;
  
  /* if (m_msp->markx & MX_HTML) fputs("<pre>\n",fd); */

  if (m_msp->std_output) {
    if (qlib==1) {
      fprintf(fd,"Query: %s\n", m_msp->tname);
    }

    if (m_msp->qdnaseq == SEQT_DNA || m_msp->qdnaseq == SEQT_RNA) {
      strncpy(tmp_str,(m_msp->qframe==1)? " (forward-only)" : "\0",sizeof(tmp_str));
      tmp_str[sizeof(tmp_str)-1]='\0';
    }
    else tmp_str[0]='\0';

    fprintf(fd,"%3d>>>%s%s\n", qlib,
	   m_msp->qtitle,
	   (m_msp->revcomp ? " (reverse complement)" : tmp_str));
    /* check for annotation */
    if (m_msp->ann_flg && m_msp->aa0a != NULL) {
      fprintf(fd,"Annotation: ");
      for (j=0; j<m_msp->n0; j++) {
	if (m_msp->aa0a[j] && m_msp->ann_arr[m_msp->aa0a[j]] != ' ' ) {
	  fprintf(fd,"|%ld:%c%c",
		 j+m_msp->q_off,m_msp->ann_arr[m_msp->aa0a[j]],ppst->sq[aa0[0][j]]);
	}
      }
      fprintf(fd,"\n");
    }

    fprintf(fd,"Library: %s%s\n", m_msp->ltitle,info_lib_range_p);

    if (m_msp->db.carry==0) {
      fprintf(fd, "  %7ld residues in %5ld sequences\n", m_msp->db.length, m_msp->db.entries);
    }
    else {
      db_tt = (double)m_msp->db.carry*(double)LONG_MAX + (double)m_msp->db.length;
      fprintf(fd, "  %.0f residues in %5ld library sequences\n", db_tt, m_msp->db.entries);
    }

  }
  else {
    if (m_msp->markx & (MX_M8OUT + MX_M8COMMENT)) {
      fprintf(fd,"# %s %s%s\n",prog_func,verstr,mp_verstr);
      fprintf(fd,"# Query: %s\n",m_msp->qtitle);
      fprintf(fd,"# Database: %s\n",m_msp->ltitle);
    }
  }
  if (m_msp->markx & MX_HTML) fputs("</pre>\n",fd);
  fflush(fd);
}

/* **************************************************************** */
/*   before showbest                                                */
/* **************************************************************** */
void print_header3(FILE *fd, int qlib, struct mngmsg *m_msp, struct pstruct *ppst) {

    if (m_msp->markx & MX_MBLAST2) {
      if (qlib == 1) {
	fprintf(fd, "\nDatabase: %s\n     %12ld sequences; %ld total letters\n\n\n",
		m_msp->ltitle, m_msp->db.entries, m_msp->db.length);
      }
      fprintf(fd, "\nQuery= %s\nLength=%d\n", m_msp->qtitle, m_msp->n0);
    }
}


/* **************************************************************** */
/* alignment tranistion                                             */
/* **************************************************************** */
void print_header4(FILE *fd, char *info_qlabel, char *argv_line, char *info_gstring3, char *info_hstring_p[2],
		   struct mngmsg *m_msp, struct pstruct *ppst) {

	if (m_msp->std_output && (m_msp->markx & (MX_AMAP+ MX_HTML + MX_M9SUMM)) && !(m_msp->markx & MX_M10FORM)) {
	  fprintf(fd,"\n>>>%s%s, %d %s vs %s library\n",
		  info_qlabel,(m_msp->revcomp ? "_rev":"\0"), m_msp->n0,
		  m_msp->sqnam,m_msp->lname);
	}

	if (m_msp->markx & MX_M10FORM) {
	  fprintf(fd,"\n>>>%s%s, %d %s vs %s library\n",
		  info_qlabel,(m_msp->revcomp ? "-":"\0"), m_msp->n0, m_msp->sqnam,
		  m_msp->lname);
	  fprintf(fd,"; pg_name: %s\n",m_msp->pgm_name);
	  fprintf(fd,"; pg_ver: %s%s\n",verstr,mp_verstr);
	  fprintf(fd,"; pg_argv: %s",argv_line);
	  fputs(info_gstring3,fd);
	  fputs(info_hstring_p[0],fd);
	  fputs(info_hstring_p[1],fd);
	}
}

void print_header4a(FILE *outfd, struct mngmsg *m_msp) {
  if (!(m_msp->markx & MX_M8OUT) && (m_msp->markx & (MX_M10FORM+MX_M9SUMM)) && m_msp->show_code != SHOW_CODE_ID) {
    fprintf(outfd,">>><<<\n");
  }
}

void print_header5(FILE *fd, int qlib, struct db_str *qtt,
		   struct mngmsg *m_msp, struct pstruct *ppst,
		   int in_mem, long tot_memK) {

  /* for MX_MBLAST2, show some statistics results */
  if (m_msp->markx & MX_MBLAST2) {
    fprintf(fd,"\n\nLambda      K     H\n");
    fprintf(fd," %6.3f  %6.3f  %6.3f\n\n",ppst->pLambda,ppst->pK,ppst->pH);
    fprintf(fd,"\nGapped\nLambda\n");
    fprintf(fd," %6.3f  %6.3f  %6.3f\n",ppst->pLambda,ppst->pK,ppst->pH);
    fprintf(fd,"\nEffective search space used: %ld\n\n",m_msp->db.entries);
  }

  if (m_msp->markx & MX_M8COMMENT) {
    fprintf(fd, "# %s processed %d queries\n",prog_func,qlib);
  }

  if ( !((m_msp->markx & MX_M8OUT) || (m_msp->markx & MX_HTML))
       && (m_msp->markx & (MX_M10FORM+MX_M9SUMM))) {
    fprintf(fd,">>>///\n");
  }

  if ( m_msp->markx & MX_HTML) fputs("<pre>",fd); 
  if (m_msp->std_output) {
    print_sum(fd, qtt, &m_msp->db, in_mem, tot_memK);}
  if ( m_msp->markx & MX_HTML) fputs("</pre>\n",fd);
#ifdef HTML_HEAD
  if (m_msp->markx & MX_HTML) fprintf(fd,"</body>\n</html>\n");
#endif

  if (m_msp->markx & MX_MBLAST2) {
      fprintf(fd,"\n  Database: %s\n",m_msp->ltitle);
      fprintf(fd,"  Number of letters in database: %ld\n",m_msp->db.length);
      fprintf(fd,"  Number of sequences in database: %ld\n",m_msp->db.entries);
      fprintf(fd,"\n\n\nMatrix: %s\n",ppst->pam_name);
      fprintf(fd,"Gap Penalties: Existence: %d, Extension: %d\n",ppst->gdelval, ppst->ggapval);
  }
}

extern int fa_max_workers;

void
print_sum(FILE *fd, struct db_str *qtt, struct db_str *ntt, int in_mem, long tot_memK)
{
  double db_tt;
  char tstr1[26], tstr2[26];
  char memstr[256];

  strncpy(tstr1,ctime(&tdstart),sizeof(tstr1));
  strncpy(tstr2,ctime(&tddone),sizeof(tstr1));
  tstr1[24]=tstr2[24]='\0';

  /* Print timing to output file as well */

  fprintf(fd, "\n%ld residues in %ld query   sequences\n", qtt->length, qtt->entries);
  if (ntt->carry == 0) 
    fprintf(fd, "%ld residues in %ld library sequences\n", ntt->length, ntt->entries);
  else {
    db_tt = (double)ntt->carry*(double)LONG_MAX + (double)ntt->length;
    fprintf(fd, "%.0f residues in %ld library sequences\n", db_tt, ntt->entries);
  }

  memstr[0]='\0';
  if (tot_memK && in_mem != 0) {
    sprintf(memstr," in memory [%ldG]",(tot_memK >> 20));
  }

#if !defined(COMP_THR) && !defined(PCOMPLIB)
  fprintf(fd," Scomplib [%s%s]\n start: %s done: %s\n",verstr,mp_verstr,tstr1,tstr2);
#endif
#if defined(COMP_THR)
  fprintf(fd," Tcomplib [%s%s] (%d proc%s)\n start: %s done: %s\n", verstr, mp_verstr,
	  fa_max_workers, memstr, tstr1,tstr2);
#endif
#if defined(PCOMPLIB)
  fprintf(fd," Pcomplib [%s%s] (%d proc%s)\n start: %s done: %s\n", verstr, mp_verstr,
	  fa_max_workers, memstr, tstr1,tstr2);
#endif
#ifndef COMP_MLIB
  fprintf(fd," Scan time: ");
  ptime(fd, tscan - tprev);
  fprintf (fd," Display time: ");
  ptime (fd, tdone - tscan);
#else
  fprintf(fd," Total Scan time: ");
  ptime(fd, ttscan);
  fprintf (fd," Total Display time: ");
  ptime (fd, ttdisp);
#endif
  fprintf (fd,"\n");
  fprintf (fd, "\nFunction used was %s [%s%s]\n", prog_func,verstr,mp_verstr);
}

extern double zs_to_Ec(double zs, long entries);

#include "aln_structs.h"

void
prhist(FILE *fd, const struct mngmsg *m_msp,
       struct pstruct *ppst, 
       struct hist_str hist,
       int nstats, int sstats,
       struct db_str ntt,
       char *stat_info2,
       char *lib_range,
       char **info_gstring2,
       char **info_hstring,
       long tdiff)
{
  int i,j,hl,hll, el, ell, ev;
  char hline[80], pch, *bp;
  int mh1, mht;
  int maxval, maxvalt, dotsiz, ddotsiz,doinset;
  double cur_e, prev_e, f_int;
  double max_dev, x_tmp;
  double db_tt;
  int n_chi_sq, cum_hl=0, max_i=0, max_dev_i;
  double zs10_off;


  if (m_msp->markx & MX_HTML) fputs("<pre>\n",fd);
  else {fprintf(fd,"\n");}
  
  if (ppst->zsflag_f < 0) {
    if (!m_msp->nohist) {
      fprintf(fd, "  %7ld residues in %5ld sequences", ntt.length,ntt.entries);
      fprintf(fd, "%s\n",lib_range);
    }
    fprintf(fd,"Algorithm: %s\nParameters: %s\n",info_gstring2[0],info_gstring2[1]);
    return;
  }

  if (nstats > 20) { 
    zs10_off = ppst->zs_off * 10.0;

    max_dev = 0.0;
    mh1 = hist.maxh-1;			/* max value for histogram */
    mht = (3*hist.maxh-3)/4 - 1;	/* x-coordinate for expansion */
    n_chi_sq = 0;

    if (!m_msp->nohist && mh1 > 0) {
      for (i=0,maxval=0,maxvalt=0; i<hist.maxh; i++) {
	if (hist.hist_a[i] > maxval) maxval = hist.hist_a[i];
	if (i >= mht &&  hist.hist_a[i]>maxvalt) maxvalt = hist.hist_a[i];
      }
      cum_hl = -hist.hist_a[0];
      dotsiz = (maxval-1)/60+1;
      ddotsiz = (maxvalt-1)/50+1;
      doinset = (ddotsiz < dotsiz && dotsiz > 2);

      if (ppst->zsflag_f>=0)
	fprintf(fd,"       opt      E()\n");
      else 
	fprintf(fd,"     opt\n");

      prev_e =  zs_to_Ec((double)(hist.min_hist-hist.histint/2)-zs10_off,hist.entries);
      for (i=0; i<=mh1; i++) {
	pch = (i==mh1) ? '>' : ' ';
	pch = (i==0) ? '<' : pch;
	hll = hl = hist.hist_a[i];
	if (ppst->zsflag_f>=0) {
	  cum_hl += hl;
	  f_int = (double)(i*hist.histint+hist.min_hist)+(double)hist.histint/2.0;
	  cur_e = zs_to_Ec(f_int-zs10_off,hist.entries);
	  ev = el = ell = (int)(cur_e - prev_e + 0.5);
	  if (hl > 0  && i > 5 && i < (90-hist.min_hist)/hist.histint) {
	    x_tmp  = fabs(cum_hl - cur_e);
	    if ( x_tmp > max_dev) {
	      max_dev = x_tmp;
	      max_i = i;
	    }
	    n_chi_sq++;
	  }
	  if ((el=(el+dotsiz-1)/dotsiz) > 60) el = 60;
	  if ((ell=(ell+ddotsiz-1)/ddotsiz) > 40) ell = 40;
	  fprintf(fd,"%c%3d %5d %5d:",
		  pch,(i<mh1)?(i)*hist.histint+hist.min_hist :
		  mh1*hist.histint+hist.min_hist,hl,ev);
	}
	else fprintf(fd,"%c%3d %5d :",
		     pch,(i<mh1)?(i)*hist.histint+hist.min_hist :
		     mh1*hist.histint+hist.min_hist,hl);

	if ((hl=(hl+dotsiz-1)/dotsiz) > 60) hl = 60;
	if ((hll=(hll+ddotsiz-1)/ddotsiz) > 40) hll = 40;
	for (j=0; j<hl; j++) hline[j]='='; 
	if (ppst->zsflag_f>=0) {
	  if (el <= hl ) {
	    if (el > 0) hline[el-1]='*';
	    hline[hl]='\0';
	  }
	  else {
	    for (j = hl; j < el; j++) hline[j]=' ';
	    hline[el-1]='*';
	    hline[hl=el]='\0';
	  }
	}
	else hline[hl] = 0;
	if (i==1) {
	  for (j=hl; j<10; j++) hline[j]=' ';
	  sprintf(&hline[10]," one = represents %d library sequences",dotsiz);
	}
	if (doinset && i == mht-2) {
	  for (j = hl; j < 10; j++) hline[j]=' ';
	  sprintf(&hline[10]," inset = represents %d library sequences",ddotsiz);
	}
	if (i >= mht&& doinset ) {
	  for (j = hl; j < 10; j++) hline[j]=' ';
	  hline[10]=':';
	  for (j = 11; j<11+hll; j++) hline[j]='=';
	  hline[11+hll]='\0';
	  if (ppst->zsflag_f>=0) {
	    if (ell <= hll) hline[10+ell]='*';
	    else {
	      for (j = 11+hll; j < 10+ell; j++) hline[j]=' ';
	      hline[10+ell] = '*';
	      hline[11+ell] = '\0';
	    }
	  }
	}

	fprintf(fd,"%s\n",hline);
	prev_e = cur_e;
      }
    }
    max_dev_i = max_i*hist.histint+hist.min_hist;
  }
  else {
    max_dev = 0.0;
    n_chi_sq = 0;
    max_i = 0;
    max_dev_i = 0;
  }

  if (ppst->zsflag_f >=0 ) {
    if (!m_msp->nohist) {
      if (ntt.carry==0) {
	fprintf(fd, "  %7ld residues in %5ld sequences", ntt.length, ntt.entries);
      }
      else {
	db_tt = (double)ntt.carry*(double)LONG_MAX + (double)ntt.length;
	fprintf(fd, "  %.0f residues in %5ld library sequences", db_tt, ntt.entries);
      }
      fprintf(fd, "%s\n",lib_range);
    }
    fprintf(fd,"Statistics: %s\n",hist.stat_info);
    if (stat_info2) {
      fprintf(fd," Statistics E2: %s\n",stat_info2);
    }

#ifdef SAMP_STATS
    fprintf(fd," statistics sampled from %ld (%d) to %ld sequences\n",
	    (hist.entries > nstats ? nstats : hist.entries),sstats, hist.entries);
#else
    fprintf(fd," statistics extrapolated from %ld to %ld sequences\n",
	    (hist.entries > nstats ? nstats : hist.entries),hist.entries);
#endif

    if (!m_msp->nohist && cum_hl > 0) {
      fprintf(fd," Kolmogorov-Smirnov  statistic: %6.4f (N=%d) at %3d\n",
	      max_dev/(float)cum_hl, n_chi_sq,max_dev_i);
    }
    if (m_msp->markx & MX_M10FORM) {
      while ((bp=strchr(hist.stat_info,'\n'))!=NULL) *bp=' ';
      if (cum_hl <= 0) cum_hl = -1;
      sprintf(info_hstring[0],"; mp_extrap: %d %ld\n; mp_stats: %s\n; mp_KS: %6.4f (N=%d) at %3d\n",
	      MAX_STATS,hist.entries,hist.stat_info,max_dev/(float)cum_hl, 
	      n_chi_sq,max_dev_i);
    }
  }

  if (m_msp->markx & MX_M10FORM) {
    if ((bp = strchr(info_gstring2[1],'\n'))!=NULL) *bp = ' ';
    sprintf(info_hstring[1],"; mp_Algorithm: %s\n; mp_Parameters: %s\n",info_gstring2[0],info_gstring2[1]);
    if (bp != NULL ) *bp = '\n';
  }

  if (ppst->other_info != NULL) {
    fputs(ppst->other_info, fd);
  }

  fprintf(fd,"Algorithm: %s\nParameters: %s\n",info_gstring2[0],info_gstring2[1]);
  
  fprintf (fd," Scan time: ");
  ptime(fd, tdiff);
  fprintf(fd,"\n");
  if (m_msp->markx & MX_HTML) {
    fputs("</pre>\n<hr>\n",fd);
  }

  fflush(fd);
}

extern char prog_name[], *verstr;

#ifdef PCOMPLIB
#include "mpi.h"
#endif

void s_abort (char *p,  char *p1)
{
  int i;

  fprintf (stderr, "\n***[%s] %s%s***\n", prog_name, p, p1);
#ifdef PCOMPLIB
  MPI_Abort(MPI_COMM_WORLD,1);
  MPI_Finalize();
#endif
  exit (1);
}

void w_abort (char *p, char *p1)
{
  fprintf (stderr, "\n***[%s] %s%s***\n\n", prog_name, p, p1);
  exit (1);
}

extern struct a_res_str *
build_ares_code(unsigned char *aa0, int n0,
		unsigned char *aa1, struct seq_record *seq,
		int frame, int *have_ares, int repeat_thresh, 
		const struct mngmsg *m_msp, struct pstruct *ppst,
		void *f_str
		);

extern struct lmf_str *
re_openlib(struct lmf_str *, int outtty);

#define MAX_BLINE 2048
#define RANLIB (m_fptr->ranlib)

extern int
re_getlib(unsigned char *, unsigned char **, 
	  int, int, int, int, int, long *, long *, 
	  struct lmf_str *m_fptr);

/* 
   pre_load_best loads a set of sequences using re_getlib

   it should be used for getting sequences for shuffling, and for showbest() if m_msg->quiet

   it both opens the m_file_p buffer, gets the bline[] descriptions,
   and reads the actual sequences.  In reading the sequences, it
   should first allocate one large buffer so that individual buffers do not need to be freed.
*/

void
pre_load_best(unsigned char *aa1save, int maxn,
	      struct beststr **bbp_arr, int nbest,
	      struct mngmsg *m_msp)
{
  int i, n1, bl_len, tmp_bline_len, l_llen;
  int seq_buf_len;
  char bline[MAX_BLINE];
  unsigned char  *seq_buf_p;
  char *bline_buf_p;

  struct beststr *bbp;
  struct lmf_str *m_fptr;

  /* 
     calculate how much room we need for sequences and blines 
  */
  
  if (m_msp->pre_load_done) return;

  seq_buf_len = 1;
  for (i=0; i<nbest; i++) {
    /* we are not (currently) allocating more than n1+1, because alignment is not vectorized,
       if it were vectorized, we would need n+16
    */
#ifdef DEBUG
    if (bbp_arr[i]->n1 != bbp_arr[i]->seq->n1) {
      fprintf(stderr,"*** [%s/compacc2.c:933] n1 (%d) != seq->n1 (%d)\n",prog_name, bbp_arr[i]->n1, bbp_arr[i]->seq->n1);
    }
#endif

    seq_buf_len += bbp_arr[i]->seq->n1 + 1;
  }

  if ((m_msp->aa1save_buf_b=(unsigned char *)calloc(seq_buf_len, sizeof(char)))==NULL) {
    fprintf(stderr, "*** error - cannot allocate space[%d] for sequence encoding\n",seq_buf_len);
    exit(1);
  }
  else {
    seq_buf_p = m_msp->aa1save_buf_b+1;		/* ensure there is an initial '\0' */
  }
  
  l_llen = m_msp->aln.llen;
  if ((m_msp->markx & MX_M9SUMM) && m_msp->show_code != SHOW_CODE_ID) {
    l_llen += 40;
    if (l_llen > 200) l_llen=200;
  }

  tmp_bline_len = sizeof(bline)-1;
  if (!(m_msp->markx & MX_M10FORM) && !m_msp->long_info) {tmp_bline_len = l_llen-5;}

  /* allocate more bline than we need for simplicity */
  if ((bline_buf_p=m_msp->bline_buf_b=(char *)calloc(nbest*tmp_bline_len, sizeof(char)))==NULL) {
    fprintf(stderr, "*** error - cannot allocate space[%d] for bline descriptions\n",nbest*tmp_bline_len);
    exit(1);
  }

  for (i=0; i<nbest; i++) {
    bbp = bbp_arr[i];

    if ((m_fptr=re_openlib(bbp->mseq->m_file_p,!m_msp->quiet))==NULL) {
      fprintf(stderr,"*** cannot re-open %s\n",bbp->mseq->m_file_p->lb_name);
      exit(1);
    }
    RANLIB(bline,tmp_bline_len,bbp->mseq->lseek,bbp->mseq->libstr,m_fptr);
    bl_len = strlen(bline);
    bbp->mseq->bline = bline_buf_p;
    bbp->mseq->bline_max = m_msp->aln.llen;
    strncpy(bbp->mseq->bline, bline, bl_len);
    bline_buf_p += bl_len+1;


    /* make sure we get annotation if present, and sequence if necessary */
    if (bbp->seq->aa1b==NULL || (m_msp->ann_flg && bbp->seq->aa1_ann==NULL)) {
      n1 = re_getlib(aa1save, m_msp->ann_flg ? &(bbp->seq->aa1_ann) : NULL, 
		     maxn,m_msp->ldb_info.maxt3, m_msp->ldb_info.l_overlap,
		     bbp->mseq->cont,m_msp->ldb_info.term_code,
		     &bbp->seq->l_offset,&bbp->seq->l_off,bbp->mseq->m_file_p);
      if (n1 != bbp->seq->n1) {
	fprintf(stderr," *** error n1[%d/%d] != n1[%d] from re_getlib() at %s [maxn:%d/maxt3:%d]\n", 
		bbp->n1, bbp->seq->n1, n1, bbp->mseq->libstr, maxn, m_msp->ldb_info.maxt3);
      }

#ifdef DEBUG
      if (adler32(1L,aa1save,n1)!=bbp->adler32_crc) {
	fprintf(stderr," *** error adler32_crc from re_getlib()\n");
      }
#endif

      if (bbp->seq->aa1b == NULL)  {
	bbp->seq->aa1b = seq_buf_p;
	memcpy(bbp->seq->aa1b, aa1save, bbp->seq->n1+1);
	seq_buf_p += bbp->seq->n1+1;
      }
    }
  }
  m_msp->pre_load_done = 1;
}

/*  merge_ares_chains()

    seeks to merge two ares chains, producing a single chain that is
    sorted by sw_score.

    Strategy -- choose the chain with the highest score, and go down
    it until the head of the other chain has higher score, then link
    the other chain to the main chain, breaking the first, and
    continue the process.

    The two pointers, max_next and alt_next, keep track of the best
    and the alternate chain
 */


#undef SHOW_MERGE_CHAIN

struct a_res_str *
merge_ares_chains(struct a_res_str *cur_ares, 
		  struct a_res_str *tmp_ares, 
		  int score_ix,
		  const char *msg)
{
  struct a_res_str *max_next, *max_ares, *alt_ares, *prev_next;

  if (!tmp_ares) return cur_ares;

#ifdef SHOW_MERGE_CHAIN
  fprintf(stderr,"cur_ares->");
  for (max_next = cur_ares; max_next; max_next = max_next->next) {
    fprintf(stderr,"%d->",max_next->rst.score[score_ix]);
  }

  fprintf(stderr,"||\n");
  fprintf(stderr,"tmp_ares->");
  for (max_next = tmp_ares; max_next; max_next = max_next->next) {
    fprintf(stderr,"%d->",max_next->rst.score[score_ix]);
  }
  fprintf(stderr,"||\n");
#endif

  /* start with the maximum score */

  if (cur_ares->rst.score[score_ix] >= tmp_ares->rst.score[score_ix]) {
    max_ares = max_next = prev_next = cur_ares;
    alt_ares = tmp_ares;
  }
  else {
    max_ares = max_next = prev_next = tmp_ares;
    alt_ares = cur_ares;
  }

  while (max_next && alt_ares) {
    /* this is guaranteed true for the first iteration */
    while (max_next && max_next->rst.score[score_ix] >= alt_ares->rst.score[score_ix]) {
      prev_next = max_next;
      max_next = max_next->next;
    }
    if (max_next==NULL) break;
    else {	/* max_next->rst.score[score_ix] no longer greater, switch
		   pointers */
      prev_next->next = alt_ares;
      alt_ares = max_next;
      max_next = prev_next->next;
    }
  }

  /* we quit whenever max_next or alt_ares == NULL; if
     (max_next==NULL), then continue adding the rest of alt_ares */

  if (max_next==NULL) {
    prev_next->next = alt_ares;
  }


#ifdef SHOW_MERGE_CHAIN
  fprintf(stderr,"[%s] merge_ares->",msg);
  for (max_next = max_ares; max_next; max_next = max_next->next) {
    fprintf(stderr,"%d->",max_next->rst.score[score_ix]);
  }
  fprintf(stderr,"||\n\n");
#endif

  return max_ares;
}

/* copies from from to to shuffling */

extern int my_nrand(int, void *);

void
shuffle(unsigned char *from, unsigned char *to, int n, void *rand_state)
{
  int i,j; unsigned char tmp;

  if (from != to) memcpy((void *)to,(void *)from,n);

  for (i=n; i>0; i--) {
    j = my_nrand(i, rand_state);
    tmp = to[j];
    to[j] = to[i-1];
    to[i-1] = tmp;
  }
  to[n] = 0;
}

/* shuffles DNA sequences as codons */
void
shuffle3(unsigned char *from, unsigned char *to, int n, void *rand_state)
{
  int i, j, i3,j3; unsigned char tmp;
  int n3;

  if (from != to) memcpy((void *)to,(void *)from,n);

  n3 = n/3;

  for (i=n3; i>0; i--) {
    j = my_nrand(i, rand_state);
    i3 = i*3;
    j3 = j*3;
    tmp = to[j3];
    to[j3] = to[i3-1];
    to[i3-1] = tmp;
    tmp = to[j3+1];
    to[j3+1] = to[i3];
    to[i3] = tmp;
    tmp = to[j3+2];
    to[j3+2] = to[i3+1];
    to[i3+1] = tmp;
  }
  to[n] = 0;
}

/* "shuffles" by reversing the sequence */
void
rshuffle(unsigned char *from, unsigned char *to, int n)
{
  unsigned char *ptr = from + n;

  while (n-- > 0) {
    *to++ = *ptr--;
  }
  *to = '\0';
}

static int ieven = 0;
/* copies from from to from shuffling, ieven changed for threads */
void
wshuffle(unsigned char *from, unsigned char *to, int n, int wsiz, void *rand_state)
{
  int i,j, k, mm; 
  unsigned char tmp, *top;

  memcpy((void *)to,(void *)from,n);
	
  mm = n%wsiz;

  if (ieven) {
    for (k=0; k<(n-wsiz); k += wsiz) {
      top = &to[k];
      for (i=wsiz; i>0; i--) {
	j = my_nrand(i, rand_state);
	tmp = top[j];
	top[j] = top[i-1];
	top[i-1] = tmp;
      }
    }
    top = &to[n-mm];
    for (i=mm; i>0; i--) {
      j = my_nrand(i, rand_state);
      tmp = top[j];
      top[j] = top[i-1];
      top[i-1] = tmp;
    }
    ieven = 0;
  }
  else {
    for (k=n; k>=wsiz; k -= wsiz) {
      top = &to[k-wsiz];
      for (i=wsiz; i>0; i--) {
	j = my_nrand(i, rand_state);
	tmp = top[j];
	top[j] = top[i-1];
	top[i-1] = tmp;
      }
    }
    top = &to[0];
    for (i=mm; i>0; i--) {
      j = my_nrand(i, rand_state);
      tmp = top[j];
      top[j] = top[i-1];
      top[i-1] = tmp;
    }
    ieven = 1;
  }
  to[n] = 0;
}

int
sfn_cmp(int *q, int *s)
{
  if (*q == *s) return *q;
  while (*q && *s) {
    if (*q == *s) return *q;
    else if (*q < *s) q++;
    else if (*q > *s) s++;
  }
  return 0;
}

#ifndef MPI_SRC
#define ESS 49
#endif

void
revcomp(unsigned char *seq, int n, int *c_nt)
{
  unsigned char tmp;
  int i, ni;

  for (i=0, ni = n-1; i< n/2; i++,ni--) {
    tmp = c_nt[seq[i]];
    seq[i] = c_nt[seq[ni]];
    seq[ni] = tmp;
  }
  if ((n%2)==1) {
    i = n/2;
    seq[i] = c_nt[seq[i]];
  }
  seq[n]=0;
}

/* check to see whether this score (or a shuff score) should
   be included in statistics */
int samp_stats_idx (int *pre_nstats, int nstats, void *rand_state) {
  int jstats = -1;


  /* this code works when every score can be used for statistics
     estimates, but fails for fasta/[t]fast[xy] where only a fraction
     of scores are used */

  if (*pre_nstats < MAX_STATS) {
    jstats = (*pre_nstats)++;
  }

  /* here, the problem is that while we may have pre_nstats
     possible samplings, in some cases (-M subsets, fasta,
     [t]fast[xy] we don't have MAX_STATS samples yet.  Until we
     have MAX_STATS, we want more.  But the stats_idx strategy
     means that there may be additional samples in the buffers
     that are not reflected in nstats.
  */

  else {
#ifdef SAMP_STATS_LESS
    /* now we have MAX_STATS samples
       we want to sample 1/2 of 60K - 120K, 1/3 of 120K - 180K, etc */
    /* check every 15K to see if we have gone past the next threshold */

    /* pre_nstats cannot be incremented before the % to ensure
       that stats_inc is incremented exactly at 60000, 120000, etc.
       use ">=" in case increment comes later
       tests suggest the first 60K are sampled about 25% more
       than the rest
    */
    if (nstats < MAX_STATS) {
      jstats = MAX_STATS - my_nrand(MAX_STATS - nstats, rand_state)-1;
    }
    else if (((*pre_nstats)++ % (MAX_STATS/4)) == 0 && 
	     *pre_nstats >= stats_inc * MAX_STATS) {
      stats_inc = (*pre_nstats / MAX_STATS) + 1;
    }
    if ((*pre_nstats % stats_inc) == 0) {
      jstats = my_nrand(MAX_STATS, rand_state);
    }
#else
    /* this sampling strategy calls my_nrand() for every
       sequence > 60K, but provides a very uniform sampling */
    jstats = my_nrand(++(*pre_nstats), rand_state);
    if (jstats >= MAX_STATS) { jstats = -1;}
#endif
  }
  return jstats;
}

char *
build_link_data(char **link_lib_file_p, 
		struct mngmsg *m_msp, struct beststr **bestp_arr,
		int debug) {
  int i, status;
  char tmp_line[MAX_SSTR];
  char link_acc_file[MAX_STR];
  int link_acc_fd;
  char *link_lib_file;
  char *link_lib_str;
  char link_script[MAX_LSTR];
  int link_lib_indirect;
  int link_lib_type;
  char *bp, *link_bp;
  FILE *link_fd=NULL;		/* file for link accessions */

#ifndef UNIX
  return NULL;
#else
  /* get two tmpfiles, one for accessions, one for library */
  link_acc_file[0] = '\0';

  if ((link_lib_file=(char *)calloc(MAX_STR,sizeof(char)))==NULL) {
    fprintf(stderr,"[build_link_data] Cannot allocate link_lib_file");
  }
  link_lib_file[0] = '\0';

  if ((bp=getenv("TMP_DIR"))!=NULL) {
    strncpy(link_acc_file,bp,sizeof(link_acc_file));
    link_acc_file[sizeof(link_acc_file)-1] = '\0';
    SAFE_STRNCAT(link_acc_file,"/",sizeof(link_acc_file));
  }

  SAFE_STRNCAT(link_acc_file,"link_acc_XXXXXX",sizeof(link_acc_file));
  link_acc_fd = mkstemp(link_acc_file);
  strncpy(link_lib_file,link_acc_file,MAX_STR);
  link_acc_file[sizeof(link_acc_file)-1] = '\0';
  SAFE_STRNCAT(link_lib_file,".lib",MAX_STR);

  /* write out accessions to link_acc_file */
  if ((link_fd =fdopen(link_acc_fd,"w"))==NULL) {
    fprintf(stderr,"Cannot open link_acc_file: %s\n",link_acc_file);
    goto no_links;
  }

  for (i=0; i<m_msp->nskip + m_msp->nshow; i++) {
    if ((bp=strchr(bestp_arr[i]->mseq->bline,' '))!=NULL) {
      *bp = '\0';
    }
    fprintf(link_fd,"%s\t%.3g\n",bestp_arr[i]->mseq->bline,bestp_arr[i]->rst.escore);
    if (bp != NULL) *bp=' ';
  }
  fclose(link_fd);

  /* build link_script link_acc_file > link_lib_file */
  /* check for indirect */
  link_bp = &m_msp->link_lname[0];
  if (*link_bp == '@') {
    link_bp++;
  }
  /* remove library type */
  if ((bp=strchr(link_bp,' '))!=NULL) {
    *bp = '\0';
    sscanf(bp+1,"%d",&link_lib_type);
  }
  else {
    link_lib_type = 0;
  }

  strncpy(link_script,link_bp,sizeof(link_script));
  link_script[sizeof(link_script)-1] = '\0';
  SAFE_STRNCAT(link_script," ",sizeof(link_script));
  SAFE_STRNCAT(link_script,link_acc_file,sizeof(link_script));
  SAFE_STRNCAT(link_script," >",sizeof(link_script));
  SAFE_STRNCAT(link_script,link_lib_file,sizeof(link_script));

  /* un-edit m_msp->link_lname */
  if (bp != NULL) *bp = ' ';

  /* run link_script link_acc_file > link_lib_file */
  status = system(link_script);
  if (!debug) {
    unlink(link_acc_file);
  }
      
  if (status == NO_FILE_EXIT) {	/* my specific return for no links */
    goto no_links;
  }

  if (status < 0 || status == 127) {
    fprintf(stderr,"*** script: %s failed\n",link_script);
    goto no_links;
  }

  if ((link_fd=fopen(link_lib_file,"r"))==NULL) {
    goto no_links;
  }
  else fclose(link_fd);

  if ((link_lib_str=(char *)calloc(MAX_STR,sizeof(char)))==NULL) {
    fprintf(stderr,"[build_link_data] Cannot allocate link_lib_str");
  }

  /* build the file string (possibly @link_lib_file libtype) */
  link_lib_str[0]='\0';
  if (m_msp->link_lname[0] == '@') {
    SAFE_STRNCAT(link_lib_str,"@",MAX_STR);
  }
  SAFE_STRNCAT(link_lib_str,link_lib_file,MAX_STR);
  if (link_lib_type > 0) {
    sprintf(tmp_line," %d",link_lib_type);
    SAFE_STRNCAT(link_lib_str,tmp_line,MAX_STR);
  }

  *link_lib_file_p = link_lib_file;
  return link_lib_str;

 no_links:
  free(link_lib_file);
  *link_lib_file_p = NULL;
  return NULL;
#endif
}

/* save_best captures much of the complexity of saving the best scores
   and appropriately sampling the scores for statistical analysis. It
   does the following:

   (1) update s_info counts for functions like fasta/x/y that don't
       optimize every score 

   (2) for every result in the buffer:
       (a) decide if it should be used for statistical sampling
       (b) if the number of samples > MAX_STATS, then run
           process_hist() and update all the zscores
       (c) reset everything for next sequence

*/

#include "thr_buf_structs.h"
#ifndef PCOMPLIB
#define RESULTS_BUF reader_buf
#define XTERNAL
#include "thr_bufs2.h"
#else
#define RESULTS_BUF worker_buf
#include "pcomp_bufs.h"
#endif

extern char *prog_func;		/* function label */
extern int fa_max_workers;
extern struct buf_head *lib_buf2_list;
#ifdef DEBUG
void check_rbuf(struct buf_head *cur_buf);
#endif
extern void get_rbuf(struct buf_head **lib_buf, int max_work_buf);
extern void put_rbuf(struct buf_head *lib_buf, int max_work_buf);
extern void wait_rbuf(int max_work_buf);
extern void rbuf_done(int nthreads);
extern void put_rbuf_done(int nthreads, struct buf_head *lib_buf, 
			  int max_work_buf);
extern int
process_hist(struct stat_str *sptr, int nstats, 
	     const struct mngmsg *m_msg,
	     struct pstruct *ppst,
	     struct hist_str *hist, void **pstat_void, struct score_count_s *s_info, int do_hist);

extern void addhistz(double, struct hist_str *); /* scaleswn.c */
void selectbestz(struct beststr **, int, int );
extern double find_z(int score, double escore, int length, double comp,void *);
extern double zs_to_E(double zs,int n1, int dnaseq, long entries, struct db_str db);
extern struct beststr **bestp_arr;	/* array of pointers */
extern int nbest;
extern int nstats, nqstats, nrstats, pre_nstats, kstats, shuff_tot, sstats;
extern double zbestcut;	/* cut off for best z-score */
extern int bestfull;	/* index for selectbest() */
extern int stats_done;	/* flag for z-value processing */
extern void *rand_state;
extern struct stat_str *stats; /* array of scores for statistics from real
			   (or shuffled) sequences*/
extern struct stat_str *qstats;	/* array of scores for shuffled query stats */
extern struct stat_str *rstats;	/* array of scores from shuffled library */

/* in the current version (fasta_35_01) save_best is used by both
   threaded and unthreaded versions */

#define COPY_RST_P(d,s) 		\
{ d->rst.score[0] = s->rst.score[0];	\
  d->rst.score[1] = s->rst.score[1];	\
  d->rst.score[2] = s->rst.score[2];	\
  d->rst.valid_stat = s->rst.valid_stat; \
  d->rst.comp = s->rst.comp;		\
  d->rst.H = s->rst.H;			\
  d->rst.escore = s->rst.escore;	\
  d->rst.segnum = s->rst.segnum;	\
  d->rst.seglen = s->rst.seglen;	\
}

void
save_best(struct buf_head *lib_bhead_p, 
	  const struct mngmsg *m_msp, struct pstruct *ppst, 
	  struct db_str *ldb, FILE *fdata,
	  struct hist_str *histp, void **pstat_voidp,
	  struct score_count_s *s_info)
{
  double zscore;
  int i_score;
  struct beststr *bbp;
  struct buf2_data_s *rbuf_dp, *lib_buf2_dp;
  struct buf2_res_s *rbuf_rp, *lib_buf2_rp;
  int i, t_best, t_rbest, t_qrbest, tm_best, t_n1, sc_ix;
  int t_valid_stat, tr_valid_stat, use_shuff, zsflag_save;
  double e_score, tm_escore, t_rescore, t_qrescore;
  int buf2_cnt;

  if (!lib_bhead_p->hdr.have_results) return;
  if ((buf2_cnt = lib_bhead_p->hdr.buf2_cnt) <= 0) return;
  lib_buf2_dp = lib_bhead_p->buf2_data;
  lib_buf2_rp = lib_bhead_p->buf2_res;

  shuff_tot += lib_bhead_p->hdr.shuff_cnt;
  s_info->s_cnt[0] += lib_bhead_p->s_cnt_info.s_cnt[0];
  s_info->s_cnt[1] += lib_bhead_p->s_cnt_info.s_cnt[1];
  s_info->s_cnt[2] += lib_bhead_p->s_cnt_info.s_cnt[2];
  s_info->tot_scores += lib_bhead_p->s_cnt_info.tot_scores;;

  sc_ix = ppst->score_ix;

  t_best = t_rbest = t_qrbest = -BIGNUM;
  tm_escore = t_rescore = t_qrescore = FLT_MAX;
  t_valid_stat = tr_valid_stat = 0;
  if (ppst->zsflag >= 10 && ppst->zsflag < 20) { use_shuff = 1;}
  else { use_shuff = 0;}

#ifdef DEBUG
  if (fdata) {
    fprintf(fdata,">save_best: %d\n",buf2_cnt);
  }
#endif
  while (buf2_cnt--) { /* count down the number of results */
    rbuf_rp = lib_buf2_rp++;	/* step through the results buffer */
    rbuf_dp = lib_buf2_dp++;	/* step through the data buffer */

    /* perhaps should use explicit flag to indicate no score */
    if (rbuf_rp->rst.score[0] == -BIGNUM) continue;

    i_score = rbuf_rp->rst.score[sc_ix];
    e_score = rbuf_rp->rst.escore;

    zscore = (double)i_score;
    if (stats_done) {
      zscore=find_z(i_score, e_score, rbuf_dp->seq->n1,(double)rbuf_rp->rst.comp,
			  *pstat_voidp);
    }

    t_n1 = rbuf_dp->seq->n1;
    if (i_score > t_best) tm_best = t_best = i_score;
    if (e_score < tm_escore) tm_escore = e_score;
    if (rbuf_rp->rst.valid_stat > t_valid_stat) {
      t_valid_stat = 1;
    }

    if (m_msp->qshuffle) {
      if (rbuf_rp->qr_score > t_qrbest)
	t_qrbest = rbuf_rp->qr_score;
      if (rbuf_rp->qr_escore < t_qrescore)
	t_qrescore = rbuf_rp->qr_escore;
      
      if (rbuf_dp->frame == m_msp->nitt1 && nqstats < m_msp->shuff_max) {
	qstats[nqstats].n1 = rbuf_dp->seq->n1;	/* save the best score */
	qstats[nqstats].comp =  rbuf_rp->rst.comp;
	qstats[nqstats].H = rbuf_rp->rst.H;
	qstats[nqstats].escore = t_qrescore;
	qstats[nqstats++].score = t_qrbest;
	t_qrbest = -BIGNUM;	/* reset t_qrbest, t_qrescore */
	t_qrescore = FLT_MAX;
      }
    }	/* m_msp->qshuffle */

    if (use_shuff) {
      /* this check is required because some sequences scheduled to be
	 used for statistics may not in fact be returning a score (if
	 they are outside the -M range, for example.
       */
      if (rbuf_rp->r_rst.score[0] == -BIGNUM) { tr_valid_stat = 0; }
      if (rbuf_rp->r_rst.valid_stat > tr_valid_stat) {
	tr_valid_stat = 1;
      }
      if (rbuf_rp->r_rst.score[sc_ix] > t_rbest) {
	t_rbest = rbuf_rp->r_rst.score[sc_ix];
	t_rescore = rbuf_rp->r_rst.escore;
      }
    }

    /* need to look for frame 0 if TFASTA, then save stats at frame 6 */
    if (fdata) {
      fprintf(fdata,
	      "%-12s %6d %d %.5f %.5f %4d %4d %4d %2d %2d %4d %4d %4d %2d %2d %5d %8lld\n",
	      rbuf_dp->mseq->libstr, rbuf_dp->seq->n1,rbuf_dp->frame,rbuf_rp->rst.comp,rbuf_rp->rst.H,
	      rbuf_rp->rst.score[0],rbuf_rp->rst.score[1],rbuf_rp->rst.score[2],
	      t_valid_stat, rbuf_rp->rst.alg_info,
	      (rbuf_rp->r_rst.score[0]<0 ? -1 : rbuf_rp->r_rst.score[0]),
	      (rbuf_rp->r_rst.score[1]<0 ? -1 : rbuf_rp->r_rst.score[1]),
	      (rbuf_rp->r_rst.score[2]<0 ? -1 : rbuf_rp->r_rst.score[2]),
	      tr_valid_stat, rbuf_rp->r_rst.alg_info,
	      rbuf_dp->stats_idx, rbuf_dp->mseq->lseek);
    }

    /* statistics done for best score of set */

    if (rbuf_dp->frame == m_msp->nitt1) {
      ldb->entries++;
      ldb->length += t_n1;
      if (ldb->length > LONG_MAX) {
	ldb->length -= LONG_MAX; ldb->carry++;
      }
    }

    if (ppst->zsflag >= 0 && rbuf_dp->frame == m_msp->nitt1) {
      /* if this sample should be used for statistics */
      if (use_shuff) t_valid_stat = tr_valid_stat;
      if (t_valid_stat) { 
	/* we've got our initial MAX_STATS values */
	if (nstats >= MAX_STATS) {
	  if (!stats_done) {
	    zsflag_save = ppst->zsflag;
	    if (ppst->zsflag > 20) {
	      ppst->zsflag -= 20;
	    }
	    ppst->zsflag_f = process_hist(stats,nstats,m_msp, ppst,
					  histp, pstat_voidp,s_info, 0);
	    ppst->zsflag = zsflag_save;
	    kstats = nstats;
	    if (ppst->zsflag >= 0) {	/* this is redundant, but rare */
	      stats_done = 1;
	      for (i=0; i< nstats; i++) {
		bestp_arr[i]->zscore = 
		  find_z(bestp_arr[i]->rst.score[ppst->score_ix],
			 bestp_arr[i]->rst.escore, bestp_arr[i]->seq->n1,
			 bestp_arr[i]->rst.comp, *pstat_voidp);
	      }
	    }
	  }
	}
	else {
	  /* this logic allows stats_idx to be over-ruled for searches
	     where every query does not generate a score */
	  rbuf_dp->stats_idx = nstats;
	  nstats++;
	}
      }

      if (rbuf_dp->stats_idx >= 0 && t_valid_stat) {
	if (rbuf_dp->stats_idx >= MAX_STATS || nstats > MAX_STATS) {
	  fprintf(stderr, " *** error *** nstats index [%d] out of range [%d,%d]\n",
		  rbuf_dp->stats_idx, nstats,MAX_STATS);
	}
	else {  /* stats_idx is in range */
	  sstats++;
	  stats[rbuf_dp->stats_idx].n1 = t_n1;
	  stats[rbuf_dp->stats_idx].comp = rbuf_rp->rst.comp;
	  stats[rbuf_dp->stats_idx].H = rbuf_rp->rst.H;
	  if (use_shuff) { /* use shuffled score */
	    stats[rbuf_dp->stats_idx].escore  = t_rescore;
	    stats[rbuf_dp->stats_idx].score = t_rbest;
	  }
	  else { /* real score, not shuffled */
	    stats[rbuf_dp->stats_idx].escore  = tm_escore;
	    stats[rbuf_dp->stats_idx].score = tm_best;
	  } 
	} /* end stats_idx in range */
      }	/* end have valid stats_idx */

      if (t_valid_stat && stats_done && histp) {
	addhistz(find_z(t_best, tm_escore, rbuf_dp->seq->n1, (double) rbuf_rp->rst.comp,
			    *pstat_voidp), histp);
      }
      /* reset best scores */
      t_best = t_rbest = -BIGNUM;
      tm_escore = t_rescore = FLT_MAX;
      t_valid_stat = tr_valid_stat = 0;
    }

    /*
    if (rbuf_rp->rst.score[ppst->score_ix] > 200) {
      fprintf(stderr, "high score[%d]: %s %d: %d\n", rbuf_dp->seq->index,
	      rbuf_dp->mseq->libstr, rbuf_dp->seq->n1, rbuf_rp->rst.score[ppst->score_ix]);
    }
    */

    if (zscore > zbestcut) {
      if (nbest >= MAX_BEST) {
	bestfull = nbest-MAX_BEST/4;
	selectbestz(bestp_arr,bestfull-1,nbest);
	zbestcut = bestp_arr[bestfull-1]->zscore;
	nbest = bestfull;
      }
      bbp = bestp_arr[nbest++];

      COPY_RST_P(bbp, rbuf_rp);

      bbp->seq = rbuf_dp->seq;
      bbp->mseq = rbuf_dp->mseq;
      bbp->n1 = rbuf_dp->seq->n1;
#ifdef DEBUG
      bbp->adler32_crc = rbuf_dp->seq->adler32_crc;
#endif
      /* rbuf_dp->best_save is set after a rbuf_dp is entered into best_str */
      if (rbuf_dp->best_save) {	
	/* a previous rbuf_dp->seq is in best_str at best_save */
	if (rbuf_dp->best_save->seq == rbuf_dp->seq) {
	  /* the best_save->seq matches the rbuf_dp->seq */
	  bbp->bbp_link = rbuf_dp->best_save;
	  /* bbp_link tells where this ->seq can be found */
	}
	else {
	  bbp->bbp_link = NULL;
	}
      }
      rbuf_dp->best_save = bbp;
      lib_bhead_p->hdr.have_best_save = 1;
      bbp->zscore = zscore;
      bbp->frame = rbuf_dp->frame;
    }
  }
}

void
save_shuf(struct buf_head *lib_bhead_p, int nitt1, int shuff_max, int sc_ix,
	  struct score_count_s *s_info)
{
  struct buf2_data_s *rbuf_dp, *lib_buf2_dp;
  struct buf2_res_s *rbuf_rp, *lib_buf2_rp;
  int t_valid_stat;
  int t_rbest;
  double t_rescore;
  int buf2_cnt, jstats;
  static int kstats=0;


  lib_buf2_dp = lib_bhead_p->buf2_data;
  lib_buf2_rp = lib_bhead_p->buf2_res;
  buf2_cnt = lib_bhead_p->hdr.buf2_cnt;
  
  s_info->s_cnt[0] += lib_bhead_p->s_cnt_info.s_cnt[0];
  s_info->s_cnt[1] += lib_bhead_p->s_cnt_info.s_cnt[1];
  s_info->s_cnt[2] += lib_bhead_p->s_cnt_info.s_cnt[2];

  s_info->tot_scores += lib_bhead_p->s_cnt_info.tot_scores;
  /* this is done because we are not using r_rst->valid_stat to limit selection of scores */
  /*   s_info->s_cnt[sc_ix] = s_info->tot_scores; */

  t_rbest = -BIGNUM;
  t_valid_stat = 0;

  while (buf2_cnt--) { /* count down the number of results */
    rbuf_dp = lib_buf2_dp++;	/* step through the results buffer */
    rbuf_rp = lib_buf2_rp++;	/* step through the results buffer */

    /* perhaps should use explicit flag to indicate no score */
    if (rbuf_rp->r_rst.score[0] == -BIGNUM) continue;

    if (rbuf_rp->r_rst.score[sc_ix] > t_rbest) {
      t_rbest = rbuf_rp->r_rst.score[sc_ix];
      t_rescore = rbuf_rp->r_rst.escore;
    }

    if (rbuf_rp->r_rst.valid_stat > t_valid_stat) {
      t_valid_stat = 1;
    }

    /* statistics done for best score of set */
    /* currently no check for rst->valid_stat, which causes
       over-estimates of shuffles */

    if (rbuf_dp->frame == nitt1) {
      if (t_valid_stat) { 
	if (nrstats < shuff_max ) { kstats = jstats = nrstats++; }
	else {	/* randomly replace */
	  jstats = my_nrand(++kstats,rand_state);
	  if (jstats >= shuff_max) goto done;
	}

	rstats[jstats].n1 = rbuf_dp->seq->n1;
	rstats[jstats].comp = rbuf_rp->r_rst.comp;
	rstats[jstats].H = rbuf_rp->r_rst.H;
	rstats[jstats].escore  = t_rescore;
	rstats[jstats].score = t_rbest;
      done:
	t_rbest = -BIGNUM;
      }
    }
  }
}

int
save_align(struct buf_head *lib_bhead_p, struct beststr **bestp_arr) 
{
  struct buf2_ares_s *rbuf_ap, *lib_buf2_ap;
  int buf2_cnt;

  if (!lib_bhead_p->hdr.have_results || lib_bhead_p->hdr.buf2_cnt <= 0) return 0;

  lib_buf2_ap = lib_bhead_p->buf2_ares;
  buf2_cnt = lib_bhead_p->hdr.buf2_cnt;
  
  while (buf2_cnt-- > 0) { /* count down the number of results */
    rbuf_ap = lib_buf2_ap++;	/* step through the results buffer */
    bestp_arr[rbuf_ap->best_idx]->have_ares = rbuf_ap->have_ares;
    bestp_arr[rbuf_ap->best_idx]->a_res = rbuf_ap->a_res;
  }

  buf2_cnt = lib_bhead_p->hdr.buf2_cnt;
  lib_bhead_p->hdr.have_results = 0;
  lib_bhead_p->hdr.buf2_cnt = 0;
  return buf2_cnt;
}

void
buf_do_work(unsigned char **aa0,  int n0,
	    struct buf_head *lib_bhead_p, 
	    struct pstruct *ppst, void **f_str) {
  
  int buf2_cnt;
  unsigned long atmp;
  struct buf2_data_s *lib_buf2_dp;
  struct buf2_res_s *lib_buf2_rp;

  lib_bhead_p->s_cnt_info.s_cnt[0] = lib_bhead_p->s_cnt_info.s_cnt[1] =
    lib_bhead_p->s_cnt_info.s_cnt[2] = lib_bhead_p->s_cnt_info.tot_scores = 0;

  buf2_cnt = lib_bhead_p->hdr.buf2_cnt;
  lib_buf2_dp = lib_bhead_p->buf2_data;
  lib_buf2_rp = lib_bhead_p->buf2_res;

  while (buf2_cnt-- > 0) {

    lib_buf2_rp->rst.score[0] =
      lib_buf2_rp->rst.score[1] =
      lib_buf2_rp->rst.score[2] = -BIGNUM;

    if (lib_buf2_dp->seq->n1 < ppst->n1_low ||
	lib_buf2_dp->seq->n1 > ppst->n1_high ) {
      /* tells save_best() there is no stats score here -- not
	 necessary as -BIGNUM indicates no score */
      lib_buf2_dp->stats_idx = -1;
      goto next_seq;
    }

#ifdef DEBUG
    if (check_seq_range(lib_buf2_dp->seq->aa1b, lib_buf2_dp->seq->n1,
			(ppst->ext_sq_set ? ppst->nsqx: ppst->nsq), "buf_do_work()")) {
      fprintf(stderr, "[%s/buf_do_work] range error at: %d/%d (n1:%d)\n",
	      prog_func,lib_bhead_p->hdr.buf2_cnt - (buf2_cnt+1),
	      lib_bhead_p->hdr.buf2_cnt, lib_buf2_dp->seq->n1);
      goto next_seq;
    };

    /* also check for adler32_crc match */
    if (lib_buf2_dp->seq->adler32_crc != (atmp=adler32(1L,lib_buf2_dp->seq->aa1b, lib_buf2_dp->seq->n1))) {
      fprintf(stderr, "[%s/buf_do_work] CRC error [%lu!=%lu] at: %d/%d (n1:%d/l_offset:%ld)\n",
	      prog_func,lib_buf2_dp->seq->adler32_crc, atmp,
	      lib_bhead_p->hdr.buf2_cnt - (buf2_cnt+1),
	      lib_bhead_p->hdr.buf2_cnt,lib_buf2_dp->seq->n1,
	      lib_buf2_dp->seq->l_offset);
      goto next_seq;
    }
#endif

    do_work (aa0[lib_buf2_dp->frame], n0,
	     lib_buf2_dp->seq->aa1b, lib_buf2_dp->seq->n1,
	     lib_buf2_dp->frame, ppst, f_str[lib_buf2_dp->frame], 0, 0,
	     &(lib_buf2_rp->rst), &(lib_bhead_p->s_cnt_info));

  next_seq:
    lib_buf2_dp++;
    lib_buf2_rp++;
  }
  lib_bhead_p->hdr.have_results = 1;
}

void
buf_do_align(unsigned char **aa0,  int n0,
	     struct buf_head *lib_bhead_p, 
	     struct pstruct *ppst, const struct mngmsg *m_msp,
	     void **f_str) {

  int buf2_cnt, i, nsq;
  struct buf2_data_s *lib_buf2_dp;
  struct buf2_res_s *lib_buf2_rp;
  struct buf2_ares_s *lib_buf2_ap;
  struct rstruct rst;

  buf2_cnt = lib_bhead_p->hdr.buf2_cnt;
  lib_buf2_dp = lib_bhead_p->buf2_data;
  lib_buf2_rp = lib_bhead_p->buf2_res;
  lib_buf2_ap = lib_bhead_p->buf2_ares;

  while (buf2_cnt-- > 0) {
    if ( m_msp->stages > 1) { 
      /* this is not typically done unless m_msp->stages > 1 */
      do_opt (aa0[lib_buf2_dp->frame], n0, lib_buf2_dp->seq->aa1b, lib_buf2_dp->seq->n1,
	      lib_buf2_dp->frame, ppst, f_str[lib_buf2_dp->frame], &rst);
      lib_buf2_rp->rst.score[2]=rst.score[2];
    }

#ifdef DEBUG
    if (lib_buf2_dp->seq->aa1b == NULL) {
      fprintf(stderr,"[compacc2.c:buf_do_align() -- null aa1b\n");
      lib_buf2_ap->a_res = NULL;
      break;
    }
    if (check_seq_range(lib_buf2_dp->seq->aa1b, lib_buf2_dp->seq->n1,
			(ppst->ext_sq_set ? ppst->nsqx: ppst->nsq), "buf_do_align()")) {
      fprintf(stderr, "[%s/buf_do_align] range error at: %d/%d (n1:%d)\n",
	      prog_func,lib_bhead_p->hdr.buf2_cnt - (buf2_cnt+1),
	      lib_bhead_p->hdr.buf2_cnt, lib_buf2_dp->seq->n1);
    };

    /* also check for adler32_crc match */
    if (lib_buf2_dp->seq->adler32_crc != adler32(1L,lib_buf2_dp->seq->aa1b, lib_buf2_dp->seq->n1)) {
      fprintf(stderr, "[%s/buf_do_align] CRC error at: %d/%d (n1:%d)\n",
	      prog_func,lib_bhead_p->hdr.buf2_cnt - (buf2_cnt+1),
	      lib_bhead_p->hdr.buf2_cnt, lib_buf2_dp->seq->n1);
    }
#endif

    lib_buf2_ap->a_res = build_ares_code(aa0[lib_buf2_dp->frame], m_msp->n0, lib_buf2_dp->seq->aa1b, lib_buf2_dp->seq,
					lib_buf2_dp->frame, &lib_buf2_ap->have_ares,
					lib_buf2_dp->repeat_thresh, m_msp, ppst, f_str[lib_buf2_dp->frame] );

    lib_buf2_dp++;
    lib_buf2_ap++;
    lib_buf2_rp++;
  }
  lib_bhead_p->hdr.have_results = 1;
}

void
buf_qshuf_work(unsigned char *aa0s,  int n0,
	       struct buf_head *lib_bhead_p, 
	       struct pstruct *ppst, void *qf_str,
	       int ix_score)
{
  int buf2_cnt;
  struct buf2_data_s *lib_buf2_dp;
  struct buf2_res_s *lib_buf2_rp;
  struct rstruct rrst;

  struct score_count_s q_scnt_info;

  buf2_cnt = lib_bhead_p->hdr.buf2_cnt;
  lib_buf2_dp = lib_bhead_p->buf2_data;
  lib_buf2_rp = lib_bhead_p->buf2_res;

  while (buf2_cnt-- > 0) {
    rrst.score[0] = rrst.score[1] = rrst.score[2] = -BIGNUM;

    if (lib_buf2_dp->seq->n1 < ppst->n1_low ||
	lib_buf2_dp->seq->n1 > ppst->n1_high ) {
      lib_buf2_dp++;
      lib_buf2_rp++;
      continue;
    }

    do_work (aa0s, n0,
	     lib_buf2_dp->seq->aa1b, lib_buf2_dp->seq->n1,
	     lib_buf2_dp->frame, ppst, qf_str, 1, 0,
	     &rrst, &q_scnt_info);
    lib_buf2_rp->qr_score = rrst.score[ix_score];
    lib_buf2_rp->qr_escore = rrst.escore;

    lib_buf2_dp++;
    lib_buf2_rp++;
  }
}

/* buf_shuf_work -- called by work_thr() and unthreaded comp_lib/main()
   **aa0/n0 query sequence
   *aa1s buffer for shuffled sequence
   *lib_bhead_p pointer to data, results
   **f_str (function internal state save)
   ix_score - index of score to be used to choose better score for itt/nframe
   *rand_state  (random state save)
*/
void
buf_shuf_work(unsigned char **aa0,  int n0, unsigned char *aa1s, 
	      struct buf_head *lib_bhead_p, 
	      struct pstruct *ppst, void **f_str,
	      int ix_score, void *rand_state)
{
  int buf2_cnt;
  int shuff_cnt;
  struct buf2_data_s *lib_buf2_dp;
  struct buf2_res_s *lib_buf2_rp;
  struct rstruct rrst;

  lib_bhead_p->s_cnt_info.s_cnt[0] = lib_bhead_p->s_cnt_info.s_cnt[1] =
    lib_bhead_p->s_cnt_info.s_cnt[2] = lib_bhead_p->s_cnt_info.tot_scores = 0;

  shuff_cnt = 0;
  buf2_cnt = lib_bhead_p->hdr.buf2_cnt;
  lib_buf2_dp = lib_bhead_p->buf2_data;
  lib_buf2_rp = lib_bhead_p->buf2_res;

  while (buf2_cnt-- > 0) {
    lib_buf2_rp->r_rst.score[0] = lib_buf2_rp->r_rst.score[1] =
      lib_buf2_rp->r_rst.score[2] = -BIGNUM;
    lib_buf2_rp->r_rst.valid_stat = 0;

    shuff_cnt++;
    if (ppst->zs_win > 0) {
      wshuffle(lib_buf2_dp->seq->aa1b,aa1s,lib_buf2_dp->seq->n1,ppst->zs_win, rand_state);
    }
    else {
      if (ppst->shuffle_dna3) {shuffle3(lib_buf2_dp->seq->aa1b,aa1s,lib_buf2_dp->seq->n1, rand_state);}
      else {shuffle(lib_buf2_dp->seq->aa1b,aa1s,lib_buf2_dp->seq->n1, rand_state);}
    }

    /* rshuffle(lib_buf2_dp->seq->aa1b,aa1s,lib_buf2_dp->seq->n1); */

#ifdef DEBUG
    if (check_seq_range(aa1s, lib_buf2_dp->seq->n1,
			(ppst->ext_sq_set ? ppst->nsqx: ppst->nsq), "buf_do_align()")) {
      fprintf(stderr, "[%s/buf_shuff_work] range error at: %d/%d (n1:%d)\n",
	      prog_func,lib_bhead_p->hdr.buf2_cnt - (buf2_cnt+1),
	      lib_bhead_p->hdr.buf2_cnt, lib_buf2_dp->seq->n1);
    };
#endif

    do_work (aa0[lib_buf2_dp->frame], n0,
	     aa1s, lib_buf2_dp->seq->n1,
	     lib_buf2_dp->frame, ppst, f_str[lib_buf2_dp->frame], 0, 1,
	     &lib_buf2_rp->r_rst, &(lib_bhead_p->s_cnt_info));
    lib_buf2_dp++;
    lib_buf2_rp++;
  }
  lib_bhead_p->hdr.shuff_cnt = shuff_cnt;
  lib_bhead_p->hdr.have_results = 1;
}

/* buf_shuf_seq is designed to:
   (1) take a list of sequences (specified by bptr[])
   (2) collect them from the database if they are not already available
   (3) send them to the threads or shuffle them directly and calculate scores
*/
void
buf_shuf_seq(unsigned char **aa0, int n0,
	     unsigned char **aa1shuff_b, unsigned char *aa1save, int maxn,
	     struct beststr **bestp_arr, int nbest,
	     struct pstruct *ppst, struct mngmsg *m_msp,
	     struct mng_thr *m_bufi_p
#if !defined(COMP_THR) && !defined(PCOMPLIB)
	     , void **f_str
#endif
	     , struct score_count_s *s_info)
{
  unsigned char *aa1shuff;
  struct beststr *bbp, **tmp_bestp;
  char l_bline[MAX_SSTR];
  int n1lib_req, shuff_mult;
  long loffset, l_off;
  int n1, itt;
  int max_do_cnt, ndiff, prev_index;
  int istats;
  int i, j;

  /* these variables track buffers of library sequences */
  int cur_buf_size, max_buf_size;
  struct buf_head *lib_bhead_p;
  struct buf2_data_s *lib_buf2_dp;
  struct buf2_res_s *lib_buf2_rp;

/* (1) get the sequences into a buffer - the sequence information is
   currently in the bestp_arr - find out how many we have, and how
   many we will need - the number to shuffle */

/* figure out how much space we need */
/* this logic overestimates the space needed when both the forward and
   reverse strand of the same sequence is in bestp_arr[].  A solution
   is to sort the bestp_arr[] by the sequence index (or just ->seq),
   so that duplicates are adjacent.  We would then get a better nbest
   that was actually ndiff.  */

  if ((tmp_bestp = (struct beststr **)calloc(nbest, sizeof(struct beststr *)))==NULL) {
    fprintf(stderr," *** %s/buf_shuf_seq() *** cannot allocate tmp_bestp[%d]\n",prog_name, nbest);
    exit(1);
  }
  for (i = 0; i < nbest; i++) {
    tmp_bestp[i] = bestp_arr[i];
  }

  sortbesti(tmp_bestp, nbest);

  prev_index = -1;
  n1lib_req = ndiff = 0;
  for (i = 0; i < nbest; i++) {
    if (tmp_bestp[i]->seq->index > prev_index) {
      prev_index = tmp_bestp[i]->seq->index;
      n1lib_req += tmp_bestp[i]->n1+ 2;
      ndiff++;
    }
  }

#if !defined(COMP_THR) && !defined(PCOMPLIB)
  if (n1lib_req >= maxn) { /* we need new space, aa1shuff is too small */
    if ((*aa1shuff_b = aa1shuff = 
	 (unsigned char *)realloc(*aa1shuff_b, n1lib_req*sizeof(char)))==NULL) {
      fprintf(stderr," *** cannot realloc aa1shuff[%d]\n",n1lib_req);
      exit(1);
    }
  }
  else { aa1shuff = *aa1shuff_b;}
  *aa1shuff = '\0';
  aa1shuff++;

#else	/* threaded */
  if (n1lib_req < 2) {
    fprintf(stderr,"[%s/buf_shuf_seq] no residues to shuffle: %d (%d)\n", 
	    prog_func,n1lib_req,ndiff);
    exit(1);
  }

  if ((*aa1shuff_b = aa1shuff = 
       (unsigned char *)calloc(n1lib_req,sizeof(char)))==NULL) {
    fprintf(stderr," *** cannot calloc aa1shuff[%d]\n",n1lib_req);
    exit(1);
  }
  *aa1shuff = '\0';
  aa1shuff++;
#endif

  shuff_mult = (m_msp->shuff_max+1)/ndiff;
  istats = 0;
      
  /* setup lib_bhead buffers for shuffle comparisons */
#if defined(COMP_THR) || defined(PCOMPLIB)	/* threaded/parallel */
  /* max_do_cnt can be smaller than max_buf2_cnt, but not larger */
  max_do_cnt = min(m_bufi_p->max_buf2_res,
		   m_msp->shuff_max / (2 * fa_max_workers));
  /* we don't have a left over one, so we need one */
  get_rbuf(&lib_bhead_p,m_bufi_p->max_work_buf);
#else	/* not threaded */
  max_do_cnt = m_bufi_p->max_buf2_res;
  lib_bhead_p = lib_buf2_list;  /* equivalent to un-threaded get_rbuf() */
#endif
  max_buf_size = n1lib_req;
  cur_buf_size = 0;
  lib_bhead_p->hdr.buf2_cnt = 0;
  lib_bhead_p->hdr.have_results = 0;
  lib_bhead_p->hdr.stop_work = 0;
  lib_bhead_p->hdr.buf2_type=BUF2_DOSHUF;
  lib_bhead_p->hdr.seq_record_continuous = 0;
  lib_buf2_dp = lib_bhead_p->buf2_data;
  lib_buf2_rp = lib_bhead_p->buf2_res;

  /* read sequences into shuffle buffer */

  for (i = 0; i < ndiff; i++) {
    bbp = tmp_bestp[i];
    if (bbp->seq->aa1b == NULL) {
      /* get the sequence */
      (bbp->mseq->m_file_p->ranlib)(l_bline, sizeof(l_bline),
				    bbp->mseq->lseek,bbp->mseq->libstr,bbp->mseq->m_file_p);
      n1 = re_getlib(aa1save,NULL, maxn,m_msp->ldb_info.maxt3,
		     m_msp->ldb_info.l_overlap,bbp->mseq->cont,m_msp->ldb_info.term_code,
		     &loffset,&l_off,bbp->mseq->m_file_p);

      /* fprintf(stderr, " %d gets %d %d\n",i,tmp_bestp[i]->seq->n1,n1); */

      memcpy(aa1shuff, aa1save, n1+1);
      bbp->seq->aa1b = aa1shuff;
      aa1shuff += n1 + 1;
    }

    /* lib_buf2_dp is used up by scores, the sequence is not sent multiple times */
    cur_buf_size += bbp->seq->n1+1;
    for (j = 0; j < shuff_mult; j++ ) {
      for (itt = m_msp->revcomp; itt <= m_msp->nitt1; itt++) {
#ifdef PCOMPLIB
	lib_buf2_dp->seq_dup = 0;	/* mark first ->seq as original, not duplicate */
#endif
	lib_buf2_dp->seq = bbp->seq;
	/* this invalidates lib_buf2_p->seq */
	lib_buf2_dp->stats_idx = istats++;
	lib_buf2_dp->frame = itt;
	lib_buf2_dp++;		/* point to next buf2 */
	lib_buf2_rp++;		/* point to next buf2 */
	lib_bhead_p->hdr.buf2_cnt++;

	if (lib_bhead_p->hdr.buf2_cnt >= max_do_cnt ||
	    cur_buf_size > max_buf_size) {
	  /* (2) send sequences for shuffling */

#if defined(COMP_THR) || defined(PCOMPLIB)	/* threaded - fill and empty buffers */
	  /* provide empty buffer to workers */
	  lib_bhead_p->hdr.aa1b_used = cur_buf_size;
	  lib_bhead_p->hdr.have_data = 1;
	  lib_bhead_p->hdr.seq_record_continuous = 0;
	  put_rbuf(lib_bhead_p,m_bufi_p->max_work_buf);
	  get_rbuf(&lib_bhead_p,m_bufi_p->max_work_buf);
#else		/* non-thread - just do the searches */
	  if (lib_bhead_p->hdr.buf2_type & BUF2_DOSHUF) {
	    buf_shuf_work(aa0,m_msp->n0, aa1save, lib_bhead_p,
			  ppst, f_str, ppst->score_ix, rand_state);
	  }
#endif
	  /* (3) save results in the rstats structure */
	  if (lib_bhead_p->hdr.buf2_cnt > 0 && lib_bhead_p->hdr.have_results) {
	    save_shuf(lib_bhead_p,m_msp->nitt1,m_msp->shuff_max,ppst->score_ix,s_info);
	  }

	  lib_bhead_p->s_cnt_info.s_cnt[0] = lib_bhead_p->s_cnt_info.s_cnt[1] =
	    lib_bhead_p->s_cnt_info.s_cnt[2] = lib_bhead_p->s_cnt_info.tot_scores = 0;

	  lib_bhead_p->hdr.buf2_cnt = 0;
	  cur_buf_size = 0;
	  lib_bhead_p->hdr.have_results = 0;
	  lib_bhead_p->hdr.buf2_type=BUF2_DOSHUF;
	  lib_bhead_p->hdr.seq_record_continuous = 0; /* seq_records are coming from bestptr in any order */
	  lib_bhead_p->hdr.stop_work = 0;
	  lib_buf2_dp = lib_bhead_p->buf2_data;
	}
      } /* for (itt .. */
    }
  }			/* done with tmp_bestp[] */

  free(tmp_bestp);

#if defined(COMP_THR) || defined(PCOMPLIB)	/* if COMP_THR/PCOMPLIB - fill and empty buffers */
  /* check last buffers for any results */
  lib_bhead_p->hdr.seq_record_continuous = 0;
  put_rbuf(lib_bhead_p,m_bufi_p->max_work_buf);
    
  /* wait for the threads to finish */

  wait_rbuf(m_bufi_p->max_work_buf);
  /*
    fprintf(stderr, " num_reader[%d]-empty[%d]: %d\tnrstats: %d\n",
    num_reader_bufs,empty_reader_bufs,
    num_reader_bufs-empty_reader_bufs, nrstats);
  */

  for (i=0; i < num_reader_bufs; i++) {
    if (RESULTS_BUF[i]->hdr.buf2_cnt > 0 && RESULTS_BUF[i]->hdr.have_results) {
      save_shuf(RESULTS_BUF[i],m_msp->nitt1, m_msp->shuff_max, ppst->score_ix, s_info);
      RESULTS_BUF[i]->hdr.buf2_cnt = RESULTS_BUF[i]->hdr.have_results = 0;
    }
  }
#else	/* just do the searches */
  /* aa1save is used for shuffles, not aa1shuf, because aa1shuf
     has library sequences */
  buf_shuf_work(aa0,m_msp->n0, aa1save, lib_bhead_p,
		ppst, f_str, ppst->score_ix, rand_state);

  save_shuf(lib_bhead_p,m_msp->nitt1,m_msp->shuff_max, ppst->score_ix, s_info);
  lib_bhead_p->hdr.buf2_cnt = lib_bhead_p->hdr.have_results = 0;
#endif
}

/* buf_align_seq is structurally almost identical to buf_shuf_seq,
   except that the appropriate sequences are pre-loaded into bbp->seq
   (and ->bline), and it gets bbp->a_res, rather than scores */

void
buf_align_seq(unsigned char **aa0, int n0,
	      struct beststr **bestp_arr, int nbest,
	      struct pstruct *ppst, struct mngmsg *m_msp,
	      struct mng_thr *m_bufi_p
#if !defined(COMP_THR) && !defined(PCOMPLIB)
	      , void **f_str
#endif
	     )
{
  struct beststr *bbp;
  int max_align_cnt;
  int i, n_pre_align;
  int cur_buf_size, max_buf_size;
  struct buf_head *lib_bhead_p;
  struct buf2_data_s *lib_buf2_dp;
  struct buf2_ares_s *lib_buf2_ap;
      
  /* setup lib_bhead buffers for alignments */
#if defined(COMP_THR) || defined(PCOMPLIB)	/* threaded */
  /* max_do_cnt can be smaller than max_buf2_res, but not larger */
#ifdef COMP_THR
  max_align_cnt = min(m_bufi_p->max_buf2_res,
		      nbest / (4 * fa_max_workers));
#else
  max_align_cnt = min(m_bufi_p->max_buf2_res, nbest / fa_max_workers);
#endif
  if (max_align_cnt < 1) max_align_cnt = 1;

  get_rbuf(&lib_bhead_p,m_bufi_p->max_work_buf);
#else	/* not threaded */
  max_align_cnt = m_bufi_p->max_buf2_res;
  lib_bhead_p = lib_buf2_list;  /* equivalent to un-threaded get_rbuf() */
#endif

  max_buf_size = lib_bhead_p->hdr.aa1b_size;
  lib_bhead_p->hdr.buf2_cnt = 0;
  lib_bhead_p->hdr.have_results = 0;
  lib_bhead_p->hdr.stop_work = 0;
  lib_bhead_p->hdr.buf2_type=BUF2_DOALIGN;
  lib_buf2_dp = lib_bhead_p->buf2_data;
  lib_buf2_ap = lib_bhead_p->buf2_ares;

  /* read sequences into align buffer */

  n_pre_align = 0;
  cur_buf_size = 0;
  for (i = 0; i < nbest; i++) {
    bbp = bestp_arr[i];

    /* this invalidates lib_buf2_p->seq */
    lib_buf2_dp->seq = bbp->seq;
    cur_buf_size += bbp->seq->n1+1;
    lib_buf2_dp->frame = bbp->frame;
    lib_buf2_dp->repeat_thresh = bbp->repeat_thresh;
#ifdef PCOMPLIB
    lib_buf2_dp->seq_dup = 0;
#endif
    lib_buf2_ap->have_ares = 0;
    lib_buf2_ap->a_res = NULL;
    lib_buf2_ap->best_idx = i;
    lib_buf2_dp++;		/* point to next buf2_data */
    lib_buf2_ap++;		/* point to next buf2_ares */
    lib_bhead_p->hdr.buf2_cnt++;

    if (lib_bhead_p->hdr.buf2_cnt >= max_align_cnt ||
	cur_buf_size >= max_buf_size - m_msp->ldb_info.maxn) {
/* (2) send sequences for alignment */
#if defined(COMP_THR) || defined(PCOMPLIB)	/* threaded - fill and empty buffers */
      /* provide empty buffer to workers */
      lib_bhead_p->hdr.seqr_cnt = lib_bhead_p->hdr.buf2_cnt;	/* for alignments, they are the same */
      lib_bhead_p->hdr.have_data = 1;
      lib_bhead_p->hdr.aa1b_used = cur_buf_size;
      lib_bhead_p->hdr.seq_record_continuous = 0;
      put_rbuf(lib_bhead_p,m_bufi_p->max_work_buf);
      get_rbuf(&lib_bhead_p,m_bufi_p->max_work_buf);
#else		/* non-thread - just do the searches */
      buf_do_align(aa0, m_msp->n0, lib_bhead_p, ppst, m_msp, f_str); 
#endif

/* (3) save alignments */
      if (lib_bhead_p->hdr.buf2_cnt > 0 && lib_bhead_p->hdr.have_results) {
	n_pre_align += save_align(lib_bhead_p,bestp_arr);
      }

      cur_buf_size = 0;
      max_buf_size = lib_bhead_p->hdr.aa1b_size;
      lib_bhead_p->hdr.buf2_cnt = 0;
      lib_bhead_p->hdr.have_results = 0;
      lib_bhead_p->hdr.buf2_type=BUF2_DOALIGN;
      lib_bhead_p->hdr.stop_work = 0;
      lib_buf2_dp = lib_bhead_p->buf2_data;
      lib_buf2_ap = lib_bhead_p->buf2_ares;
    }
  }			/* done with bestp_arr[] */

#if defined(COMP_THR) || defined(PCOMPLIB)	/* if COMP_THR - fill and empty buffers */
  /* check last buffers for any results */
  lib_bhead_p->hdr.seqr_cnt = lib_bhead_p->hdr.buf2_cnt;	/* for alignments, they are the same */
  lib_bhead_p->hdr.have_data = 1;
  lib_bhead_p->hdr.aa1b_used = cur_buf_size;
  lib_bhead_p->hdr.seq_record_continuous = 0;
  put_rbuf(lib_bhead_p,m_bufi_p->max_work_buf);
    
  /* wait for the threads to finish */

  wait_rbuf(m_bufi_p->max_work_buf);

  for (i=0; i < num_reader_bufs; i++) {
    if (RESULTS_BUF[i]->hdr.buf2_cnt > 0 && RESULTS_BUF[i]->hdr.have_results) {
      n_pre_align += save_align(RESULTS_BUF[i],bestp_arr);
      RESULTS_BUF[i]->hdr.buf2_cnt = RESULTS_BUF[i]->hdr.have_results = 0;
    }
  }
#else	/* just do the searches */
  buf_do_align(aa0, m_msp->n0, lib_bhead_p, ppst, m_msp, f_str); 
  n_pre_align += save_align(lib_bhead_p,bestp_arr);
  lib_bhead_p->hdr.buf2_cnt = lib_bhead_p->hdr.have_results = 0;
#endif

  if (n_pre_align != nbest) {
    fprintf(stderr,"*** error n_pre_align:%d != nbest: %d\n",n_pre_align, nbest);
  }
  for (i=0; i < nbest; i++) {
    if (bestp_arr[i]->a_res == NULL) {
      fprintf(stderr, "*** error - have NULL a_res: %d\n",i);
    }
  }
}

int
check_seq_range(unsigned char *aa1b, int n1, int nsq, char *str) {
  int i, range_error;
  unsigned char *aa1p;

  range_error = 0;
  for (aa1p = aa1b, i=0; i < n1; i++, aa1p++) {
    if (*aa1p > nsq) {
      range_error = 1;
      /*      fprintf(stderr, "%s seq %d (%c) out of range at %d\n",
	      str, *aa1p, *aa1p,i);
      */
      *aa1p = 1;
    }
  }
  return range_error;
}

