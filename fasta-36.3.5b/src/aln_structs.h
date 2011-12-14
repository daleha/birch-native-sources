/* $Id: aln_structs.h 685 2011-04-18 18:56:17Z wrp $ */
/* $Revision: 685 $  */

#ifndef A_STRUCT
#define A_STRUCT

struct a_struct {
  int smin0;		/* coordinate of display start in seqc0 */
  int smin1;		/* coordinate of display start in seqc1 */
  int amin0, amax0;	/* coordinate of alignment start in seqc0 */
  int amin1, amax1;	/* coordinate of alignment start in seqc1 */
  int calc_last_set; 	/* boolean that indicates structure was set by
			   calc_code, calc_cons, etc */

  int llen;
  int llcntx, llcntx_set, showall;

  int qlrev, qlfact;
  int llrev, llfact, llmult;
  int frame;

  int nident, nsim, npos, nmismatch, lc, ngap_q, ngap_l, nfs;	/* number of identities, gaps in q, l */
  long q_start_off, q_end_off;	/* used in -m 9 calculate for full segment offset */
  long l_start_off, l_end_off;	
  long q_offset, l_offset;	/* offsets that include everything */
  long d_start0,d_stop0;
  long d_start1,d_stop1;
};

struct a_res_str {
  int sw_score;		/* do_walign() score */
  struct rstruct rst;
  int min0, max0;	/* boundaries of alignment in aa0 */
  int min1, max1;	/* boundaries of alignment in aa1 */
  int v_start, v_len;	/* virtual start, length */
  int *res;		/* encoded alignment */
  int nres;		/* length of decoded alignment */
  int mres;		/* length of encoding in res[] */
  int n1;		/* length of library sequence used for this (sub) alignment */
  struct a_res_str *next;	/* pointer to next alignment */

  /* encoded alignment/annotation information */
  char *aln_code;
  int aln_code_n;
  char *ann_code;
  int ann_code_n;
  struct a_struct aln;
};
#endif
