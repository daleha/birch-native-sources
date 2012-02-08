
/* copyright (c) 2010 William R. Pearson and the U. of Virginia */

/* $Id: build_ares.c 625 2011-03-23 17:21:38Z wrp $ */
/* $Revision: 625 $ */

/* build_ares_code is called by showbest() (in threaded/serial code) or
   p2_workcomp in PCOMPLIB code to produce the cur_ares-> chain that
   is displayed in showbest().

   For PCOMPLIB, the cur_ares->chain is passed to bbp-a_res by
   do_stage2(), where it is available to showbest();

   By using this code, the a_res chain used in either mode will be the
   same, so the code required to display an a_res should be the same.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "structs.h"
#include "param.h"

/* #include "mm_file.h" */
#include "best_stats.h"

#include "drop_func.h"

extern void calc_coord(int n0, int n1, long qoffset, long loffset,
		      struct a_struct *aln);

extern void calc_astruct(struct a_struct *aln_p, struct a_res_str *a_res_p);

/* in build_ares_code, *aa1 is separate from *seq because *seq has
   permanent information about aa1, but aa1 may be temporary */

struct a_res_str *
build_ares_code(unsigned char *aa0, int n0, 
		unsigned char *aa1,
		struct seq_record *seq,
		int frame, int *have_ares, int repeat_thresh,
		const struct mngmsg *m_msp, struct pstruct *ppst,
		void *f_str
		)
{
  unsigned char *aa1_ann;
  int n1;
  struct rstruct rst;
  struct a_res_str *my_ares_p, *cur_ares_p;
  struct a_struct *aln_p;
  long loffset;		/* loffset is offset from beginning of real sequence */
  long l_off;		/* l_off is the the virtual coordinate of residue 1 */
  int seqc_max, annc_max;
  char *seq_code, *ann_code;
  int seq_code_len, ann_code_len;

  n1 = seq->n1;
  aa1_ann = seq->aa1_ann;
  loffset = seq->l_offset;
  l_off = seq->l_off;

  if (! (*have_ares & 0x1)) {	/* we don't have an a_res, and we need one */

    my_ares_p = do_walign(aa0, n0, aa1, n1, frame, 
			  repeat_thresh, ppst, f_str,
			  have_ares);
  }
  else {	/* we already have the a_res */
    pre_cons(aa1,n1,frame,f_str);
  }

  /* here, we need to loop through all the alignments, and produce
     the statistics/codes for each */

  for (cur_ares_p = my_ares_p; cur_ares_p != NULL; cur_ares_p = cur_ares_p->next) {

    seqc_max = my_ares_p->nres + 4*m_msp->aln.llen+4;
    cur_ares_p->aln_code = seq_code = NULL;
    cur_ares_p->aln_code_n = seq_code_len = 0;
    cur_ares_p->ann_code = NULL;
    cur_ares_p->ann_code_n = 0;

    aln_p = &cur_ares_p->aln;

    /* this sets a number of constants, from the alignment function
       and frame, and only needs to be called once */
    aln_func_vals(frame, aln_p);

    if ((m_msp->show_code & SHOW_CODE_ALIGN) == SHOW_CODE_ALIGN) {
      cur_ares_p->aln_code = seq_code=(char *)calloc(seqc_max,sizeof(char));
      /* if we have an annotation string, allocate space for the
	 encoded annotation */
      if (m_msp->ann_arr[0] != '\0') {
	/* the annotation encoding can be considerably longer than
	   the alignment encoding */
	annc_max = 4*seqc_max;
	cur_ares_p->ann_code = ann_code=(char *)calloc(annc_max,sizeof(char));
      }
      else {
	ann_code = NULL;
	annc_max = 0;
      }

      if (seq_code != NULL) {
	calc_astruct(aln_p, cur_ares_p);

	/* we need this for offset information for calc_code, but it is
	   incomplete so we must do it again */
	calc_coord(m_msp->n0,n1,
		   m_msp->q_offset + (m_msp->q_off-1) + (m_msp->sq0off-1),
		   loffset + (l_off-1) + (m_msp->sq1off-1),
		   aln_p);

	aln_p->lc=calc_code(aa0, m_msp->n0,
			    aa1,n1, 
			    aln_p,cur_ares_p,
			    ppst,seq_code,seqc_max,
			    m_msp->ann_arr,
			    m_msp->aa0a, aa1_ann,
			    ann_code, annc_max,
			    f_str, m_msp->show_code);

	cur_ares_p->aln_code_n = seq_code_len = strlen(seq_code);
	if (seq_code[1] == '0' && seq_code[0] == '=') {
	  fprintf(stderr," code begins with 0: %s\n", seq_code);
	}

	if (ann_code != NULL) ann_code_len = strlen(ann_code);
	else ann_code_len = 0;
	cur_ares_p->ann_code_n = ann_code_len;
      }
    }
    else {
      aln_p->lc=calc_id(aa0,m_msp->n0,aa1,n1,
			aln_p, cur_ares_p,
			ppst,f_str);
    }
    /* this should be all the information we need on the alignment */
  } /* end for (cur_ares_p;) */
  return my_ares_p;
}
