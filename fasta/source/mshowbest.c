
/* copyright (c) 1996, 1997, 1998, 1999 William R. Pearson and the
   U. of Virginia */

/* $Id: mshowbest.c 846 2011-10-18 18:50:49Z wrp $ */
/* $Revision: 846 $  */

/* 18-Oct-2011 -- remove calculation of m_msp->nshow from showbest();
   m_msp->nshow is set by caller */

/* 2-April-2009 changes to simplify interactive display logic.  Coming
   into showbest(), things are interactive (quiet==0) or use
   m_msg.nshow */

/*   29-Oct-2003 - changes so that bbp->mseq->cont < 0 => aa1 sequence is
     already in aa1, no re_openlib or re_getlib required
*/

/*   14-May-2003 Changes to use a more consistent coordinate numbering
     system for displays.  aln->d_start[01] is now consistently used
     to report the start of the alignment in all functions, and
     mshowbest.c has been modified to use d_start[01] instead of
     d_start[01]-1.  aln->min[01] now starts at 0 for all functions;
     instead of 1 for some functions (dropnfa.c, dropgsw.c, dropfs2.c
     earlier).
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "structs.h"
#include "param.h"
#include "mm_file.h"
#include "best_stats.h"

#define MAX_BLINE 256

/* function calls necessary to re_getlib() the sequence and, do
   alignments, if necessary
*/

#define RANLIB (m_fptr->ranlib)


int
re_getlib(unsigned char *, unsigned char **, 
	  int, int, int, int, int, long *, long *, 
	  struct lmf_str *m_fptr);

#include "drop_func.h"

struct a_res_str *
build_ares_code(unsigned char *aa0, int n0,
		unsigned char *aa1, struct seq_record *seq,
		int frame, int *have_ares, int repeat_thresh, 
		struct mngmsg *m_msp, struct pstruct *ppst,
		void *f_str
		);

struct lmf_str *re_openlib(struct lmf_str *, int outtty);

extern void calc_coord(int n0, int n1, long qoffset, long loffset,
		      struct a_struct *aln);

extern void calc_astruct(struct a_struct *aln_p, struct a_res_str *a_res_p);
extern double find_z(int score, double escore, int length, double comp,void *);
extern double zs_to_E(double zs, int n1, int dnaseq, long db_size, struct db_str db);
extern double zs_to_bit(double, int, int);
extern int E1_to_s(double e_val, int n0, int n1, int db_size, void *pu);

void header_aux(FILE *);
void show_aux(FILE *, struct beststr *);
void w_abort (char *p, char *p1);

extern double zs_to_bit(double, int, int);

/* showbest() shows a list of high scoring sequence descriptions, and
   their rst.scores.  If -m 9, then an additional complete set of
   alignment information is provided.

   If PCOMPLIB or m_msg.quiet then the number of high scores to be
   shown is pre-determined by m_msg.mshow before showbest is called.

   The comp_lib2.c version re_getlib()'s the sequence for its
   discription, and then does another alignment for -m 9 (Thus, it
   needs an f_str.  The PCOMPLIB version has everything available in
   beststr before showbest() is called.
*/

void showbest (FILE *fp, unsigned char **aa0, unsigned char *aa1save, int maxn,
	       struct beststr **bptr,int nbest, int qlib, 
	       struct mngmsg *m_msp,
	       struct pstruct *ppst, struct db_str db,
	       char **info_gstring2,
	       void **f_str
)
{
  unsigned char *aa1;
  int ntmp = 0;
  char bline[MAX_BLINE], fmt[40], pad[MAX_BLINE], fmt2[40], rline[40];
  char l_name[128];
  int istart = 0, istop, ib;
  int nshow;		/* number of sequences shown before prompt,
			   and ultimately displayed */
  int first_line;
  int quiet;
  int r_margin;
  struct beststr *bbp;
  int n1tot;
  char *bp, *bline_p;
  char rel_label[12];
  char tmp_str[20], *seq_code, *ann_code;
  int seq_code_len, ann_code_len;
  long loffset;		/* loffset is offset from beginning of real sequence */
  long l_off;		/* l_off is the the virtual coordinate of residue 1 */
  int n1, ranlib_done;
  struct rstruct rst;
  int l_score0, ngap;
  double lzscore, lzscore2, lbits;
  float percent, gpercent;
  struct a_struct *aln_p;
  struct a_res_str *cur_ares_p;
  struct rstruct *rst_p;
  int gi_num;
  char html_pre_E[120], html_post_E[120];

  struct lmf_str *m_fptr;

  rel_label[0]='\0';

  quiet = m_msp->quiet;

  if (m_msp->aln.llen > MAX_BLINE) m_msp->aln.llen = MAX_BLINE;

  if (ppst->zsflag < 0) r_margin = 10;
  else if (ppst->zsflag>=0  && m_msp->srelv > 1 ) r_margin = 19;
  else r_margin = 10;

  if (m_msp->markx & MX_M9SUMM && m_msp->show_code == SHOW_CODE_ID) {
#ifdef SHOWSIM
    r_margin += 15;
#else
    r_margin += 10;
#endif
  }
  else if (m_msp->markx & MX_MBLAST2) {
    r_margin -= 10;
  }

  if (m_msp->markx & MX_HTML) {
    strncpy(html_pre_E,"<font color=\"darkred\">",sizeof(html_pre_E));
    strncpy(html_post_E,"</font>",sizeof(html_post_E));

  }
  else {
    html_pre_E[0] = html_post_E[0] = '\0';
  }

  if (m_msp->nframe < 0) {
    sprintf(fmt,"%%-%ds (%%4d)",m_msp->aln.llen-r_margin);
  }
  else {
    sprintf(fmt,"%%-%ds (%%4d)",m_msp->aln.llen-(r_margin+4));
  }
  sprintf(fmt2,"%%-%ds",m_msp->aln.llen-r_margin+8);

  memset(pad,' ',m_msp->aln.llen-(r_margin+6));
  pad[m_msp->aln.llen-(r_margin+12)]='\0';

  nshow = min(m_msp->nshow,nbest);

  if ((bp = strchr (m_msp->qtitle, '\n')) != NULL) *bp = '\0';
  if (m_msp->markx & MX_M8OUT) {
    if ((bp = strchr (m_msp->qtitle, ' ')) != NULL) *bp = '\0';
  }

/*   fprintf (fp, "%3d %s\n", qlib,m_msp->qtitle); */

  if (m_msp->markx & MX_HTML) fprintf(fp,"<pre>");

  /* **************************************************************** */
  /* done with display format */
  /* **************************************************************** */

  /* **************************************************************** */
  /* prompt for number of best scores if quiet == 0 */
  /* **************************************************************** */

  if (quiet == 0) {	/* interactive */
    nshow = min(m_msp->nshow, nbest);
    printf(" How many scores would you like to see? [%d] ",nshow);
    fflush(stdout);
    if (fgets(rline,20,stdin)==NULL) exit(0);
    if (rline[0]!='\n' && rline[0]!=0) sscanf(rline,"%d",&nshow);
    if (nshow > nbest) nshow=nbest;
    if (nshow<=0) nshow = min(20,nbest);
  }

  /* display number of hits for -m 8C (Blast Tab-commented format) */
  if (m_msp->markx & MX_M8COMMENT) {
    /* line below copied from BLAST+ output */
    fprintf(fp,"# Fields: query id, subject id, %% identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score\n");
    fprintf(fp,"# %d hits found\n",nshow);
  }

  /* **************************************************************** */
  /* have number of scores in interactive or quiet mode */
  /* display "The best scores are" */
  /* **************************************************************** */

  if (m_msp->markx & MX_MBLAST2) {
    fprintf(fp, "%81s\n"," Score     E");
    fprintf(fp, "Sequences producing significant alignments:                          (Bits)  Value\n\n");
  }
  else if (!(m_msp->markx & MX_M8OUT)) {
    if (ppst->zsflag >= 0) {
      if (m_msp->z_bits==1) {/* show bit score */
	fprintf(fp,"\nThe best%s scores are:%s%s bits %sE(%ld)%s",
		rel_label,pad,m_msp->label,html_pre_E,ppst->zdb_size,html_post_E);
	if (ppst->zsflag > 20) {
	  fprintf(fp," E2()");
	}
      }
      else {/* show z-score */
	fprintf(fp,"\nThe best%s scores are:%s%s z-sc %sE(%ld)%s",
		rel_label,pad,m_msp->label,html_pre_E,ppst->zdb_size,html_post_E);
	if (ppst->zsflag > 20) {
	  fprintf(fp," E2()");
	}
      }
      header_aux(fp);
      if (m_msp->markx & MX_M9SUMM) {
	if (m_msp->show_code == SHOW_CODE_ID) {
#ifdef SHOWSIM
	  fprintf(fp," %%_id  %%_sim  alen");
#else
	  fprintf(fp," %%_id  alen");
#endif
	}
	else {
	  if (m_msp->markx & MX_HTML && m_msp->show_code !=1) { fprintf(fp,"<!-- ");}
#ifndef SHOWSIM
	  fprintf(fp,"\t%%_id  %%_gid %4s  alen  an0  ax0  pn0  px0  an1  ax1 pn1 px1 gapq gapl  fs ",m_msp->f_id1);
#else
	  fprintf(fp,"\t%%_id  %%_sim %4s  alen  an0  ax0  pn0  px0  an1  ax1 pn1 px1 gapq gapl  fs ",m_msp->f_id1);
#endif
	}
	if (m_msp->show_code == SHOW_CODE_ALIGN) { fprintf(fp," aln_code"); }
	if (m_msp->markx & MX_HTML && m_msp->show_code!=1) { fprintf(fp," -->");}
      }
      fprintf(fp,"\n");
    }
    else {
      fprintf(fp,"\nThe best%s scores are:%s%s",rel_label,pad,m_msp->label);
      header_aux(fp);
      if (m_msp->markx & MX_M9SUMM) {
	if (m_msp->show_code == SHOW_CODE_ID) {
#ifdef SHOWSIM
	  fprintf(fp," %%_id  %%_sm  alen");
#else
	  fprintf(fp," %%_id  alen");
#endif
	}
	else {
#ifndef SHOWSIM
	  fprintf(fp,"\t%%_id  %%_gid %4s  alen  an0  ax0  pn0  px0  an1  ax1 pn1 px1 gapq gapl  fs ",m_msp->f_id1);
#else
	  fprintf(fp,"\t%%_id  %%_sim %4s  alen  an0  ax0  pn0  px0  an1  ax1 pn1 px1 gapq gapl  fs ",m_msp->f_id1);
#endif	/* SHOWSIM */
	}
      }
      if (m_msp->show_code == SHOW_CODE_ALIGN) {	fprintf(fp," aln_code"); }
      fprintf(fp,"\n");
    }
  }	/* !(m_msp->markx & MX_M8OUT) */

  istart = 0;
l1:
  istop = min(nshow, nbest);

  for (ib=istart; ib<istop; ib++) {
    bbp = bptr[ib];
    if (ppst->do_rep) {
      bbp->repeat_thresh = 
	min(E1_to_s(ppst->e_cut_r, m_msp->n0, bbp->seq->n1,ppst->zdb_size, m_msp->pstat_void),
	    bbp->rst.score[ppst->score_ix]);
    }

#ifdef DEBUG
    if (bbp->seq->n1 != bbp->n1 ) {
      fprintf(stderr, " *** lib len error [%d!=%d] *** %s score %d\n",
	      bbp->seq->n1,bbp->n1, bbp->mseq->libstr, bbp->rst.score[0]);
    }
#endif

    /* this gets us a valid bline[] and the library for searching if necessary
       do not read if we have a long enough bline or we don't need a sequence 
    */
    if (bbp->mseq->bline != NULL && bbp->mseq->bline_max >= m_msp->aln.llen) {
      ranlib_done = 0;

      /* copy m_msp->aln.llen, not llen-r_margin, because the r_margin
	 will be set later, possibly after the gi|12345 is removed */
      strncpy(bline,bbp->mseq->bline,m_msp->aln.llen);
      bline[m_msp->aln.llen]='\0';
    }
    else {
      if ((m_fptr=re_openlib(bbp->mseq->m_file_p,!m_msp->quiet))==NULL) {
	fprintf(stderr,"*** cannot re-open %s\n",bbp->mseq->m_file_p->lb_name);
	exit(1);
      }
      RANLIB(bline,m_msp->aln.llen,bbp->mseq->lseek,bbp->mseq->libstr,m_fptr);
      ranlib_done = 1;
    }

  /* l_name is used to build an HTML link from the bestscore line to
     the alignment.  It can also be used to discriminate multiple hits
     from the same long sequence.  This requires that fast_pan use -m 6. */

    strncpy(l_name,bline,sizeof(l_name)); /* get rid of text after second "|" */
    l_name[sizeof(l_name)-1]='\0';
    if ((bp=strchr(l_name,' '))!=NULL) *bp=0;
    if ((bp=strchr(&l_name[3],'|'))!=NULL) *bp='\0';
    if (m_msp->nframe > 2) sprintf(&l_name[strlen(l_name)],"_%d",bbp->frame+1);
    else if (m_msp->nframe > 0 && bbp->frame == 1)
      SAFE_STRNCAT(l_name,"_r",sizeof(l_name));
    if (bbp->mseq->cont-1 > 0) {
      sprintf(tmp_str,":%d",bbp->mseq->cont-1);
      SAFE_STRNCAT(l_name,tmp_str,sizeof(l_name));
    }

    /* get a valid cur_ares_p chain and put it in bbp->ares */
    if (m_msp->stages>1 || m_msp->markx & MX_M9SUMM) {	/* we need a sequence */
      if (bbp->seq->aa1b == NULL || (m_msp->ann_flg && bbp->seq->aa1_ann==NULL)) {
	if (!ranlib_done) {	/* we didn't open the library already */
	  if ((m_fptr=re_openlib(bbp->mseq->m_file_p,!m_msp->quiet))==NULL) {
	    fprintf(stderr,"*** cannot re-open %s\n",bbp->mseq->m_file_p->lb_name);
	    exit(1);
	  }
	  RANLIB(bline,m_msp->aln.llen,bbp->mseq->lseek,bbp->mseq->libstr,m_fptr);
	  ranlib_done = 1;
	}
	n1 = re_getlib(aa1save,
		       m_msp->ann_flg ? &(bbp->seq->aa1_ann) : NULL, 
		       maxn,m_msp->ldb_info.maxt3,
		       m_msp->ldb_info.l_overlap,bbp->mseq->cont,m_msp->ldb_info.term_code,
		       &bbp->seq->l_offset,&bbp->seq->l_off,bbp->mseq->m_file_p);


	aa1 = aa1save;
      }
      else {
	n1 = bbp->seq->n1;
	aa1 = bbp->seq->aa1b;
      }

      if (n1 != bbp->n1) {
	fprintf(stderr," *** sequence length conflict %d != %d: %s\n", n1, bbp->n1, bline);
	continue;
      }

      if ( m_msp->stages > 1 && bbp->rst.score[2] == -BIGNUM) { 
	/* this is not typically done unless m_msp->stages > 1 */
	do_opt (aa0[bbp->frame], m_msp->n0, aa1, n1, bbp->frame, ppst, f_str[bbp->frame], &rst);
	bbp->rst.score[2]=rst.score[2];
      }

      if (!bbp->have_ares & 0x1) {
	bbp->a_res = build_ares_code(aa0[bbp->frame], m_msp->n0, aa1, bbp->seq,
				     bbp->frame, &bbp->have_ares,
				     bbp->repeat_thresh, m_msp, ppst, f_str[bbp->frame] );
      }
    }	/* end stages > 1 || MX_M9SUMM9 */

    n1tot = (bbp->mseq->n1tot_p) ? *bbp->mseq->n1tot_p : bbp->seq->n1;

    bline_p = bline;
    if ((m_msp->markx & (MX_HTML + MX_M8OUT+ MX_MBLAST)) && !strncmp(bline,"gi|",3)) {
      bline_p = strchr(bline+4,'|')+1;
      *(bline_p-1) = 0;
      gi_num = atoi(bline+3);
    }

    if (m_msp->markx & MX_M8OUT) {
      if ((bp=strchr(bline_p,' '))!=NULL) *bp = '\0';
    }
    else {
      bline_p[m_msp->aln.llen-r_margin]='\0';
      /* check for translated frame info */
      if (m_msp->nframe > -1) bline_p[m_msp->aln.llen-(r_margin+4)]='\0';
    }
    /* now its time to report the summary numbers for all the alignments */

    /* in the next loop, cur_ares_p could be NULL if we haven't done do_walign() */
    cur_ares_p = bbp->a_res;

    first_line = 1;
    do {
      /* if cur_res_p != NULL, then we get rst from a_res->rst
	 Otherwise, it comes from bbp->rst
      */

      if (cur_ares_p && !first_line) {
	rst_p = &cur_ares_p->rst;
      }
      else {
	rst_p = &bbp->rst;
      }

      n1 = bbp->seq->n1;
      l_score0 = rst_p->score[ppst->score_ix];
      lzscore = find_z(l_score0, rst_p->escore, n1, rst_p->comp, m_msp->pstat_void);
      if (ppst->zsflag > 20) {
	lzscore2 = find_z(l_score0, rst_p->escore, n1, rst_p->comp, m_msp->pstat_void2);
      }
      lbits = zs_to_bit(lzscore, m_msp->n0, n1);

      /* *********************************** */
      /* standard "The best scores are" here */
      /* *********************************** */

      if (!(m_msp->markx & (MX_M8OUT + MX_MBLAST2))) {
	if (first_line) {
	  first_line = 0;
	  fprintf (fp, fmt,bline_p,n1tot);
	  if (m_msp->nframe > 2) fprintf (fp, " [%d]", bbp->frame+1);
	  else if (m_msp->nframe >= 0) fprintf(fp," [%c]",(bbp->frame > 0 ?'r':'f'));
	}
	else {
	  fprintf (fp, fmt2,"\n+-");
	}

	if (m_msp->srelv == 1) fprintf (fp, " %4d", rst_p->score[ppst->score_ix]);
	else {
	  if (m_msp->srelv-1 > 0) fprintf (fp, " %4d", rst_p->score[0]);
	  if (m_msp->srelv-1 > 1 || m_msp->stages>1)
	    fprintf (fp, " %4d", rst_p->score[1]);
	  fprintf (fp, " %4d", rst_p->score[ppst->score_ix]);
	}

	if (ppst->zsflag>=0) { 
	  if (m_msp->z_bits==1) {
	    fprintf (fp, " %.1f %s%7.2g%s",lbits,html_pre_E,
		     zs_to_E(lzscore, n1, ppst->dnaseq, ppst->zdb_size, m_msp->db),
		     html_post_E);
	    if (ppst->zsflag > 20) {
	      fprintf (fp, " %7.2g",zs_to_E(lzscore2, n1, ppst->dnaseq, ppst->zdb_size, m_msp->db));
	    }
	  }
	  else {
	    fprintf (fp, " %.1f %s%7.2g%s",lzscore,html_pre_E,
		     zs_to_E(lzscore, n1, ppst->dnaseq, ppst->zdb_size, m_msp->db),
		     html_post_E);
	    if (ppst->zsflag > 20) {
	      fprintf (fp, " %7.2g",zs_to_E(lzscore2, n1, ppst->dnaseq, ppst->zdb_size, m_msp->db));
	    }
	  }
	}
	show_aux(fp,bbp);
      }
      else if (m_msp->markx & MX_M8OUT) {	/* MX_M8OUT -- provide query, library */
	if (first_line) {first_line = 0;}
	fprintf (fp,"%s\t%s",m_msp->qtitle,bline_p);
      }
      else if (m_msp->markx & MX_MBLAST2) {	/* blast "Sequences producing" */ 
	if (first_line) {first_line = 0;}
	fprintf (fp,"%-67s %6.1f    %.1g", bline_p, lbits,
		    zs_to_E(lzscore,n1,ppst->dnaseq,ppst->zdb_size,m_msp->db));
      }

      if (m_msp->markx & MX_M9SUMM) {
	loffset = bbp->seq->l_offset;
	l_off = bbp->seq->l_off;
	aln_p = &cur_ares_p->aln;
	seq_code = cur_ares_p->aln_code;
	seq_code_len = cur_ares_p->aln_code_n;
	ann_code = cur_ares_p->ann_code;
	ann_code_len = cur_ares_p->ann_code_n;

	if (aln_p->lc > 0) percent = (100.0*(float)aln_p->nident)/(float)aln_p->lc;
	else percent = -100.00;

	ngap = cur_ares_p->aln.ngap_q + cur_ares_p->aln.ngap_l;
#ifndef SHOWSIM
	if (aln_p->lc-ngap > 0) gpercent = (100.0*(float)aln_p->nident)/(float)(aln_p->lc-ngap);
	else gpercent = -100.00;
#else
	if (aln_p->lc-ngap > 0) gpercent = (100.0*(float)cur_ares_p->aln.nsim)/(float)(aln_p->lc);
	else gpercent = -100.00;
#endif	/* SHOWSIM */

	if (m_msp->show_code != SHOW_CODE_ID) {	/* show more complete info than just identity */

	  /*  	calc_astruct(aln_p, cur_ares_p); -- this function
		should not be used after calc_code or any other
		alignment that calculates amax0/amax1 */

	  /* we need the coordinates for annotated SHOW_CODE_ALIGN */
	  calc_coord(m_msp->n0,bbp->seq->n1,
		     m_msp->q_offset + (m_msp->q_off-1) + (m_msp->sq0off-1),
		     loffset + (l_off-1) + (m_msp->sq1off-1),
		     aln_p);

	  if (m_msp->markx & MX_HTML) fprintf(fp,"<!-- ");
	  /*            %_id  %_sim s-w alen an0  ax0  pn0  px0  an1  ax1  pn1  px1 gapq gapl fs  */
	  /*                    alignment    min  max            min  max */
	  /*                    sequence coordinate    min  max            min  max */
	  if (!(m_msp->markx & MX_M8OUT)) {
	    fprintf(fp,"\t%5.3f %5.3f %4d %4d %4ld %4ld %4ld %4ld %4ld %4ld %4ld %4ld %3d %3d %3d",
		    percent/100.0,gpercent/100.0, 
		    cur_ares_p->sw_score,
		    aln_p->lc,
		    aln_p->d_start0,aln_p->d_stop0,
		    aln_p->q_start_off, aln_p->q_end_off,
		    aln_p->d_start1,aln_p->d_stop1,
		    aln_p->l_start_off, aln_p->l_end_off,
		    aln_p->ngap_q,aln_p->ngap_l,aln_p->nfs);
	    if ((m_msp->show_code & SHOW_CODE_ALIGN) == SHOW_CODE_ALIGN
		&& seq_code_len > 0 && seq_code != NULL) {
	      fprintf(fp,"\t%s",seq_code);
	      if (ann_code_len > 0 && ann_code != NULL) {
		fprintf(fp,"\t%s",ann_code);
	      }
	    }
	  }
	  else {	/* MX_M8OUT -- blast order, tab separated */
	    fprintf(fp,"\t%.2f\t%d\t%d\t%d\t%ld\t%ld\t%ld\t%ld\t%.2g\t%.1f\n",
		    percent,aln_p->lc,aln_p->nmismatch,
		    aln_p->ngap_q + aln_p->ngap_l+aln_p->nfs,
		    aln_p->d_start0, aln_p->d_stop0,
		    aln_p->d_start1, aln_p->d_stop1,
		    zs_to_E(lzscore,n1,ppst->dnaseq,ppst->zdb_size,m_msp->db),
		    lbits);
	  }
	}
	else {	/* !SHOW_CODE */
#ifdef SHOWSIM
	  fprintf(fp," %5.3f %5.3f %4d", percent/100.0,(float)aln_p->nsim/(float)aln_p->lc,aln_p->lc);
#else
	  fprintf(fp," %5.3f %4d", percent/100.0,aln_p->lc);
#endif
	}
      }
    } while ( cur_ares_p && (cur_ares_p = cur_ares_p->next));

    if (m_msp->markx & MX_HTML) fprintf(fp," <a href=\"#%s\">align</a>",l_name);
    if (!(m_msp->markx & MX_M8OUT)) fprintf(fp, "\n");
    fflush(fp);
  }

  if (quiet==0) {
    printf(" More scores? [0] ");
    fflush(stdout);
    if (fgets(rline,20,stdin)==NULL) exit(0);
    ntmp = 0;
    if (rline[0]!='\n' && rline[0]!=0) sscanf(rline,"%d",&ntmp);
    if (ntmp<=0) ntmp = 0;
    if (ntmp>0) {
      istart = istop;
      nshow = min(nshow+ntmp, nbest);
      goto l1;
    }
  }	/* end of for (ib) loop */

  if (m_msp->markx & MX_MBLAST2) {fprintf(fp, "\n\n");}

  m_msp->nshow = nshow;	/* save the number of hits displayed for showalign */

  if (m_msp->markx & MX_HTML) fprintf(fp,"</pre><hr>\n");
}
