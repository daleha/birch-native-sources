
/* thr_bufs.h - structures for passing buffers of sequences to threads */

/* copyright (c) 2007 William R. Pearson and the U. of Virginia */

/* $Id: thr_buf_structs.h 743 2011-05-18 23:34:34Z wrp $ */
/* $Revision: 743 $  */

#include <sys/types.h>

struct thr_str {
  int worker;
  void *status;
  int max_work_buf;
  int qframe;
  struct pstruct *ppst;
  const struct mngmsg *m_msp;
  int qshuffle;
  int nrshuffle;
  unsigned char *aa0;	/* query sequence */
  int n0;		/* query sequence length */
  int nm0;
  int max_tot;
  char info_lib_range[MAX_SSTR];
  void **f_str_ap;
};


/* this structure passes library sequences to the worker threads
   and returns scores */

struct buf2_hdr_s {
  int buf2_cnt;		/* number of buf2 records */
  int buf2_type;	/* search, shuffle, etc */
  int stop_work;
  int have_data;
  int have_results;
  int have_best_save;
  int shuff_cnt;
  int worker_idx;
  struct seq_record *seq_b;	/* pointer to contig. actual sequence data */
  struct mseq_record *mseq_b;	/* pointer to contig. seq. meta-info */
  struct seqr_chain *my_chain;	/* pointer to the seqr_chain providing data for this buffer */
  int seqr_cnt;		/* count of seq_records, which can be less than buf2_cnt */
  int seq_record_continuous;	/* contiguous seq_record/aa1b */
  unsigned char *aa1b_start;	/* pointer to contiguous aa1b buffer */
  int aa1b_size;	/* allocated size of aab buffer -- needed for PCOMPLIB */
  int aa1b_used;	/* used size of aab buffer -- needed for PCOMPLIB */
  int my_id;
  int my_worker;
};


/* this structure contains a single sequence record, with all the
   information necessary to calculate a score */

struct buf2_data_s {
  struct seq_record *seq;	/* pointer to sequence */
  struct mseq_record *mseq;	/* pointer to sequence meta data */
  int frame;			/* query frame used for returning results, indexes into aa0[], f_str[]
				   also used in best_stats.h bbp->frame */
  int repeat_thresh;		/* threshold for sub-alignment */
  int stats_idx;		/* where to save for statistics */
  int seq_dup;			/* duplicate entry for alternate frame */
  struct beststr *best_save;	/* if beststr is pointing to this record, where is it saved? */
};

struct buf2_res_s {
  struct rstruct rst;
  int qr_score;
  struct rstruct r_rst;
  double qr_escore;
};

struct buf2_ares_s {
  int have_ares;
  struct a_res_str *a_res;
  int best_idx;
};

#define BUF2_DOWORK 0x1
#define BUF2_DOSHUF 0x2
#define BUF2_DOOPT 0x4
#define BUF2_DOALIGN 0x8

struct buf_head {
  struct buf2_hdr_s hdr;		/* meta-information */
  struct buf2_data_s *buf2_data;	/* input (sequence) datat */
  struct buf2_res_s *buf2_res; 		/* score type results */
  struct buf2_ares_s *buf2_ares;	/* alignment results */
  struct score_count_s s_cnt_info;	/* statitics info */
};
