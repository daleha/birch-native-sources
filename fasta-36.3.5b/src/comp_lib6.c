/* copyright (c) 1996, 1997, 1998, 1999, 2002  William R. Pearson and the
   U. of Virginia */

/*  $Id: comp_lib6.c 852 2011-10-27 15:14:53Z wrp $ */
/*  $Revision: 852 $  */

/*
 * Jan 17, 2007 - remove #ifdef PRSS - begin better statistics in place
 * for small libraries, related libraries
 *
 * Concurrent read version
 *
 *	Feb 20, 1998 modifications for prss3
 *
 *	December, 1998 - DNA searches are now down with forward and reverse
 *			 strands
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#include <limits.h>
#include <float.h>
#include <math.h>

#ifdef UNIX
#include <unistd.h>
#include <sys/types.h>
#include <signal.h>
#endif

#ifdef MPI_SRC
#include "mpi.h"
#endif

#include "defs.h"
#include "mm_file.h"

#include "best_stats.h"			/* defines beststr */
#include "structs.h"		/* mngmsg, libstruct */

#include "thr_buf_structs.h"
#include "drop_func.h"

#define XTERNAL
#include "uascii.h"

#ifdef PCOMPLIB
char *mp_verstr="version 36.3.3preload Feb, 2011 MPI";
#else
char *mp_verstr="version 36.3.3preload Feb, 2011";
#endif

/********************************/
/* global variable declarations */
/********************************/

extern int fa_max_workers;

/********************************/
/* extern variable declarations */
/********************************/
extern char *prog_func;		/* function label */
extern char *verstr, *iprompt0, *iprompt1, *iprompt2, *refstr;

/********************************/
/*extern function declarations  */
/********************************/
struct lmf_str *open_lib(struct lib_struct *lib_p, int dnaseq, int *sascii, int quiet);

int re_getlib(unsigned char *, unsigned char *, 
	      int, int, int, int, int, long *, long *, 
	      struct lmf_str *m_fptr);

int closelib(struct lmf_str *m_fptr);

void *my_srand();
unsigned int my_nrand(int, void *);

struct seqr_chain *
new_seqr_chain(int max_seq_len, int aa1b_size, struct seqr_chain *old_chain, int maxn);
void end_seqr_chain(struct seqr_chain *last_seqr);
struct seqr_chain *
read_seqr_chain(struct mng_thr *m_bufi, struct lib_struct *lib_list_p,
		struct mngmsg *m_msp, struct pstruct *ppst, int fa_workers);

struct seq_record *
new_sequence_p(struct mseq_record **cur_mseq_p, struct seq_record *prev_seq_p,
	       struct seqr_chain **cur_seqr_chain, int max_seq_cnt, int maxn);

void reset_seqr_chain(struct seqr_chain *seqr_base);
struct seq_record *
next_sequence_p(struct mseq_record **, struct seqr_chain *);

void init_aa0(unsigned char **aa0, int n0, int nm0,
	      unsigned char **aa0s, unsigned char **aa1s, 
	      int qframe, int qshuffle_flg, int max_tot,
	      struct pstruct *ppst, void **f_str, void **qf_str,
	      void *my_rand_state);

extern int ann_scan(unsigned char *, int, unsigned char **, int);
extern int scanseq(unsigned char *seq, int n, char *str);
extern void re_ascii(int *qascii, int *sascii);
extern int recode(unsigned char *seq, int n, int *qascii, int nsq);
extern void revcomp(unsigned char *seq, int n, int *c_nt);

extern void init_ascii(int is_ext, int *sascii, int nsq, int is_dna);
extern void validate_novel_aa(int *sascii, int nsq, int is_dna);
extern void qshuffle(unsigned char *aa0, int n0, int nm0, void *);
extern void free_pam2p(int **);

#ifdef DEBUG
int check_seq_range(unsigned char *aa1b, int n1, int nsq, char *);
#endif

/* initialize environment (doinit.c) */
extern void initenv (int argc, char **argv, struct mngmsg *m_msg,
		     struct pstruct *ppst, unsigned char **aa0);

/* print timing information */
extern void ptime (FILE *, long);

#ifdef COMP_MLIB 
#define QGETLIB (q_file_p->getlib)
#endif

#define GETLIB (m_file_p->getlib)

int samp_stats_idx (int *pre_nstats, int nstats, void *rand_state);

void
save_best(struct buf_head *lib_buf, const struct mngmsg *, struct pstruct *ppst,
	  struct db_str *, FILE *fdata, struct hist_str *, void **,
	  struct score_count_s *);

void
save_shuf(struct buf_head *lib_buf, int nitt, int shuff_max, int score_ix,
	  struct score_count_s *);

int
save_align(struct buf_head *lib_bhead_p, struct beststr **bestp_arr);

void preserve_seq(struct buf2_data_s *, struct seq_record *, struct mseq_record *, struct beststr *);
void
buf_do_work(unsigned char **aa0, int n0, struct buf_head *lib_bhead_p,
	    struct pstruct *ppst, void **f_str);
void
buf_qshuf_work(unsigned char *aa0s, int n0, struct buf_head *lib_bhead_p,
	       struct pstruct *ppst, void *qf_str, int score_ix);
void
buf_shuf_work(unsigned char **aa0, int n0, unsigned char *aa1s,
	      struct buf_head *lib_bhead_p, struct pstruct *ppst,
	      void **f_str, int score_ix, void *rand_state);
void
buf_shuf_seq(unsigned char **aa0, int n0,
	     unsigned char **aa1_shuff_b, unsigned char *aa1save, int maxn,
	     struct beststr **bestp_arr, int nbest, 
	     struct pstruct *pst, struct mngmsg *m_msp,
	     struct mng_thr *m_thr_p
#if !defined(COMP_THR) && !defined(PCOMPLIB)
	     , void **f_str
#endif
	     , struct score_count_s *s_info);

void
buf_align_seq(unsigned char **aa0, int n0,
	      struct beststr **bestp_arr, int nbest,
	      struct pstruct *ppst, struct mngmsg *m_msp,
	      struct mng_thr *m_bufi_p
#if !defined(COMP_THR) && !defined(PCOMPLIB)
	      , void **f_str
#endif
	      );

void
buf_do_align(unsigned char **aa0,  int n0,
	     struct buf_head *lib_bhead_p, 
	     struct pstruct *ppst, const struct mngmsg *m_msp,
	     void **f_str);

struct buf_head *
alloc_comp_bufs (struct mng_thr *m_bufi_p, struct mngmsg *m_msp,
		 int ave_seq_len);

/* statistics functions */
extern int
process_hist(struct stat_str *sptr, int nstats, 
	     const struct mngmsg *m_msg,
	     struct pstruct *ppst,
	     struct hist_str *hist, void **pstat_void, struct score_count_s *s_info, int do_hist);

extern void addhistz(double, struct hist_str *); /* scaleswn.c */
void selectbestz(struct beststr **, int, int );
extern double find_z(int score, double escore, int length, double comp,void *);
extern double zs_to_E(double zs,int n1, int dnaseq, long entries, struct db_str db);

void last_stats(const unsigned char *, int, 
		struct stat_str *sptr, int nstats,
		struct beststr **bestp_arr, int nbest,
		const struct mngmsg *m_msg, struct pstruct *ppst, 
		struct hist_str *histp, void *);

int last_calc( unsigned char **aa0, unsigned char *aa1, int maxn,
	       struct beststr **bestp_arr, int nbest,
	       const struct mngmsg *m_msg, struct pstruct *ppst, 
	       void **f_str, void *rs_str);

void scale_scores(struct beststr **bestp_arr, int nbest,
		  struct db_str,struct pstruct *ppst, void *);

int E1_to_s(double e_val, int n0, int n1, int db_size, void *pu);

extern int shuffle(unsigned char *, unsigned char *, int, void *);
extern int shuffle3(unsigned char *, unsigned char *, int, void *);
extern int rshuffle(unsigned char *, unsigned char *, int);
extern int wshuffle(unsigned char *, unsigned char *, int, int, void *);

extern void set_db_size(int, struct db_str *, struct hist_str *);

/* pre-alignment */
extern void 
pre_load_best(unsigned char *aa1, int maxn,struct beststr **bbp_arr, int nbest,
	      struct mngmsg *m_msp);

extern struct a_res_str *
build_ares_code(unsigned char *aa0, int n0,
		unsigned char *aa1, struct seq_record *seq,
		int frame, int *have_ares, int repeat_thresh, 
		const struct mngmsg *m_msp, struct pstruct *ppst,
		void *f_str
		);

/* display functions */
extern void
showbest (FILE *fp, unsigned char **aa0, unsigned char *aa1, int maxn,
	  struct beststr **bestp_arr, int nbest,
	  int qlib, struct mngmsg *m_msg,struct pstruct *ppst,
	  struct db_str db, char **gstring2p, void **f_str);

extern void
showalign (FILE *fp, unsigned char **aa0, unsigned char *aa1, int maxn,
	   struct beststr **bestp_arr, int nbest, int qlib, 
	   const struct mngmsg *m_msg, const struct pstruct *ppst,
	   char **gstring2p, void **f_str, struct mng_thr *m_bufi_p);

/* misc functions */
void h_init(struct pstruct *, struct mngmsg *, char *);		/* doinit.c */
void last_init(struct mngmsg *, struct pstruct *); /* initfa/sw.c */
void last_params(unsigned char *, int, struct mngmsg *, struct pstruct *);
int validate_params(const unsigned char *, int, const struct mngmsg *,
		    const struct pstruct *,
		    const int *lascii, const int *pascii);

void s_abort(char *, char *);		/* compacc.c */

/* initfa/sw.c */
void resetp(struct mngmsg *, struct pstruct *); 

void gettitle(char *, char *, int);	/* nxgetaa.c */
void lib_choice(char *lname, int nln, char *flstr, int ldnaseq);	/* lib_sel.c */
struct lib_struct *lib_select(char *lname, struct mngmsg *m_msp);	/* lib_sel.c */
void query_parm(struct mngmsg *, struct pstruct *); /* initfa/sw.c */
void selectbestz(struct beststr **, int, int);

/* compacc.c */
void prhist(FILE *, const struct mngmsg *, struct pstruct *, struct hist_str hist,
	    int nstats, int sstats, struct db_str, char *, char *, char **, char **);
void printsum(FILE *, struct db_str db);
int reset_maxn(struct mngmsg *, int, int);	/* set m_msg.maxt, maxn from maxl */

FILE *outfd;			/* Output file */

/* this information is global for fsigint() */
extern long s_time();			/* fetches time */
long tstart, tscan, tprev, tdone;	/* Timing */
#ifdef COMP_MLIB
long ttscan, ttdisp;
#endif
long tdstart, tddone;

static struct db_str qtt = {0l, 0l, 0};

#if defined(COMP_THR) || defined(PCOMPLIB)
/***************************************/
/* thread global variable declarations */
/***************************************/

/* functions for getting/sending buffers to threads (thr_sub.c) */
#ifndef PCOMPLIB
extern void init_thr(int , struct thr_str *, const struct mngmsg *, struct pstruct *,
		     unsigned char *, struct mng_thr *m_bufi_p);
extern void start_thr(void);
#define RESULTS_BUF reader_buf
#else
extern void init_thr(int , char *, const struct mngmsg *, struct pstruct *,  unsigned char *, struct mng_thr *m_bufi_p);
extern void work_comp(int);
#define RESULTS_BUF worker_buf
#endif
extern void get_rbuf(struct buf_head **lib_buf, int max_work_buf);
extern void put_rbuf(struct buf_head *lib_buf, int max_work_buf);
extern void wait_rbuf(int max_work_buf);
extern void rbuf_done(int nthreads);
extern void put_rbuf_done(int nthreads, struct buf_head *lib_buf, 
			  int max_work_buf);
#ifndef PCOMPLIB
#undef XTERNAL
#include "thr_bufs2.h"
#else
#include "pcomp_bufs.h"
#endif
#endif

struct buf_head *lib_buf2_list;

/* these variables must be global for comp_thr.c so that save_best()
   can use them */
static struct beststr **bestp_arr;	/* array of pointers */
static int nbest;	/* number of best scores */

static struct stat_str *stats; /* array of scores for statistics from real
			     (or shuffled) sequences*/
static struct stat_str *qstats;	/* array of scores for shuffled query stats */
static struct stat_str *rstats;	/* array of scores from shuffled library */

  /* these variables are global so they can be set both by the main()
     program and save_best() in threaded mode.
  */
static int nstats, nqstats, nrstats, pre_nstats, kstats, shuff_tot, sstats;
static double zbestcut;		/* cut off for best z-score */
static int bestfull;		/* index for selectbest() */
static int stats_done=0;	/* flag for z-value processing */
static void *rand_state;

static int seq_index=0;

void fsigint();

/* **************************************************************** */
/* start of main() program                                          */
/* **************************************************************** */
int
main (int argc, char *argv[]) 
{
  unsigned char *aa0[6], *aa0s;
  unsigned char *aa1save;	/* aa1shuff and aa1save must be distinct */
  unsigned char *aa1shuff, *aa1shuff_b=NULL;	/* for new unthreaded version */
  int n1;
  struct seqr_chain *lib_seqr_chain;

  int j;
  int shuff_mult, n1lib_req, jstats;
  double zs_off_save;
  int stats_inc=1;
  struct beststr *bbp;
  int n_sig;

  struct a_res_str *next_ares_p, *cur_ares_p; /* used to free-up old a_res */
  int n_pre_align;

  /* status/parameter information */
  char info_lib_range[MAX_FN];
  char *info_lib_range_p;
  char info_pgm_abbr[MAX_SSTR];
  char info_qlabel[MAX_FN];
  char *info_gstring2p[2];
  char info_gstring3[MAX_STR];
  char *info_hstring_p[2];

#ifdef COMP_MLIB
  char q_bline[MAX_STR];
  fseek_t qseek;
  int qlib;
  struct lib_struct *q_lib_p;
  struct lmf_str *q_file_p;
  int sstart, sstop, is;
#endif

  long next_q_offset;
  int lstart, lstop;
  int id;
  struct lib_struct *lib_list_p, *cur_lib_p;

  int t_best, t_rbest, t_qrbest;	/* best score of two/six frames */
  double t_escore, t_rescore, t_qrescore; /* best evalues of two/six frames */
  double db_tt;
  int utmp;		/* user input tmp */

  struct pstruct pst;
  void *f_str[6], *qf_str;	/* different f_str[]'s for forward,reverse */
  int have_f_str=0;

  /* these variables track buffers of library sequences */
  struct buf_head *lib_bhead_p, *t_lib_bhead_p;
  struct buf2_data_s *lib_buf2_dp;
  struct seq_record *current_seq_p;
  struct mseq_record *current_mseq_p;

  long ntbuff;			/* current length (in residues/bytes) of buffer */
  struct mng_thr m_bufi;	/* has max_work_buf, max_buf2_cnt */
  int ave_seq_len, buf_siz;
  /*   int empty_reader_bufs; */
#ifdef COMP_THR
  /*   int t_reader_buf_readp; */
  struct thr_str *work_info;
#endif
#ifdef MPI_SRC
  int mpi_tid;
#endif
  /* end of library sequence buffers */

  struct mngmsg m_msg;		/* Message from host to manager */
  struct hist_str hist2;	/* hist str for zsflag > 2 */
  int zsflag_save;		/* save zsflag > 20 */
  int itt;			/* index into library names */
  char rline[MAX_FN];
  char argv_line[MAX_STR];
  int t_quiet;

  int i;

  FILE *fdata=NULL;		/* file for full results */
  struct beststr *best;		/* array of best scores */

  /* save sequence meta info for sequences that are not currently available */
  struct seq_record *best_seqs;
  struct mseq_record *best_mseqs;

  int leng;			/* leng is length of the descriptive line */
  int maxn;			/* size of the library sequence examined */
  int qlcont;			/* continued query sequence */
  char *bp;			/* general purpose string ptr */
  
  /* this is necessary because of an SGI Irix 64 issue */
  info_gstring2p[0] = calloc(MAX_STR,sizeof(char));
  info_gstring2p[1] = calloc(MAX_STR,sizeof(char));
  info_hstring_p[0] = calloc(MAX_STR,sizeof(char));
  info_hstring_p[1] = calloc(MAX_STR,sizeof(char));

  if ((bp = strrchr(argv[0],'/'))!=NULL) {
    strncpy(m_msg.pgm_name,bp+1,sizeof(m_msg.pgm_name));
  }
  else {
    strncpy(m_msg.pgm_name,argv[0],sizeof(m_msg.pgm_name));
  }

  /* Initialization */

  m_msg.s_info.s_cnt[0] = m_msg.s_info.s_cnt[1] =
    m_msg.s_info.s_cnt[2] =   m_msg.s_info.tot_scores = 0;
  m_msg.ss_info.s_cnt[0] = m_msg.ss_info.s_cnt[1] =
    m_msg.ss_info.s_cnt[2] =   m_msg.ss_info.tot_scores = 0;

#ifndef SHOW_HELP
#if defined(UNIX)
  m_msg.quiet= !isatty(1);
#else
  m_msg.quiet = 0;
#endif
#else
  m_msg.quiet = 1;
#endif

#ifdef MPI_SRC
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&mpi_tid);
  if (mpi_tid > 0) {
    work_comp(mpi_tid); 
    MPI_Finalize();
    exit(0);
  }
#endif

#ifdef PGM_DOC
  argv_line[0]='\0';
  for (i=0; i<argc; i++) {
    strncat(argv_line," ",sizeof(argv_line)-strlen(argv_line)-1);
    if (strchr(argv[i],' ')) {
      strncat(argv_line,"\"",sizeof(argv_line)-strlen(argv_line)-1);
      strncat(argv_line,argv[i],sizeof(argv_line)-strlen(argv_line)-1);
      strncat(argv_line,"\"",sizeof(argv_line)-strlen(argv_line)-1);
    }
    else {
      strncat(argv_line,argv[i],sizeof(argv_line)-strlen(argv_line)-1);
    }
  }
  argv_line[sizeof(argv_line)-1]='\0';
#endif

  /* first initialization routine - nothing is known */
  h_init(&pst, &m_msg, info_pgm_abbr);
  
  m_msg.db.length = m_msg.ldb.length = qtt.length = 0l;
  m_msg.db.entries = m_msg.db.carry = 
    m_msg.ldb.entries = m_msg.ldb.carry = qtt.entries = qtt.carry = 0;
  m_msg.pstat_void = m_msg.pstat_void2 = NULL;
  m_msg.hist.entries = 0;

  f_str[5] = f_str[4] = f_str[3] = f_str[2] = f_str[1] = f_str[0] = NULL;
  aa0[0] = NULL;
  rand_state = my_srand();

  /* second initialiation - get commmand line arguments */
  initenv (argc, argv, &m_msg, &pst,&aa0[0]);

#ifdef PGM_DOC
  if (!(m_msg.markx & MX_M8OUT)) fprintf(stdout, "#%s\n",argv_line);
#endif

  if (m_msg.markx & MX_M11OUT) {
    fprintf(stdout, "#:lav\n\nd {\n   \"%s\"\n}\n",argv_line+1);
  }

#ifndef PCOMPLIB
#ifdef COMP_THR
  if ((work_info=
       (struct thr_str *)calloc(fa_max_workers,sizeof(struct thr_str)))==NULL) {
    fprintf(stderr, " cannot allocate work_info[%d]\n",fa_max_workers);
    exit(1);
  }
#else
  fa_max_workers = 1;
#endif
#endif

  ttscan = ttdisp = 0;
  tstart = tscan = s_time();
  tdstart = time(NULL);

  /* Allocate space for the query and library sequences */
  /* pad aa0[] with an extra SEQ_PAD chars for ALTIVEC padding */
  if (aa0[0]==NULL) {
    if ((aa0[0] = (unsigned char *)malloc((m_msg.max_tot+1+SEQ_PAD)*sizeof(unsigned char)))
	== NULL)
      s_abort ("Unable to allocate query sequence", "");
    *aa0[0]=0;
    aa0[0]++;
  }
  aa0[5]=aa0[4]=aa0[3]=aa0[2]=aa0[1]=aa0[0];

  if ((aa1save = (unsigned char *)malloc((m_msg.max_tot+1)*sizeof (char))) == NULL) {
    s_abort ("Unable to allocate library overlap", "");
  }
  *aa1save=0;
  aa1save++;

  if (m_msg.markx & MX_HTML) {
#ifdef HTML_HEAD    
    fprintf(stdout,"<html>\n<head>\n<title>%s Results</title>\n</head>\n<body>\n",prog_func);
#endif
    fprintf(stdout,"<pre>\n");
  }

  if (m_msg.std_output) {
    fprintf(stdout,"%s\n",iprompt0);
    fprintf(stdout," version %s\nPlease cite:\n %s\n",mp_verstr,refstr);
  }
  if (m_msg.markx & MX_HTML) fputs("</pre>\n",stdout);

  /* get query library name if not in argv[1] */
  if (m_msg.tname[0] == '\0') {
      if (m_msg.quiet == 1)
	s_abort("Query sequence undefined","");
    l1:	fputs (iprompt1, stdout);
      fflush  (stdout);
      if (fgets (m_msg.tname, MAX_FN, stdin) == NULL)
	s_abort ("Unable to read query library name","");
      m_msg.tname[MAX_FN-1]='\0';
      if ((bp=strchr(m_msg.tname,'\n'))!=NULL) *bp='\0';
      if (m_msg.tname[0] == '\0') goto l1;
  }

  /* **************************************************************** */
  /* (1) open the query library; 
     (2) get a sequence;
     (3) check for annotations */

  /* we need a q_lib_p before opening the library */
  if ((q_lib_p = (struct lib_struct *)calloc(1,sizeof(struct lib_struct)))==NULL) {
    s_abort(" cannot allocate q_lib_p","");
  }
  else {
    q_lib_p->file_name = m_msg.tname;
  }

  /* Open query library */
  if ((q_file_p= open_lib(q_lib_p, m_msg.qdnaseq,qascii,!m_msg.quiet))==NULL) {
    s_abort(" cannot open library ",m_msg.tname);
  }
  /* Fetch first sequence */
  qlib = 0;
  m_msg.q_offset = next_q_offset = 0l;
  qlcont = 0;
  m_msg.n0 = 
    QGETLIB (aa0[0], MAXTST, m_msg.qtitle, sizeof(m_msg.qtitle),
	     &qseek, &qlcont,q_file_p,&m_msg.q_off);
  if ((bp=strchr(m_msg.qtitle,' '))!=NULL) *bp='\0';
  strncpy(info_qlabel,m_msg.qtitle,sizeof(info_qlabel));
  if (bp != NULL) *bp = ' ';
  info_qlabel[sizeof(info_qlabel)-1]='\0';

  /* if annotations are included in sequence, remove them */
  if (m_msg.ann_flg) {
    m_msg.n0 = ann_scan(aa0[0],m_msg.n0,&m_msg.aa0a,m_msg.qdnaseq);
  }

  /* if protein and ldb_info.term_code set, add '*' if not there */
  if (m_msg.ldb_info.term_code && !(m_msg.qdnaseq==SEQT_DNA || m_msg.qdnaseq==SEQT_RNA) &&
      aa0[0][m_msg.n0-1]!='*') {
    aa0[0][m_msg.n0++]='*';
    aa0[0][m_msg.n0]=0;
  }

  /* check for subset */
  if (q_file_p->opt_text[0]!='\0') {
    if (q_file_p->opt_text[0]=='-') {
      sstart=0; sscanf(&q_file_p->opt_text[1],"%d",&sstop);
    }
    else {
      sscanf(&q_file_p->opt_text[0],"%d-%d",&sstart,&sstop);
      sstart--;
      if (sstop <= 0 ) sstop = BIGNUM;
    }

    for (id=0,is=sstart; is<min(m_msg.n0,sstop); ) {
      aa0[0][id++]=aa0[0][is++];
    }
    aa0[0][id]=0;
    m_msg.n0 = min(m_msg.n0,sstop)-sstart;
    m_msg.q_off += sstart;
  }

  /* check to see if query has been segmented */
  if (qlcont) {
    next_q_offset = m_msg.q_offset + m_msg.n0 - m_msg.q_overlap;
  }
  else {
    next_q_offset = 0l;
  }

  /* this probably cannot happen any more */
  if (m_msg.n0 > MAXTST) {
    fprintf(stderr," sequence truncated to %d\n %s\n",MAXTST,m_msg.sqnam);
    if (m_msg.std_output) 
      fprintf(stdout," sequence truncated to %d\n %s\n",MAXTST,m_msg.sqnam);
    aa0[0][MAXTST]='\0';
    m_msg.n0=MAXTST;
  }

  /* check for protein/DNA alphabet type */
  if (m_msg.qdnaseq == SEQT_UNK) {
  /* cannot change the alphabet mapping if a matrix has been set */
  /* do automatic sequence recognition,but only for sequences > 20 residues */
    if ( !pst.pam_set && m_msg.n0 > 20 &&
	(float)scanseq(aa0[0],m_msg.n0,"ACGTUNacgtun")/(float)m_msg.n0 >0.85) {
      pascii = nascii;
      m_msg.qdnaseq = SEQT_DNA;
    }
    else {	/* its protein */
      pascii = aascii;
      m_msg.qdnaseq = SEQT_PROT;
    }
    /* modify qascii to use encoded version 
       cannot use memcpy() because it loses annotations 
    */
    re_ascii(qascii,pascii);
    init_ascii(pst.ext_sq_set,qascii,pst.nsq,m_msg.qdnaseq);
    validate_novel_aa(qascii, pst.nsq, m_msg.qdnaseq);
    m_msg.n0 = recode(aa0[0],m_msg.n0,qascii, pst.nsqx);
  }

   /* check sequence length -- cannot do before now because query
      alphabet may change */
  if (m_msg.n0 <= 0)
    s_abort ("Query sequence length <= 0: ", m_msg.tname);

   /* reset algorithm parameters for alphabet */
  resetp (&m_msg, &pst);

#ifndef COMP_MLIB
  gettitle(m_msg.tname,m_msg.qtitle,sizeof(m_msg.qtitle));
  if (m_msg.tname[0]=='-' || m_msg.tname[0]=='@') {
    strncmp(m_msg.tname,m_msg.qtitle,sizeof(m_msg.tname));
    if ((bp=strchr(m_msg.tname,' '))!=NULL) *bp='\0';
  }
#endif

   /* get library file names from argv[2] or by prompting */

  if (strlen (m_msg.lname) == 0) {
    if (m_msg.quiet == 1) s_abort("Library name undefined","");
    lib_choice(m_msg.lname,sizeof(m_msg.lname),m_msg.flstr, m_msg.ldb_info.ldnaseq);
  }
  
   /* map library abbreviation to file name(s) */
  lib_list_p = lib_select(m_msg.lname, &m_msg);

  /* Get additional parameters here */
  if (!m_msg.quiet) query_parm (&m_msg, &pst);
  
  last_init(&m_msg, &pst);

  /* Allocate space for saved scores */
  if ((best = 
       (struct beststr *)calloc((MAX_BEST+1),sizeof(struct beststr)))==NULL)
    s_abort("Cannot allocate best struct","");
  if ((bestp_arr = 
       (struct beststr **)malloc((MAX_BEST+1)*sizeof(struct beststr *)))==NULL)
    s_abort("Cannot allocate bestp_arr","");

  /* initialize high score boundary */
  bestp_arr[0] = &best[0];
  best[0].rst.score[0]=best[0].rst.score[1]=best[0].rst.score[2]= INT_MAX;
  best[0].rst.escore=FLT_MIN;	/* for E()-values, lower is best */
  best[0].zscore=FLT_MAX;	/* for Z-scores, bigger is best */
  best++; bestp_arr++;

  /* save best score sequence info */
  if ((best_seqs =
       (struct seq_record *)calloc((MAX_BEST+1),sizeof(struct seq_record)))==NULL)
    s_abort("Cannot allocate best_seqs","");

  if ((best_mseqs =
       (struct mseq_record *)calloc((MAX_BEST+1),sizeof(struct mseq_record)))==NULL)
    s_abort("Cannot allocate best_seqs","");

  /* allocate space for sampled scores */
  if ((stats =
       (struct stat_str *)calloc(MAX_STATS,sizeof(struct stat_str)))==NULL)
    s_abort ("Cannot allocate stats struct","");

  /* allocate space for shuffled library scores */
  if ((rstats =
       (struct stat_str *)calloc(m_msg.shuff_max,sizeof(struct stat_str)))==NULL)
    s_abort ("Cannot allocate rstats struct","");

#ifdef UNIX
  /* set up signals now that input is done */
  signal(SIGHUP,SIG_IGN);
#endif

  /* **************************************************************** */
  /* begin setting things up for threads */
  /* **************************************************************** */
  /* 
     This section defines m_bufi.max_buf2_cnt, the average number of entries
     per buffer, and m_bufi.max_work_buf, the total number of buffers

     Use a 2 Mbyte (DEF_WORKER_BUF) buffer for each worker.  For
     proteins, that means 5,000 sequences of length 400 (average).
     For DNA, that means 2,000 sequences of length 1000.

     To accommodate larger libraries in memory, use more buffers, not
     bigger buffers.

     Once m_bufi.max_buf2_cnt/max_work_buf and ave_seq_len are set,
     allocate all the communication buffers with alloc_comp_bufs();
  */

  if (m_msg.ldb_info.ldnaseq== SEQT_DNA) {
    ave_seq_len = AVE_NT_LEN;
    m_bufi.max_buf2_cnt = DEF_WORKER_BUF/AVE_NT_LEN;
  }
  else {
    ave_seq_len = AVE_AA_LEN;
    m_bufi.max_buf2_cnt = DEF_WORKER_BUF/AVE_AA_LEN;
  }

  /* however - buffer sizes should be a function of the number of
     workers so that all the workers are kept busy.  Assuming a 10,000
     entry library is the smallest we want to schedule, then

  if (max_buf2_cnt > 10000/fa_max_workers) 
    max_buf2_cnt = 10000/(2*fa_max_workers);
  */

  m_bufi.max_buf2_cnt /= m_msg.thr_fact;

  /* finally, max_buf2_cnt should be mod 6 for tfasta/s/f */
  m_bufi.max_buf2_cnt -= (m_bufi.max_buf2_cnt % 6);

  /* max_work_buf is the number of buffers - if the worker buffers are
     small, then make lots more buffers */

#ifdef PCOMPLIB	/* PCOMPLIB -- one buffer per worker */
  m_bufi.max_work_buf = fa_max_workers;
#else	/* !PCOMPLIB */
  m_bufi.max_work_buf = (DEF_WORKER_BUF * 2 * fa_max_workers)/(ave_seq_len * m_bufi.max_buf2_cnt);
  if (m_bufi.max_work_buf < 2*fa_max_workers) m_bufi.max_work_buf = 2*fa_max_workers;

  m_bufi.max_work_buf -= (m_bufi.max_work_buf%fa_max_workers);
#ifndef COMP_THR
  /* if not threaded, only one (larger) buffer */
  m_bufi.max_buf2_cnt *= m_bufi.max_work_buf;
  m_bufi.max_work_buf = 1;
#endif	/* !COMP_THR */
#endif

  /* allocate lib_buf2_lib[] and the associated data and results buffers,
     as well as reader_buf[] and worker_buf[] */

  lib_buf2_list = alloc_comp_bufs(&m_bufi, &m_msg, ave_seq_len);

  /* initialization of global variables for threads/buffers */

#if defined(COMP_THR) || defined(PCOMPLIB)
#ifdef DEBUG
  /*   fprintf(stderr," max_work_buf: %d\n", m_bufi.max_work_buf); */
#endif
  num_reader_bufs = m_bufi.max_work_buf;
#endif

#ifdef COMP_THR
  num_worker_bufs = 0;
  reader_done = 0;
  worker_buf_workp = 0;
  worker_buf_readp = 0;
  reader_buf_workp = 0;
  reader_buf_readp = 0;

  start_thread = 1;	/* keeps threads from starting */
#endif

  /* Label the output */
  if ((bp = (char *) strchr (m_msg.lname, ' ')) != NULL) *bp = '\0';
  if (m_msg.ltitle[0] == '\0') {
    strncpy(m_msg.ltitle,m_msg.lname,sizeof(m_msg.ltitle));
    m_msg.ltitle[sizeof(m_msg.ltitle)-1]='\0';
  }

  if (m_msg.dfile[0]) {
    fdata=fopen(m_msg.dfile,"w");
    fprintf(fdata, "#%s\n",argv_line);
  }

#ifdef COMP_MLIB
  if (m_msg.std_output) {
    printf("Query: %s\n", m_msg.tname);
  /*   if (m_msg.nln > 0) printf("searching %s library\n\n",m_msg.lbnames[0]); */
  }
#endif

  /* pre-load the library into a [m]seq_record array */

  m_msg.db.length = 0l;
  m_msg.db.entries = m_msg.db.carry = 0;

  /* also sets ldb_info.l_overlap, use a fixed 150 residue overlap */
  m_msg.ldb_info.maxn = maxn = reset_maxn(&m_msg, 150, m_msg.max_tot);
  pst.maxlen = maxn;

  seq_index = 0;

  outfd = stdout;  

  /* open the files lin lib_list_p, start reading the sequences */
  lib_seqr_chain = read_seqr_chain(&m_bufi, lib_list_p, &m_msg, &pst, fa_max_workers);

  /* main loop for doing a search, getting the next query */
  while(1) {

    /* Initialize bestp_arr */
    for (nbest = 0; nbest < MAX_BEST; nbest++)
      bestp_arr[nbest] = &best[nbest];
    nbest = 0;

    qlib++;
    stats_done = 0;

    zbestcut = -FLT_MAX;
    nstats = nrstats = pre_nstats = shuff_tot = sstats = 0;

    /* get the last parameters */
    last_params(aa0[0],m_msg.n0, &m_msg, &pst);

    if (!validate_params(aa0[0],m_msg.n0, &m_msg, &pst,
			 lascii, pascii)) {
      fprintf(stderr," *** ERROR *** validate_params() failed:\n -- %s\n", argv_line);
      exit(1);
    }

    /*
      if our function returns approximate E()-scores (FASTS, FASTM),
      we do not need to work with raw scores and later calculate
      z-scores.  When approx. E()-scores are calculated, we still need
      various statistics structures, but we can get them immediately.
      In this case, find_z() must produce a z_score (large positive
      is good) from an e_score.
    */

    if (m_msg.escore_flg) {
      pst.zsflag_f = process_hist(stats,nstats,&m_msg,&pst,
				  &m_msg.hist,&m_msg.pstat_void,&m_msg.s_info,0);
      stats_done=1;
    }

#ifdef COMP_THR    
    init_thr(fa_max_workers, work_info, &m_msg, &pst, aa0[0], &m_bufi);
#endif
#ifdef PCOMPLIB
    info_lib_range_p = &info_lib_range[0];
    init_thr(fa_max_workers, info_lib_range_p, &m_msg, &pst, aa0[0], &m_bufi);
#endif

    /* allocate space for shuffled query scores (if needed */
    if (m_msg.qshuffle && qstats==NULL) {
      if ((qstats =
	   (struct stat_str *)calloc(m_msg.shuff_max+1,sizeof(struct stat_str)))==NULL)
	s_abort ("Cannot allocate qstats struct","");
    }
    nqstats = 0;

    if (m_msg.markx & MX_HTML) fputs("<pre>\n",stdout);

    /* rline[] is a tmp string */
    if (m_msg.qdnaseq == SEQT_DNA || m_msg.qdnaseq == SEQT_RNA) {
      strncpy(rline,(m_msg.qframe==1)? " (forward-only)" : "\0",sizeof(rline));
      rline[sizeof(rline)-1]='\0';
    }
    else rline[0]='\0';

    leng = (int)strlen(m_msg.qtitle);
    if (leng > 50) leng -= 10;

    sprintf (&m_msg.qtitle[leng], " - %d %s", m_msg.n0, m_msg.sqnam);
    m_msg.seqnm = 0;


    if (m_msg.std_output) {
#ifdef COMP_MLIB
      printf("%3d>>>%s%s\n", qlib,
	     m_msg.qtitle,
	     (m_msg.revcomp ? " (reverse complement)" : rline));
#else
      printf("%.50s: %d %s%s\n %s\n",
	     info_qlabel, m_msg.n0, m_msg.sqnam,
	     (m_msg.revcomp ? " (reverse complement)" : rline));
#endif
      /* check for annotation */
      if (m_msg.ann_flg && m_msg.aa0a != NULL) {
	printf("Annotation: ");
	for (j=0; j<m_msg.n0; j++) {
	  if (m_msg.aa0a[j] && m_msg.ann_arr[m_msg.aa0a[j]] != ' ' ) {
	    printf("|%ld:%c%c",
		   j+m_msg.q_off,m_msg.ann_arr[m_msg.aa0a[j]],pst.sq[aa0[0][j]]);
	  }
	}
	printf("\n");
      }
      printf("Library: %s", m_msg.ltitle);
      fflush(stdout);
    }
    else {
      if (m_msg.markx & (MX_M8OUT + MX_M8COMMENT)) {
	printf("# %s %s\n",prog_func,verstr);
	printf("# Query: %s\n",m_msg.qtitle);
	printf("# Database: %s\n",m_msg.ltitle);
      }
    }

    fflush(outfd);

    tprev = s_time();
  
    if (fdata) fprintf(fdata,">>>%ld %3d\t%-50s\n",qtt.entries,m_msg.n0,m_msg.qtitle);

    qtt.length += m_msg.n0;
    qtt.entries++;

#if !defined(COMP_THR) && !defined(PCOMPLIB)
    have_f_str=1;

    init_aa0(aa0, m_msg.n0, m_msg.nm0, &aa0s, &aa1shuff,
	     m_msg.qframe, m_msg.qshuffle, m_msg.max_tot,
	     &pst,  &f_str[0], &qf_str, rand_state);
    aa1shuff_b = aa1shuff-1;

    /* label library size limits */
    if (pst.n1_low > 0 && pst.n1_high < BIGNUM) 
      sprintf(info_lib_range," (range: %d-%d)",pst.n1_low,pst.n1_high);
    else if (pst.n1_low > 0) 
      sprintf(info_lib_range," (range: >%d)",pst.n1_low);
    else if (pst.n1_high < BIGNUM)
      sprintf(info_lib_range," (range: <%d)",pst.n1_high);
    else
      info_lib_range[0]='\0';
    info_lib_range[sizeof(info_lib_range)-1]='\0';
    info_lib_range_p = info_lib_range;

    lib_bhead_p = lib_buf2_list;	/* equivalent to un-threaded get_rbuf() */
#else	/* COMP_THR/PCOMPLIB */
#ifndef PCOMPLIB
    start_thr();
#endif
    /* now open the library and start reading */
    /* get a buffer and fill it up */
    get_rbuf(&lib_bhead_p,m_bufi.max_work_buf);
    /*     empty_reader_bufs++; */
#endif

    lib_bhead_p->hdr.buf2_cnt = 0;
    lib_bhead_p->hdr.have_results = 0;
    lib_bhead_p->hdr.stop_work = 0;
    lib_bhead_p->hdr.buf2_type=BUF2_DOWORK;
    lib_bhead_p->hdr.aa1b_used = 0;
    /* set shuffle flag for -z > 10 */
    if (pst.zsflag >= 10 && pst.zsflag < 20) {
      lib_bhead_p->hdr.buf2_type |= BUF2_DOSHUF;
    }

    lib_buf2_dp = lib_bhead_p->buf2_data;
    ntbuff = 0;

    /* start the search */

    while ((current_seq_p = next_sequence_p(&current_mseq_p, lib_seqr_chain))) {

#ifdef DEBUG
      if (current_seq_p->aa1b[-1] != '\0') {
	fprintf(stderr," invalid current_seq_p->aa1b[-1] = %d\n",current_seq_p->aa1b[-1]);
      }
#endif

      /*
      if (current_seq_p->n1 == 0 || current_seq_p->aa1b == NULL) {
	fprintf(stderr, " *** ERROR: empty seq_record (next_sequence_p) ***\n");
	continue;
      }
      */

      lib_buf2_dp->seq = current_seq_p;
      lib_buf2_dp->mseq = current_mseq_p;
      n1 = current_seq_p->n1;

      /* check to see whether this score (or a shuff score) should
	 be included in statistics */
      jstats = samp_stats_idx(&pre_nstats, nstats, rand_state);

      ntbuff += n1+1;
      lib_bhead_p->hdr.aa1b_used += n1+1;

#ifdef PCOMPLIB
      lib_buf2_dp->seq_dup = 0;	/* mark first ->seq as original, not duplicate */
#endif
      for (itt=m_msg.revcomp; itt<=m_msg.nitt1; itt++) {
	lib_buf2_dp->frame = itt;
	lib_buf2_dp->stats_idx = jstats;
	lib_buf2_dp++;		/* point to next buf2 */
	lib_bhead_p->hdr.buf2_cnt++;

	/* point to the current sequence */
#ifdef PCOMPLIB
	lib_buf2_dp->seq_dup = 1;	/* mark duplicates */
#endif
	lib_buf2_dp->seq = current_seq_p;
	lib_buf2_dp->mseq = current_mseq_p;
      } /* for (itt .. */

	/* if the buffer is filled */
      if (lib_bhead_p->hdr.buf2_cnt >= m_bufi.max_buf2_cnt ||
	  lib_bhead_p->hdr.aa1b_used + maxn > m_bufi.seq_buf_size) {

#if defined(COMP_THR) || defined(PCOMPLIB)  /* if COMP_THR/PCOMPLIB - fill and empty buffers */
	/* provide filled buffer to workers */
	lib_bhead_p->hdr.have_data = 1;
	put_rbuf(lib_bhead_p,m_bufi.max_work_buf);
	get_rbuf(&lib_bhead_p,m_bufi.max_work_buf);	/* get an empty buffer to fill */
#else	/* just do the searches */
	if (lib_bhead_p->hdr.buf2_type & BUF2_DOWORK) {
	  buf_do_work(aa0, m_msg.n0, lib_bhead_p, &pst, f_str);
	  if (m_msg.qshuffle)
	    buf_qshuf_work(aa0s,m_msg.n0, lib_bhead_p, &pst, qf_str, pst.score_ix);
	}
	if (lib_bhead_p->hdr.buf2_type & BUF2_DOSHUF) {
	  buf_shuf_work(aa0,m_msg.n0, aa1shuff, lib_bhead_p,
			&pst, f_str, pst.score_ix,rand_state);
	}
#endif
	/* "empty" buffers have results that must be processed */
	if (lib_bhead_p->hdr.buf2_cnt && lib_bhead_p->hdr.have_results) {
	  save_best(lib_bhead_p,&m_msg, &pst, &m_msg.ldb, fdata,
		    &m_msg.hist, &m_msg.pstat_void, &m_msg.s_info);

	  /* this section of code is only used for re-cycled buffers */
	  /* no need to preserve seq with comp_lib6.c, which keeps all seq in memory */
	  /*
	  if (lib_bhead_p->hdr.have_best_save) {
	    lib_buf2_dp = lib_bhead_p->buf2_data;
	    while (lib_bhead_p->hdr.buf2_cnt--) {
	      if (lib_buf2_dp->best_save != NULL) {
		preserve_seq(lib_buf2_dp, best_seqs, best_mseqs, best);
	      }
	      lib_buf2_dp->best_save = NULL;
	      lib_buf2_dp++;
	    }
	    lib_bhead_p->hdr.have_best_save = 0;
	  }
	  */
	}

	/* now the buffer is truly empty, fill it up */
	lib_bhead_p->hdr.buf2_cnt = 0;
	lib_bhead_p->hdr.have_results = 0;
	lib_bhead_p->hdr.stop_work = 0;
	lib_bhead_p->hdr.aa1b_used = 0;

	/* specify the kinds of work required */
	lib_bhead_p->hdr.buf2_type = BUF2_DOWORK;
	if (pst.zsflag >= 10  && pst.zsflag < 20) {
	  lib_bhead_p->hdr.buf2_type |= BUF2_DOSHUF;
	}
	/* make certain we have a valid data pointer */
	lib_buf2_dp = lib_bhead_p->buf2_data;
      }
    } /* end next_sequence() loop */

#if defined(COMP_THR) || defined(PCOMPLIB)
    /* send off the final data buffer */
    lib_bhead_p->hdr.have_data = 1;	/* ignored if buf2_cnt <= 0 */
    put_rbuf(lib_bhead_p,m_bufi.max_work_buf);

#ifndef PCOMPLIB
    info_lib_range_p = work_info[0].info_lib_range;
#endif

    /* wait for the threads to finish */
    wait_rbuf(m_bufi.max_work_buf);
    /* wait_rbuf(m_bufi.max_work_buf - empty_reader_bufs); */

    /* save the final results */
    for (i=0; i < num_reader_bufs; i++) {
      save_best(RESULTS_BUF[i],&m_msg, &pst, &m_msg.ldb, fdata,
		&m_msg.hist, &m_msg.pstat_void, &m_msg.s_info);
      RESULTS_BUF[i]->hdr.buf2_cnt = RESULTS_BUF[i]->hdr.have_results = 0;
    }
#else	/* !defined(COMP_THR || PCOMPLIB) */
    if (lib_bhead_p->hdr.buf2_type & BUF2_DOWORK) {
      buf_do_work(aa0, m_msg.n0, lib_bhead_p, &pst, f_str);
      if (m_msg.qshuffle)
	buf_qshuf_work(aa0s,m_msg.n0, lib_bhead_p, &pst, qf_str, pst.score_ix);
    }

    if (lib_bhead_p->hdr.buf2_type & BUF2_DOSHUF) {
      buf_shuf_work(aa0, m_msg.n0, aa1shuff, lib_bhead_p,
		    &pst, f_str, pst.score_ix, rand_state);
    }

    save_best(lib_bhead_p,&m_msg, &pst, &m_msg.ldb, fdata,
	      &m_msg.hist, &m_msg.pstat_void, &m_msg.s_info);
    lib_bhead_p->hdr.buf2_cnt = lib_bhead_p->hdr.have_results = 0;
#endif

    m_msg.nbr_seq = m_msg.db.entries;
    get_param(&pst, info_gstring2p,info_gstring3, &m_msg.s_info);

    /* *************************** */
    /* analyze the last results    */
    /* *************************** */
    
#ifndef SAMP_STATS
    if (!stats_done && nstats > 0) {
#endif
      /* we ALWAYS do this if SAMP_STATS, because the statistics may have changed */
      /* the new incremental sampling produces an nstats that is much
	 too large */

      zsflag_save = pst.zsflag;
      if (pst.zsflag > 20) {
	pst.zsflag -= 20;
      }
      pst.zsflag_f = process_hist(stats,nstats,&m_msg, &pst,&m_msg.hist,
				  &m_msg.pstat_void, &m_msg.s_info, stats_done);
      pst.zsflag = zsflag_save;

      if (m_msg.pstat_void != NULL) {
	stats_done = 1;
	for (i = 0; i < nbest; i++) {
	  bestp_arr[i]->zscore =
	    find_z(bestp_arr[i]->rst.score[pst.score_ix],
		       bestp_arr[i]->rst.escore, bestp_arr[i]->seq->n1, 
		       bestp_arr[i]->rst.comp, m_msg.pstat_void);
	}
#ifndef SAMP_STATS
      }
      else pst.zsflag = -1;
#endif
    }

    /* **************************************************************** */
    /* do shuffles if too few library sequences or for second estimate  */
    /* **************************************************************** */

    /* if there are not many scores, produce better statistics by shuffling */
    /* but only if statistics are enabled (7-July-2008) */
    if (pst.zsflag != 3 && pst.zsflag > -1 && 
	nbest > 0 && nbest < m_msg.shuff_max) {

      buf_shuf_seq(aa0, m_msg.n0, &aa1shuff_b, aa1save, maxn,
		   bestp_arr, nbest, &pst, &m_msg, &m_bufi
#if !defined(COMP_THR) && !defined(PCOMPLIB)
		   , f_str
#endif
		   , &m_msg.ss_info);

      /* (4) analyze rstats */
      if (pst.zsflag < 10) pst.zsflag += 10;
      if (pst.zsflag > 20) pst.zsflag -= 10;

      pst.zsflag_f = process_hist(rstats,nrstats,&m_msg, &pst,&m_msg.hist,
				  &m_msg.pstat_void,&m_msg.ss_info,0);
    }

    /* **************************************************************** */
    /* done with shuffling for small sample size                        */
    /* **************************************************************** */

    if (!pst.zdb_size_set) pst.zdb_size = m_msg.ldb.entries;

#if defined(COMP_THR) || defined(PCOMPLIB)
    /* before I call last_calc/showbest/showalign, I need init_work() to
       get an f_str. This duplicates some code above, which is used in
       the non-threaded version.

       I have tried to get an f_str from one of the threads, but on
       some architectures, that f_str is not available to the main thread.
    */

    if (!have_f_str) {
      init_work(aa0[0],m_msg.n0,&pst,&f_str[0]);
      /* f_str[0] = work_info[0].f_str_ap[0]; */
      have_f_str = 1;
      f_str[5] = f_str[4] = f_str[3] = f_str[2] = f_str[1] =  f_str[0];

      if (m_msg.qframe == 2) {
	if ((aa0[1]=(unsigned char *)calloc((size_t)m_msg.n0+2+SEQ_PAD,
					    sizeof(unsigned char)))==NULL) {
	  fprintf(stderr," cannot allocate aa0[1][%d] for alignments\n",
		  m_msg.n0+2+SEQ_PAD);
	}
	*aa0[1]='\0';
	aa0[1]++;
	memcpy(aa0[1],aa0[0],m_msg.n0+1);
	/* for ALTIVEC/SSE2, must pad with 16 NULL's, but not necessary after calloc() */
	for (id=0; id<SEQ_PAD; id++) {aa0[1][m_msg.n0+id]=0;}

	revcomp(aa0[1],m_msg.n0,&pst.c_nt[0]);
	init_work(aa0[1],m_msg.n0,&pst,&f_str[1]);
	/* f_str[1] = work_info[0].f_str_ap[1]; */
      }
    }
#endif

    /* now we have one set of scaled scores for in bestp_arr  -
       for FASTS/F, we need to do some additional processing */

    if (!m_msg.qshuffle) {
      last_stats(aa0[0], m_msg.n0, stats,nstats, bestp_arr,nbest,
		 &m_msg, &pst, &m_msg.hist, &m_msg.pstat_void);
    }
    else {
      last_stats(aa0[0], m_msg.n0,
		 qstats,nqstats, bestp_arr,nbest, &m_msg, &pst, 
		 &m_msg.hist, &m_msg.pstat_void);
    }

    /* here is a contradiction: if pst.zsflag < 0, then m_msg.pstat_void
       should be NULL; if it is not, then process_hist() has been called */
    if (pst.zsflag < 0 && m_msg.pstat_void != NULL) pst.zsflag = 1;

    if (m_msg.last_calc_flg) {
      /* last_calc may need coefficients from last_stats() */
      nbest = last_calc(aa0, aa1save, maxn, bestp_arr, nbest, &m_msg, &pst,
			f_str, m_msg.pstat_void);
    }

    /* in addition to scaling scores, this sorts bestp_arr[nbest] */
    scale_scores(bestp_arr,nbest,m_msg.db, &pst,m_msg.pstat_void);

    /* For large databases, we have good zscores for all the MAX_BEST
       sequences.  Thus, we can sort the scores, and get a list of all
       sequences with scores better than the E() threshold.  If zsflag
       > 20, then we should shuffle those guys to get an alternative
       estimate of lambda and K */
    
    if (pst.zsflag > 20 && nbest >= m_msg.shuff_max) {
      n_sig = nbest;
      for (i=0; i<nbest; i++) {
	if (bestp_arr[i]->rst.escore > pst.e_cut) {
	  n_sig = i;
	  break;
	}
      }

      /* if there are no significant hits, shuffle the top 10 */
      if (n_sig < 10) n_sig = 10;

      /* check to see how many significant sequences there are, and
	 ensure that every sequence is shuffled at least 5 times */
      if (n_sig * 5 > m_msg.shuff_max) {
	m_msg.shuff_max = n_sig*5;
	if (m_msg.shuff_max > MAX_BEST/2) m_msg.shuff_max = MAX_BEST/2;
	if ((rstats = (struct stat_str *)realloc(rstats, m_msg.shuff_max * sizeof(struct stat_str)))==NULL) {
	  fprintf(stderr, " *** Cannot reallocate rstats[%d] ***\n",m_msg.shuff_max);
	  exit(1);
	}
      }

      buf_shuf_seq(aa0, m_msg.n0, &aa1shuff_b, aa1save, maxn,
		   bestp_arr, n_sig, &pst, &m_msg, &m_bufi
#if !defined(COMP_THR) && !defined(PCOMPLIB)
		   , f_str
#endif
		   , &m_msg.ss_info);

      zs_off_save = pst.zs_off;
      /* ensure that hist2.hist_a is initialized properly */
      hist2.hist_a = NULL;
      pst.zsflag_f = process_hist(rstats,nrstats,&m_msg, &pst,&hist2,
				  &m_msg.pstat_void2,&m_msg.ss_info, 0);
      pst.zs_off = zs_off_save;

      for (i=0; i<nbest; i++) {
	bestp_arr[i]->zscore2 =
	  find_z(bestp_arr[i]->rst.score[pst.score_ix],
		     bestp_arr[i]->rst.escore,
		     bestp_arr[i]->seq->n1,bestp_arr[i]->rst.comp,m_msg.pstat_void2);
      }
    }

    get_param(&pst, info_gstring2p,info_gstring3, &m_msg.s_info);

    /* **************************************************************** */
    /* label Library: output                                            */
    /* **************************************************************** */

    tscan = s_time();
    if (m_msg.std_output) {
      printf("%s\n",info_lib_range_p);

      if (m_msg.db.carry==0) {
	fprintf(stdout, "  %7ld residues in %5ld sequences\n", m_msg.db.length, m_msg.db.entries);
      }
      else {
	db_tt = (double)m_msg.db.carry*(double)LONG_MAX + (double)m_msg.db.length;
	fprintf(stdout, "  %.0f residues in %5ld library sequences\n", db_tt, m_msg.db.entries);
      }

      prhist (stdout, &m_msg, &pst, m_msg.hist, nstats, sstats, m_msg.ldb,
	      (pst.zsflag > 20? hist2.stat_info:NULL),info_lib_range_p,
	      info_gstring2p, info_hstring_p);

      printf (" Scan time: ");
      ptime(stdout,tscan-tprev);
      printf ("\n");
    }
#ifdef COMP_MLIB
    ttscan += tscan-tprev;
#endif

    /* fprintf(stderr, " *** shuff_tot: %d ***\n",shuff_tot); */

  l3:
    if (!m_msg.quiet) {
      printf("Enter filename for results [%s]: ", m_msg.outfile);
      fflush(stdout);
    }

    rline[0]='\0';
    if (!m_msg.quiet && fgets(rline,sizeof(rline),stdin)==NULL) goto end_l;
    if ((bp=strchr(rline,'\n'))!=NULL) *bp = '\0';
    if (rline[0]!='\0') strncpy(m_msg.outfile,rline,sizeof(m_msg.outfile));
    if (m_msg.outfile[0]!='\0') {
      if (!m_msg.outfd) {
	if ((outfd=fopen(m_msg.outfile,"w"))==NULL) {
	  fprintf(stderr," could not open %s\n",m_msg.outfile);
	  if (!m_msg.quiet) goto l3;
	  else goto l4;
	}
	m_msg.outfd = outfd;
      }
      else {
	outfd = m_msg.outfd;
      }

#ifdef PGM_DOC
      fputs(argv_line,outfd);
      fputc('\n',outfd);
#endif  
      fputs(iprompt0,outfd);
      fprintf(outfd," %s%s\n",mp_verstr,refstr);

      fprintf(outfd,"Query: %s%s, %d %s\n",
	      info_qlabel, (m_msg.revcomp ? "-" : "\0"), m_msg.n0,
	      m_msg.sqnam);

#ifdef COMP_MLIB
      fprintf(outfd,"%3d>>>%s - %d %s%s\n", qlib,
	   m_msg.qtitle, m_msg.n0, m_msg.sqnam,
	   (m_msg.revcomp ? " (reverse complement)" : rline));
#else
      fprintf(outfd,"%.50s: %d %s%s\n %s\n",
	   info_qlabel, m_msg.n0, m_msg.sqnam,
	   (m_msg.revcomp ? " (reverse complement)" : rline));
#endif
      fprintf(outfd, "Library: %s%s\n", m_msg.ltitle,info_lib_range_p);
      if (m_msg.db.carry==0) {
	fprintf(outfd, "  %7ld residues in %5ld sequences\n", m_msg.db.length, m_msg.db.entries);
      }
      else {
	db_tt = (double)m_msg.db.carry*(double)LONG_MAX + (double)m_msg.db.length;
	fprintf(outfd, "  %.0f residues in %5ld library sequences\n", db_tt, m_msg.db.entries);
      }

      prhist(outfd, &m_msg, &pst,m_msg.hist, nstats, sstats, m_msg.db,
	     (pst.zsflag > 20? hist2.stat_info:NULL), info_lib_range_p,
	     info_gstring2p, info_hstring_p);

      fprintf (outfd, " Scan time: ");
      ptime(outfd,tscan-tprev);
      fprintf (outfd,"\n");
    }

  l4:   
    if (m_msg.markx & MX_HTML) {
      fputs("</pre>\n<p>\n<hr>\n<p>\n",outfd);
    }

    /* code from p2_complib.c to pre-calculate -m 9 alignment info -
       requires -q with -m 9 */

    /* find the lowest scoring alignment to be displayed */
    /* the logic has been modified to preserve m_msg.mshow, so that
       m_msg.nshow always provides the number of alignments to be
       displayed in quiet mode */

    /* skip entries if -F e_low specified */
    if (pst.zsflag >= 0) {
      for (i=0; i<nbest && bestp_arr[i]->rst.escore < m_msg.e_low; i++) {};
      m_msg.nskip = i;
      if (pst.do_rep) {
	for (i=m_msg.nskip; i < nbest && bestp_arr[i]->rst.escore < m_msg.e_cut; i++) {
	  bestp_arr[i]->repeat_thresh = 
	    min(E1_to_s(pst.e_cut_r, m_msg.n0, bestp_arr[i]->seq->n1,
			pst.zdb_size, m_msg.pstat_void),bestp_arr[i]->rst.score[pst.score_ix]);
	}
      }
    }
    else {
      m_msg.nskip = 0;
      for (i=0; i < nbest; i++) {
	bestp_arr[i]->repeat_thresh = bestp_arr[i]->rst.score[pst.score_ix];
      }
    }

    if (m_msg.quiet || m_msg.markx & MX_M9SUMM) {

      /* to determine how many sequences to re-align (either for
	 do_opt() or calc_id() we need to modify m_msg.mshow to get
	 the correct number of alignments */

      if (pst.zsflag >= 0) {	/* do we have e_values? */
	for (i=m_msg.nskip; i<nbest && bestp_arr[i]->rst.escore < m_msg.e_cut; i++) {}

	if (m_msg.mshow_flg != 1) { m_msg.nshow = min(i - m_msg.nskip, nbest-m_msg.nskip);}
	else { m_msg.nshow = min(m_msg.mshow, i-m_msg.nskip);}
      }

      if (m_msg.nshow <= 0) { /* no results to display */
	if (m_msg.std_output) fprintf(outfd,"!! No sequences with E() < %0.5g\n",m_msg.e_cut);
	m_msg.nshow = 0;
	goto end_l;
      }
    }

    if (m_msg.quiet && ((m_msg.stages > 1) || (m_msg.markx & MX_M9SUMM))) {
      /* pre-load sequence data for  alignments for showbest, showalign */
      pre_load_best(aa1save, maxn, &bestp_arr[m_msg.nskip], m_msg.nshow, &m_msg);


      buf_align_seq(aa0, m_msg.n0, &bestp_arr[m_msg.nskip], m_msg.nshow - m_msg.nskip,
		    &pst, &m_msg, &m_bufi
#if !defined(COMP_THR) && !defined(PCOMPLIB)
		    , f_str
#endif
		    );
    }

    if (m_msg.do_showbest) {
      showbest (stdout, aa0, aa1save, maxn, &bestp_arr[m_msg.nskip], nbest-m_msg.nskip,
		qtt.entries, &m_msg, &pst,m_msg.db, info_gstring2p, f_str);

      if (outfd != stdout) {
	t_quiet = m_msg.quiet;
	m_msg.quiet = -1;	/* should guarantee 1..nbest shown */
	showbest (outfd, aa0, aa1save, maxn, &bestp_arr[m_msg.nskip], nbest-m_msg.nskip,
		  qtt.entries, &m_msg, &pst, m_msg.db, info_gstring2p, f_str);
	m_msg.quiet = t_quiet;
      }
    }

    /* m_msg.ashow can be -1 or > 0 to show results */
    if (m_msg.nshow > 0 && m_msg.ashow != 0) {
      rline[0]='N';
      if (!m_msg.quiet){
	printf(" Display alignments also? (y/n) [n] "); fflush(stdout);
	if (fgets(rline,sizeof(rline),stdin)==NULL) goto end_l;
      }
      else rline[0]='Y';

      if (toupper((int)rline[0])=='Y') {
	if (!m_msg.quiet && m_msg.do_showbest) {
	  printf(" number of alignments [%d]? ",m_msg.nshow);
	  fflush(stdout);
	  if (fgets(rline,sizeof(rline),stdin)==NULL) goto end_l;
	  if (rline[0]!=0) sscanf(rline,"%d",&utmp);
	  if (utmp == 0) utmp = -1;
	  m_msg.ashow=min(utmp,m_msg.nshow);
	}

	if (m_msg.std_output && m_msg.markx & (MX_AMAP+ MX_HTML + MX_M9SUMM) && !(m_msg.markx & MX_M10FORM)) {
	  fprintf(outfd,"\n>>>%s%s, %d %s vs %s library\n",
		  info_qlabel,(m_msg.revcomp ? "_rev":"\0"), m_msg.n0,
		  m_msg.sqnam,m_msg.lname);
	}

	if (m_msg.markx & MX_M10FORM) {
	  fprintf(outfd,"\n>>>%s%s, %d %s vs %s library\n",
		  info_qlabel,(m_msg.revcomp ? "-":"\0"), m_msg.n0, m_msg.sqnam,
		  m_msg.lname);
	  fprintf(outfd,"; pg_name: %s\n",argv[0]);
	  fprintf(outfd,"; pg_ver: %s\n",mp_verstr);
	  fprintf(outfd,"; pg_argv:");
	  for (i=0; i<argc; i++)
	    fprintf(outfd," %s",argv[i]);
	  fputc('\n',outfd);
	  fputs(info_gstring3,outfd);
	  fputs(info_hstring_p[0],outfd);
	  fputs(info_hstring_p[1],outfd);
	}

	showalign (outfd, aa0, aa1save, maxn, &bestp_arr[m_msg.nskip], nbest-m_msg.nskip,
		   qtt.entries, &m_msg, &pst, info_gstring2p, f_str, &m_bufi);

	fflush(outfd);
      }
    }

  end_l:

    if (!(m_msg.markx & MX_M8OUT) && (m_msg.markx & (MX_M10FORM+MX_M9SUMM))) {
      fprintf(outfd,">>><<<\n");
    }

    if (fdata) {
      fprintf(fdata,"/** Algorithm : %s  **/\n",info_gstring2p[0]);
      fprintf(fdata,"/** Parameters : %s  **/\n",info_gstring2p[1]);
      fprintf(fdata,"%3ld%-50s\n",qtt.entries-1,m_msg.qtitle);
      fflush(fdata);
    }
    
#if defined(COMP_THR) || defined(PCOMPLIB)
    rbuf_done(fa_max_workers);
#endif

#if defined(COMP_THR) || defined(PCOMPLIB)
    for (i=0; i<m_bufi.max_work_buf; i++) {
      lib_buf2_list[i].hdr.buf2_cnt=
	lib_buf2_list[i].hdr.have_results=
	lib_buf2_list[i].hdr.have_best_save = 
	lib_buf2_list[i].hdr.aa1b_used = 0;
    }

    num_reader_bufs = m_bufi.max_work_buf;
#endif

#if defined(COMP_THR)
    num_worker_bufs = 0;
    reader_done = 0;
    reader_wait = 1;
    worker_buf_workp = 0;
    worker_buf_readp = 0;
    reader_buf_workp = 0;
    reader_buf_readp = 0;

    start_thread = 1;	/* stop thread from starting again */
#endif


    /* clean up best_seqs */
    memset(best_seqs,0,(MAX_BEST+1)*sizeof(struct seq_record));

    /* re-initialize lib_buf2_list buffers */
    for (lib_bhead_p = lib_buf2_list; 
	 lib_bhead_p < lib_buf2_list+m_bufi.max_work_buf; ) {

      /* this wipes out lib_bhead_p->hdr.buf2[0].seq, .mseq */
      memset(lib_bhead_p->buf2_data,0,(size_t)(m_bufi.max_buf2_cnt+1)*sizeof(struct buf2_data_s));
      /* replace it */
      lib_bhead_p++->hdr.have_results = 0;
    }

    /* re-initialize library counts */
    m_msg.ldb.length = 0l;
    m_msg.ldb.entries = m_msg.ldb.carry = 0;

    /* clean up alignment encodings */
    /* needs to deallocate aln_code, ann_code */
    /* bestp_arr[i]->a_res is a pointer, so it must be free()'d */
    for (i=m_msg.nskip; i < m_msg.nskip+m_msg.nshow; i++) {
      if (bestp_arr[i]->have_ares & 0x2) {
	cur_ares_p = bestp_arr[i]->a_res;
	while (cur_ares_p) {
	  if (cur_ares_p->aln_code) free(cur_ares_p->aln_code);
	  if (cur_ares_p->ann_code) free(cur_ares_p->ann_code);
	  if (bestp_arr[i]->have_ares & 0x1 && cur_ares_p->res) free(cur_ares_p->res);
	  next_ares_p = cur_ares_p->next;
	  free(cur_ares_p);
	  cur_ares_p = next_ares_p;
	}
	bestp_arr[i]->a_res = NULL;
      }
      bestp_arr[i]->have_ares = 0;
      bestp_arr[i]->mseq->bline = NULL;
      bestp_arr[i]->mseq->bline_max = 0;
    }

#ifdef DEBUG
    /* check to see if there are ANY un-reset have_ares */
    for (i=0; i< nbest; i++) {
      if (bestp_arr[i]->have_ares) {
	fprintf(stderr," Un-reset have_ares[%d]: %d\n",i,bestp_arr[i]->have_ares);
	bestp_arr[i]->have_ares = 0;
      }
    }
#endif

    if (m_msg.qframe == 2) free(aa0[1]-1);

    if (have_f_str) {
      have_f_str = 0;
      if (f_str[1]!=f_str[0]) {
	close_work (aa0[1], m_msg.n0, &pst, &f_str[1]);
      }
      close_work (aa0[0], m_msg.n0, &pst, &f_str[0]);
    }
#if !defined(COMP_THR) && !defined(PCOMPLIB)
    if (m_msg.qshuffle) close_work (aa0s, m_msg.n0, &pst, &qf_str);
#endif
    if (pst.pam_pssm) {
      free_pam2p(pst.pam2p[0]);
      free_pam2p(pst.pam2p[1]);
    }

    if (aa1shuff_b != NULL) {
      free(aa1shuff_b);
      aa1shuff_b = NULL;
    }

    if (m_msg.aa1save_buf_b != NULL) {
      free(m_msg.aa1save_buf_b);
      m_msg.aa1save_buf_b = NULL;
    }

    if (m_msg.bline_buf_b != NULL) {
      free(m_msg.bline_buf_b);
      m_msg.bline_buf_b = NULL;
    }

    for (cur_lib_p=lib_list_p; cur_lib_p != NULL; cur_lib_p = cur_lib_p->next) {
      if (cur_lib_p->m_file_p !=NULL) {
	closelib(cur_lib_p->m_file_p);
      }
    }

    tddone = time(NULL);
    tdone = s_time();
    fflush(outfd);

    ttdisp += tdone-tscan;

    /* reset pst parameters to original */
    pst.zsflag = m_msg.zsflag;
    pst.zsflag2 = m_msg.zsflag2;
    pst.n1_low = m_msg.n1_low;
    pst.n1_high = m_msg.n1_high;

    /*     maxn = m_msg.max_tot; */

    m_msg.q_offset = next_q_offset;

    m_msg.n0 = 
      QGETLIB (aa0[0], MAXTST, m_msg.qtitle, sizeof(m_msg.qtitle),
	       &qseek, &qlcont,q_file_p,&m_msg.q_off);
    if (m_msg.n0 <= 0) break;
    if ((bp=strchr(m_msg.qtitle,' '))!=NULL) *bp='\0';
    strncpy(info_qlabel, m_msg.qtitle,sizeof(info_qlabel));
    if (bp != NULL) *bp=' ';
    info_qlabel[sizeof(info_qlabel)-1]='\0';


    if (m_msg.ann_flg) {
      m_msg.n0 = ann_scan(aa0[0],m_msg.n0,&m_msg.aa0a,m_msg.qdnaseq);
    }

    if (m_msg.ldb_info.term_code && m_msg.qdnaseq==SEQT_PROT &&
	aa0[0][m_msg.n0-1]!=m_msg.ldb_info.term_code) {
      aa0[0][m_msg.n0++]=m_msg.ldb_info.term_code;
      aa0[0][m_msg.n0]=0;
    }
    if (m_msg.outfd) {fputc('\n',stdout);}

    /* reset the seqr_chain for next query */
    reset_seqr_chain(lib_seqr_chain);
  }	/* end of while(1) for multiple queries */

#ifdef PCOMPLIB
  /* tell workers to quit */
  init_thr(fa_max_workers, NULL, NULL, NULL, NULL, NULL);
#endif

  if (m_msg.markx & MX_M8COMMENT) {
    printf("# %s processed %d queries\n",prog_func,qlib);
  }

  if ( !(m_msg.markx & MX_M8OUT) && (m_msg.markx & (MX_M10FORM+MX_M9SUMM))) {
    fprintf(outfd,">>>///\n");
  }

  tdone = s_time();
  if ( m_msg.markx & MX_HTML) fputs("<p><pre>\n",outfd); 
  if (m_msg.std_output) {
    printsum(outfd, m_msg.db);
  }
  if ( m_msg.markx & MX_HTML) fputs("</pre>\n",outfd);
#ifdef HTML_HEAD
  if (m_msg.markx & MX_HTML) fprintf(outfd,"</body>\n</html>\n");
#endif
  if (m_msg.std_output && outfd!=stdout) printsum(stdout,m_msg.db);

#ifdef PCOMPLIB
#ifdef MPI_SRC
  MPI_Finalize();
#endif
#endif
  exit(0);
}
/* **************************************************************** */
/* end of main() program                                            */
/* **************************************************************** */

void
printsum(FILE *fd, struct db_str ntt)
{
  double db_tt;
  char tstr1[26], tstr2[26];

  strncpy(tstr1,ctime(&tdstart),sizeof(tstr1));
  strncpy(tstr2,ctime(&tddone),sizeof(tstr1));
  tstr1[24]=tstr2[24]='\0';

  /* Print timing to output file as well */
  fprintf(fd, "\n\n%ld residues in %ld query   sequences\n", qtt.length, qtt.entries);
  if (ntt.carry == 0) 
    fprintf(fd, "%ld residues in %ld library sequences\n", ntt.length, ntt.entries);
  else {
    db_tt = (double)ntt.carry*(double)LONG_MAX + (double)ntt.length;
    fprintf(fd, "%.0f residues in %ld library sequences\n", db_tt, ntt.entries);
  }

#if !defined(COMP_THR) && !defined(PCOMPLIB)
  fprintf(fd," Scomplib [%s]\n start: %s done: %s\n",mp_verstr,tstr1,tstr2);
#endif
#if defined(COMP_THR)
  fprintf(fd," Tcomplib [%s] (%d proc)\n start: %s done: %s\n", mp_verstr,
    fa_max_workers,tstr1,tstr2);
#endif
#if defined(PCOMPLIB)
  fprintf(fd," Pcomplib [%s] (%d proc)\n start: %s done: %s\n", mp_verstr,
    fa_max_workers,tstr1,tstr2);
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
  fprintf (fd, "\nFunction used was %s [%s]\n", prog_func,mp_verstr);
}

void fsigint()
{
  struct db_str db;

  db.entries = db.length = db.carry = 0;
  tdone = s_time();
  tddone = time(NULL);

  printf(" /*** interrupted ***/\n");
  if (outfd!=stdout) fprintf(outfd,"/*** interrupted ***/\n");
  fprintf(stderr,"/*** interrupted ***/\n");

  printsum(stdout,db);
  if (outfd!=stdout) printsum(outfd,db);

  exit(1);
}

/* save the sequence meta-data for a sequence if we need to re-use the
   rbuf/wbuf buffer */

/* 20-May-2010 -- I think preserve seq is called far too often,
   essentially after every save_best() that adds results to beststr.
   It should only be called if the lib_buf2_dp is about to be reused.
   A better strategy would be to have a large set of seq_record[], and
   only call preserve seq when they are used up.
*/

void preserve_seq(struct buf2_data_s *lib_buf2_dp,
		  struct seq_record *best_seqs,
		  struct mseq_record *best_mseqs,
		  struct beststr *best) {
  struct seq_record *dest_seq_p, *saved_seq_p;
  struct mseq_record *dest_mseq_p, *saved_mseq_p;
  struct beststr *next_bbp;

  /* preserve seq is only called when lib_buf2_dp->best_save != NULL,
     so there is always something to save */
  saved_seq_p = lib_buf2_dp->best_save->seq;
  saved_mseq_p = lib_buf2_dp->best_save->mseq;

  dest_seq_p = &best_seqs[lib_buf2_dp->best_save - best];
  dest_mseq_p = &best_mseqs[lib_buf2_dp->best_save - best];

  lib_buf2_dp->best_save->seq = dest_seq_p;
  lib_buf2_dp->best_save->mseq = dest_mseq_p;

  for (next_bbp = lib_buf2_dp->best_save->bbp_link;
       (next_bbp != NULL) && (next_bbp->seq == saved_seq_p)
	 && (next_bbp->n1 == saved_seq_p->n1);
       next_bbp = next_bbp->bbp_link) {
    next_bbp->seq = dest_seq_p;
    next_bbp->mseq = dest_mseq_p;
  }

  memcpy(dest_seq_p,lib_buf2_dp->seq,sizeof(struct seq_record));
  memcpy(dest_mseq_p,lib_buf2_dp->mseq,sizeof(struct mseq_record));
  dest_seq_p->aa1b = NULL;
}

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

/* save_best captures much of the complexity of saving the best scores
   and appropriately sampling the scores for statistical analysis. It
   does the following:

   (1) update s_info counts for functions like fasta/x/y that don't
       optimize every score 

   (2) for every result in the buffer:
       (a) decide if it should be used for statistical sampling
       (b) if the number of samples > MAX_STATS, then run
           process_hist() and update all the zscores
       (c) if the zscore > zbestcut, save the results (and
           preserve_seq)
       (d) reset everything for next sequence

*/

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
	      rbuf_dp->mseq->libstr,
	      rbuf_dp->seq->n1,rbuf_dp->frame,rbuf_rp->rst.comp,rbuf_rp->rst.H,
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
	    stats_done = 1;
	    for (i=0; i< nstats; i++) {
	      bestp_arr[i]->zscore = 
		find_z(bestp_arr[i]->rst.score[ppst->score_ix],
			   bestp_arr[i]->rst.escore, bestp_arr[i]->seq->n1,
			   bestp_arr[i]->rst.comp, *pstat_voidp);
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

      if (t_valid_stat && stats_done) {
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
/*
      if (rbuf_dp->best_save) {	
	if (rbuf_dp->best_save->seq == rbuf_dp->seq) {
	  bbp->bbp_link = rbuf_dp->best_save;
	}
	else {
	  bbp->bbp_link = NULL;
	}
      }
      rbuf_dp->best_save = bbp;
*/
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
      lib_buf2_dp++;
      lib_buf2_rp++;
      continue;
    }

#ifdef DEBUG
    if (lib_buf2_dp->seq->aa1b[0] == '\0') {
      fprintf(stderr,"[%s/buf_do_work] -- invalid aa1 NULL at: %s frame: %d\n",
	      prog_func,lib_buf2_dp->mseq->libstr, lib_buf2_dp->frame);
    }

    if (check_seq_range(lib_buf2_dp->seq->aa1b, lib_buf2_dp->seq->n1,
			(ppst->ext_sq_set ? ppst->nsqx: ppst->nsq), "buf_do_work()")) {
      fprintf(stderr, "[%s/buf_do_work] range error at: %d/%d (n1:%d)\n",
	      prog_func,lib_bhead_p->hdr.buf2_cnt - (buf2_cnt+1),
	      lib_bhead_p->hdr.buf2_cnt, lib_buf2_dp->seq->n1);
    };

    /* also check for adler32_crc match */
    if (lib_buf2_dp->seq->adler32_crc != (atmp=adler32(1L,lib_buf2_dp->seq->aa1b, lib_buf2_dp->seq->n1))) {
      fprintf(stderr, "[%s/buf_do_work] CRC error [%lu!=%lu] at: %d/%d (n1:%d/l_offset:%ld)\n",
	      prog_func,lib_buf2_dp->seq->adler32_crc, atmp,
	      lib_bhead_p->hdr.buf2_cnt - (buf2_cnt+1),
	      lib_bhead_p->hdr.buf2_cnt,lib_buf2_dp->seq->n1,
	      lib_buf2_dp->seq->l_offset);
    }
#endif

    do_work (aa0[lib_buf2_dp->frame], n0,
	     lib_buf2_dp->seq->aa1b, lib_buf2_dp->seq->n1,
	     lib_buf2_dp->frame, ppst, f_str[lib_buf2_dp->frame], 0, 0,
	     &(lib_buf2_rp->rst), &(lib_bhead_p->s_cnt_info));

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

void
buf_shuf_work(unsigned char **aa0,  int n0,
	      unsigned char *aa1s,
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

    if ((lib_buf2_dp->stats_idx < 0) || lib_buf2_dp->seq->n1 < ppst->n1_low ||
	lib_buf2_dp->seq->n1 > ppst->n1_high ) {
      lib_buf2_dp++;
      lib_buf2_rp++;
      continue;
    }

    shuff_cnt++;
    if (ppst->zs_win > 0) 
      wshuffle(lib_buf2_dp->seq->aa1b,aa1s,lib_buf2_dp->seq->n1,ppst->zs_win, rand_state);
    else {
      if (ppst->shuffle_dna3) {	shuffle3(lib_buf2_dp->seq->aa1b,aa1s,lib_buf2_dp->seq->n1, rand_state);}
      else { shuffle(lib_buf2_dp->seq->aa1b,aa1s,lib_buf2_dp->seq->n1, rand_state);}
    }

    /* rshuffle(lib_buf2_dp->seq->aa1b,aa1s,lib_buf2_dp->seq->n1); */

#ifdef DEBUG
    if (check_seq_range(aa1s, lib_buf2_dp->seq->n1,
			(ppst->ext_sq_set ? ppst->nsqx: ppst->nsq), "buf_do_align()")) {
      fprintf(stderr, "[%s/buf_do_shuff] range error at: %d/%d (n1:%d)\n",
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
  struct beststr *bbp;
  char l_bline[MAX_SSTR];
  int n1lib_req, shuff_mult;
  long loffset, l_off;
  int n1, itt;
  int max_do_cnt;
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
      n1lib_req = 0;
      for (i = 0; i < nbest; i++) {
	n1lib_req += bestp_arr[i]->n1+ 2;
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

#else
      if (n1lib_req < 2) {
	fprintf(stderr,"[%s/buf_shuf_seq] no residues to shuffle: %d (%d)\n", 
		prog_func,n1lib_req,nbest);
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

      shuff_mult = (m_msp->shuff_max / nbest) + 1;
      istats = 0;
      
      /* setup lib_bhead buffers for shuffle comparisons */
#if defined(COMP_THR) || defined(PCOMPLIB)	/* threaded/parallel */
      /* max_do_cnt can be smaller than max_buf2_cnt, but not larger */
      max_do_cnt = min(m_bufi_p->max_buf2_cnt,
		       m_msp->shuff_max / (2 * fa_max_workers));
      /* we don't have a left over one, so we need one */
      get_rbuf(&lib_bhead_p,m_bufi_p->max_work_buf);
#else	/* not threaded */
      max_do_cnt = m_bufi_p->max_buf2_cnt;
      lib_bhead_p = lib_buf2_list;  /* equivalent to un-threaded get_rbuf() */
#endif
      max_buf_size = lib_bhead_p->hdr.aa1b_size;
      cur_buf_size = 0;
      lib_bhead_p->hdr.buf2_cnt = 0;
      lib_bhead_p->hdr.have_results = 0;
      lib_bhead_p->hdr.stop_work = 0;
      lib_bhead_p->hdr.buf2_type=BUF2_DOSHUF;
      lib_bhead_p->hdr.seq_record_continuous = 0;
      lib_buf2_dp = lib_bhead_p->buf2_data;
      lib_buf2_rp = lib_bhead_p->buf2_res;

      /* read sequences into shuffle buffer */

      for (i = 0; i < nbest; i++) {
	bbp = bestp_arr[i];
	if (bbp->seq->aa1b == NULL) {
	  /* get the sequence */
	  (bbp->mseq->m_file_p->ranlib)(l_bline, sizeof(l_bline),
				       bbp->mseq->lseek,bbp->mseq->libstr,bbp->mseq->m_file_p);
	  n1 = re_getlib(aa1save,NULL, maxn,m_msp->ldb_info.maxt3,
			 m_msp->ldb_info.l_overlap,bbp->mseq->cont,m_msp->ldb_info.term_code,
			 &loffset,&l_off,bbp->mseq->m_file_p);

	  /* fprintf(stderr, " %d gets %d %d\n",i,bestp_arr[i]->seq->n1,n1); */

	  memcpy(aa1shuff, aa1save, n1+1);
	  bbp->seq->aa1b = aa1shuff;
	  aa1shuff += n1 + 1;
	}

	for (j = 0; j < shuff_mult; j++ ) {
	  for (itt = m_msp->revcomp; itt <= m_msp->nitt1; itt++) {
	    /* this invalidates lib_buf2_p->seq */
	    lib_buf2_dp->seq = bbp->seq;
	    cur_buf_size += bbp->seq->n1+1;
	    lib_buf2_dp->stats_idx = istats++;
	    lib_buf2_dp->frame = itt;
	    lib_buf2_dp++;		/* point to next buf2 */
	    lib_buf2_rp++;		/* point to next buf2 */
	    lib_bhead_p->hdr.buf2_cnt++;

	    if (lib_bhead_p->hdr.buf2_cnt >= max_do_cnt ||
		cur_buf_size >= max_buf_size - m_msp->ldb_info.maxn) {
/* (2) send sequences for shuffling */
#if defined(COMP_THR) || defined(PCOMPLIB)	/* threaded - fill and empty buffers */
	      /* provide empty buffer to workers */
	      lib_bhead_p->hdr.have_data = 1;
	      put_rbuf(lib_bhead_p,m_bufi_p->max_work_buf);
	      get_rbuf(&lib_bhead_p,m_bufi_p->max_work_buf);
#else		/* non-thread - just do the searches */
	      if (lib_bhead_p->hdr.buf2_type & BUF2_DOSHUF) {
		buf_shuf_work(aa0,m_msp->n0, aa1save,lib_bhead_p,
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
      }			/* done with bestp_arr[] */

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
  /* max_do_cnt can be smaller than max_buf2_cnt, but not larger */
#ifdef COMP_THR
  max_align_cnt = min(m_bufi_p->max_buf2_cnt,
		      nbest / (4 * fa_max_workers));
#else
  max_align_cnt = min(m_bufi_p->max_buf2_cnt, nbest / fa_max_workers);
#endif
  if (max_align_cnt < 1) max_align_cnt = 1;

  get_rbuf(&lib_bhead_p,m_bufi_p->max_work_buf);
#else	/* not threaded */
  max_align_cnt = m_bufi_p->max_buf2_cnt;
  lib_bhead_p = lib_buf2_list;  /* equivalent to un-threaded get_rbuf() */
#endif

  max_buf_size = lib_bhead_p->hdr.aa1b_size;
  lib_bhead_p->hdr.buf2_cnt = 0;
  lib_bhead_p->hdr.have_results = 0;
  lib_bhead_p->hdr.stop_work = 0;
  lib_bhead_p->hdr.buf2_type=BUF2_DOALIGN;
  lib_bhead_p->hdr.seq_record_continuous = 0;
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
      lib_bhead_p->hdr.have_data = 1;
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
      lib_bhead_p->hdr.seq_record_continuous = 0;
      lib_bhead_p->hdr.stop_work = 0;
      lib_buf2_dp = lib_bhead_p->buf2_data;
      lib_buf2_ap = lib_bhead_p->buf2_ares;
    }
  }			/* done with bestp_arr[] */

#if defined(COMP_THR) || defined(PCOMPLIB)	/* if COMP_THR - fill and empty buffers */
  /* check last buffers for any results */
  lib_bhead_p->hdr.have_data = 1;
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

struct buf_head *
alloc_comp_bufs (struct mng_thr *m_bufi_p, struct mngmsg *m_msp,
		 int ave_seq_len) {
  struct buf_head *lib_buf2_list;
  int i, buf_siz;

  /* first allocate space for buffer headers */
  if ((lib_buf2_list = 
       (struct buf_head *)calloc((size_t)(m_bufi_p->max_work_buf),
				 sizeof(struct buf_head))) == NULL) {
    fprintf(stderr," cannot allocate lib_buf2_list[%d]\n", m_bufi_p->max_work_buf);
    exit(1);
  }

#if defined(COMP_THR) || defined(PCOMPLIB)
  if ((worker_buf = 
       (struct buf_head **)calloc((size_t)(m_bufi_p->max_work_buf),
				  sizeof(struct buf_head *))) == NULL) {
    fprintf(stderr," cannot allocate **worker_buf[%d]\n", m_bufi_p->max_work_buf);
    exit(1);
  }
#endif

#ifdef COMP_THR
  if ((reader_buf = 
       (struct buf_head **)calloc((size_t)(m_bufi_p->max_work_buf),
				  sizeof(struct buf_head *))) == NULL) {
    fprintf(stderr," cannot allocate **reader_buf[%d]\n", m_bufi_p->max_work_buf);
    exit(1);
  }
#endif

  /* allocate space for library buffers and results */

  /* there are four structures/buffers used to keep track of
     sequences/results:

     (1) lib_buf2_list[] is a buf_head array, with:

     	buf2_hdr_s hdr - 
	   that stores whether the results are ready, number of
	   results available, and information about the seq_record
	   buffer and aa1b buffer

     (2) buf2_data[] -> seq_record/mseq_record arrays

     (3) buf2_res[] -> score results

     (4) buf2_ares[] -> alignment results
  */

  buf_siz = max(m_bufi_p->max_buf2_cnt*ave_seq_len, m_msp->max_tot * 4);
  if (buf_siz < m_msp->max_tot) buf_siz = m_msp->max_tot;
  m_bufi_p->seq_buf_size = buf_siz;

  for (i=0; i<m_bufi_p->max_work_buf; i++) {
    /* allocate max_buf2_cnt buf2_str's into each buf2 */

    if ((lib_buf2_list[i].buf2_ares =
	 (struct buf2_ares_s *)calloc((size_t)(m_bufi_p->max_buf2_cnt+1),
				   sizeof(struct buf2_ares_s)))
        ==NULL) {
      fprintf(stderr," cannot allocate buffer struct %d %d\n",
	      i,m_bufi_p->max_buf2_cnt+1);
      exit(1);
    }

    if ((lib_buf2_list[i].buf2_res =
	 (struct buf2_res_s *)calloc((size_t)(m_bufi_p->max_buf2_cnt+1),
				   sizeof(struct buf2_res_s)))
        ==NULL) {
      fprintf(stderr," cannot allocate buffer struct %d %d\n",
	      i,m_bufi_p->max_buf2_cnt+1);
      exit(1);
    }

    if ((lib_buf2_list[i].buf2_data =
	 (struct buf2_data_s *)calloc((size_t)(m_bufi_p->max_buf2_cnt+1),
				   sizeof(struct buf2_data_s)))
        ==NULL) {
      fprintf(stderr," cannot allocate buffer struct %d %d\n",
	      i,m_bufi_p->max_buf2_cnt+1);
      exit(1);
    }

    lib_buf2_list[i].hdr.seq_record_continuous = 0;	/* should be 1 */

    lib_buf2_list[i].hdr.aa1b_used = 0;
    lib_buf2_list[i].hdr.aa1b_size = buf_siz;

    /* make certain there is a '\0' at the beginning */

    lib_buf2_list[i].hdr.have_results=0;

#if defined(COMP_THR) || defined(PCOMPLIB)
    RESULTS_BUF[i] = &lib_buf2_list[i];
#endif
  }
  return lib_buf2_list;
}

struct seqr_chain *
new_seqr_chain(int max_seq_cnt, int aa1b_size, struct seqr_chain *old_seqr_chain, int maxn) {
  struct seqr_chain *my_seqr_chain;

  /* allocate the chain */
  if ((my_seqr_chain = (struct seqr_chain *)calloc(1,sizeof(struct seqr_chain)))==NULL) {
    fprintf(stderr," Cannot allocate library seqr_chain\n");
    exit(1);
  }
  if (old_seqr_chain) { old_seqr_chain->next = my_seqr_chain;}

  my_seqr_chain->cur_seqr_p = my_seqr_chain;
  my_seqr_chain->max_seq_cnt = max_seq_cnt;
  my_seqr_chain->cur_seq_cnt = 0;

  /* now allocate the seq_record, mseq_record buffers, and the space for the sequences */
  if ((my_seqr_chain->seqr_base = (struct seq_record *)
       calloc((size_t)max_seq_cnt,sizeof(struct seq_record))) == NULL) {
    fprintf(stderr," cannot allocate seq_record buffer %d\n",max_seq_cnt);
    exit(1);
  }

  if ((my_seqr_chain->mseqr_base = (struct mseq_record *)
       calloc((size_t)max_seq_cnt,sizeof(struct mseq_record))) == NULL) {
    fprintf(stderr," cannot allocate mseq_record buffer %d\n",max_seq_cnt);
    exit(1);
  }

  /* try to adjust aa1b_size to a sensible value based on past history */
  if (old_seqr_chain != NULL) {
    if (old_seqr_chain->aa1b_next+maxn < old_seqr_chain->aa1b_size / 2) {
      aa1b_size /= 2;
      aa1b_size = max(aa1b_size, maxn * 2);
    }
    else if (old_seqr_chain->aa1b_next > old_seqr_chain->aa1b_size - maxn) {
      aa1b_size *= 2;
    }
  }

  /* finally, we need a big sequence buffer */
  if ((my_seqr_chain->aa1b_base =
       (unsigned char *)calloc((size_t)(aa1b_size+1),sizeof(unsigned char)))
      ==NULL) {
    fprintf(stderr," cannot allocate buffer %d\n",aa1b_size+1);
    exit(1);
  }
  my_seqr_chain->aa1b_base++;
  my_seqr_chain->aa1b_size = aa1b_size;

  return my_seqr_chain;
}

/* make certain that max_seq_cnt reflects the number of sequences in
   seqr_chain if terminated by aa1b_size or end of library */
void end_seqr_chain(struct seqr_chain *last_seqr) {
  /* the last seqrecord_p is always one too far */
  last_seqr->cur_seq_cnt--;
  last_seqr->max_seq_cnt = last_seqr->cur_seq_cnt;
}

struct seq_record *
new_sequence_p(struct mseq_record **cur_mseq_p,  struct seq_record *prev_seq_p,
	       struct seqr_chain **cur_seqr_chain, int max_seq_cnt, int maxn) {

  struct seq_record *next_seq_p;

  /* update ->aa1b_next */
  if (prev_seq_p != NULL) {
#ifdef DEBUG
    if (prev_seq_p->n1 <= 0) {
      fprintf(stderr,"*** ERROR: prev_seq_p->n1 <=0: %d\n",prev_seq_p->n1);
    }
#endif
    (*cur_seqr_chain)->aa1b_next += prev_seq_p->n1 + 1;
  }

  /* is there room in cur_seqr_chain ?? */
  if ((*cur_seqr_chain)->cur_seq_cnt >= (*cur_seqr_chain)->max_seq_cnt
      || (*cur_seqr_chain)->aa1b_next + maxn >= (*cur_seqr_chain)->aa1b_size) {
    /* out of room in seqr_chain or aa1b_base */

    if ((*cur_seqr_chain)->aa1b_next + maxn >= (*cur_seqr_chain)->aa1b_size) {
      end_seqr_chain(*cur_seqr_chain);
    }

    *cur_seqr_chain = new_seqr_chain(max_seq_cnt, (*cur_seqr_chain)->aa1b_size,
				     *cur_seqr_chain, maxn);
    prev_seq_p = NULL;
  }

  /* is there room in aa1b_base ?? */
  if ((*cur_seqr_chain)->aa1b_next + maxn < (*cur_seqr_chain)->aa1b_size) {

    *cur_mseq_p = (*cur_seqr_chain)->mseqr_base + (*cur_seqr_chain)->cur_seq_cnt;
    next_seq_p = (*cur_seqr_chain)->seqr_base + (*cur_seqr_chain)->cur_seq_cnt;
    next_seq_p->aa1b = (*cur_seqr_chain)->aa1b_base + (*cur_seqr_chain)->aa1b_next;

#ifdef DEBUG
    if (next_seq_p->aa1b[-1]!='\0') {
      fprintf(stderr,"*** ERROR: new_sequence_p() missing NULL [-1]: %d\n",next_seq_p->aa1b[-1]);
    }
#endif

    (*cur_seqr_chain)->cur_seq_cnt++;
    return next_seq_p;
  }

  return NULL;
};


struct seqr_chain *
read_seqr_chain(struct mng_thr *m_bufi_p, struct lib_struct *lib_list_p,
		struct mngmsg *m_msp, struct pstruct *ppst, int fa_workers)
{

  struct seqr_chain *my_seqr_chain=NULL, *seqr_base;
  struct lib_struct *cur_lib_p;
  struct seq_record *current_seq_p=NULL;
  struct mseq_record *current_mseq_p;

  int lcont, ocont, maxt;	/* continued sequence */
  int igncnt=0;			/* count for ignoring sequences warning */
  struct lmf_str *m_file_p;
  unsigned char *aa1ptr, *aa1, *aa1save;
  char *bp;
  int n1;
  int i, id;

  char libstr[MAX_UID];		/* required because getlib() does not return lcont > 0 */
  int n_libstr=MAX_UID;			/* length of libstr */

  int *n1tot_ptr=NULL, *n1tot_cur;
  int n1tot_cnt=0;
  int n1tot_v;
  long loffset;

  if ((aa1save = (unsigned char *)calloc((m_msp->ldb_info.maxn+1),sizeof (char))) == NULL) {
    s_abort ("Unable to allocate library overlap", "");
  }
  *aa1save= '\0';
  aa1save++;

  seqr_base = my_seqr_chain = 
    new_seqr_chain(m_bufi_p->max_buf2_cnt*fa_workers,
		   (m_bufi_p->seq_buf_size+1)*fa_workers, NULL, m_msp->ldb_info.maxn);

  current_seq_p = new_sequence_p(&current_mseq_p, NULL, &my_seqr_chain, 
				 m_bufi_p->max_buf2_cnt*fa_workers, m_msp->ldb_info.maxn);

  /* scan through the libraries */
  for (cur_lib_p = lib_list_p; cur_lib_p; cur_lib_p=cur_lib_p->next) {

    if ((cur_lib_p->m_file_p = m_file_p=
	 open_lib(cur_lib_p, m_msp->ldb_info.ldnaseq, lascii, !m_msp->quiet))
	==NULL) {
      fprintf(stderr," cannot open library %s\n",cur_lib_p->file_name);
      continue;
    }

    loffset = 0l;
    lcont = 0;
    ocont = 0;
    n1tot_v = n1tot_cnt = 0;
    n1tot_cur = n1tot_ptr = NULL;

    /* get next buffer to read into */
    maxt = m_msp->ldb_info.maxn;

    aa1ptr = aa1 = current_seq_p->aa1b;

    while ((n1=GETLIB(aa1ptr, maxt, libstr, n_libstr,
		      &(current_mseq_p->lseek), &lcont, m_file_p,
		      &(current_seq_p->l_off)))>=0) {

#ifdef DEBUG
      /* check for out of range sequence */
      for (id=0; id<n1; id++) {
	if (aa1[id] > ppst->nsq_e) {
	  fprintf(stderr," *** ERROR *** %s[%d] = %d > %d out of range\n",libstr, id, aa1[id], ppst->nsq_e);
	  aa1[id] = 1;
	}
      }
#endif

      current_seq_p->index = seq_index;
      current_mseq_p->index = seq_index++;
      current_mseq_p->m_file_p = (void *)m_file_p;
      current_mseq_p->cont = ocont+1;
      current_seq_p->l_offset = loffset;

      if ((bp=strchr(libstr,' '))!=NULL) *bp='\0';
      strncpy(current_mseq_p->libstr,libstr,MAX_UID);	/* get old libstr for lcont>0 */

      /* add termination code for FASTX/FASTY if necessary */
      if (m_msp->ldb_info.term_code && !lcont &&
	  m_msp->ldb_info.ldnaseq==SEQT_PROT && aa1ptr[n1-1]!=m_msp->ldb_info.term_code) {
	aa1ptr[n1++]=m_msp->ldb_info.term_code;
	aa1ptr[n1]=0;
      }

      /* update n1 after possible changes */
      current_seq_p->n1 = n1;

#ifdef DEBUG
      if (n_libstr <= MAX_UID) {if ((bp=strchr(current_mseq_p->libstr,' '))!=NULL) *bp='\0';}
      if (aa1[-1]!='\0' || aa1ptr[n1]!='\0') {
	fprintf(stderr,"%s: aa1[%d] at %ld:%lld  missing NULL boundaries: %d %d\n",
		current_mseq_p->libstr,n1, m_msp->db.entries+1,current_mseq_p->lseek,
		aa1[-1],aa1ptr[n1]);
      }
#endif

      /* check for a continued sequence and provide a pointer to 
	 the n1_tot array if lcont || ocont */
      n1tot_v += n1;
      if (lcont && !ocont) {	/* get a new pointer */
	if (n1tot_cnt <= 0) {
	  if ((n1tot_ptr=calloc(1000,sizeof(int)))==NULL) {
	    fprintf(stderr," cannot allocate n1tot_ptr\n");
	    exit(1);
	  }
	  else {n1tot_cnt=1000;}
	}
	n1tot_cnt--;
	n1tot_cur = n1tot_ptr++;
      }
      current_mseq_p->n1tot_p = n1tot_cur;

      m_msp->db.entries++;
      m_msp->db.length += n1;
      if (m_msp->db.length > LONG_MAX) {
	m_msp->db.length -= LONG_MAX; m_msp->db.carry++;
      }

      /* don't count long sequences more than once */
      if (aa1!=aa1ptr) {	/* this is a continuation */
	current_seq_p->n1 = n1 += m_msp->ldb_info.l_overlap;	/* corrected 28-June-2008 */
	m_msp->db.entries--;
      }

      if ( n1 <= 1) {
	/*	if (igncnt++ <10)
		fprintf(stderr,"Ignoring: %s\n",current_seq_p->libstr);
	*/
	goto loop2;
      }

#ifdef DEBUG
      current_seq_p->adler32_crc =
	current_mseq_p->adler32_crc = adler32(1L,current_seq_p->aa1b,current_seq_p->n1);

      /* This finds most reasons for core dumps */
      if (ppst->debug_lib)
	for (i=0; i<n1; i++) {
	  if (aa1[i]>ppst->nsq || aa1[i] <= 0) {
	    fprintf(stderr,
		    "%s residue[%d/%d] %d range (%d) lcont/ocont: %d/%d\n%s\n",
		    current_mseq_p->libstr,i,current_seq_p->n1,aa1[i],ppst->nsq,lcont,ocont,aa1ptr+i);
	    aa1[i]=0;
	    n1=i-1;
	    break;
	  }
	}
#endif

      if (lcont) {
	memcpy(aa1save,&aa1[n1-m_msp->ldb_info.l_overlap],m_msp->ldb_info.l_overlap);
      }

      current_seq_p = new_sequence_p(&current_mseq_p, current_seq_p, &my_seqr_chain,
				     m_bufi_p->max_buf2_cnt*fa_workers, m_msp->ldb_info.maxn);

      aa1 = current_seq_p->aa1b;

    loop2: 
      if (lcont) {
	maxt = m_msp->ldb_info.maxt3;
	memcpy(aa1,aa1save,m_msp->ldb_info.l_overlap);
	aa1ptr= &aa1[m_msp->ldb_info.l_overlap];	/* aa1ptr is where the next GETLIB sequence goes */
					/* aa1 is the beginning of the sequence for do_work() */
	loffset += n1 - m_msp->ldb_info.l_overlap;	/* this must be n1, which is the old value, not current_seq_p->n1 */ 
	ocont = lcont;
      }
      else {
	maxt = m_msp->ldb_info.maxn;
	aa1ptr=aa1;
	if (ocont) *n1tot_cur = n1tot_v;
	ocont = 0;
	loffset = 0l;
	n1tot_v = 0;
	n1tot_cur = NULL;
      }
    } /* end while((n1=getlib())) */
  } /*  end cur_lib_p */

  end_seqr_chain(my_seqr_chain);
  free(aa1save-1);
  return seqr_base;
};

/* called before next_sequence_p to ensure that sequences start from the beginning */
void
reset_seqr_chain(struct seqr_chain *seqr_base) {
  struct seqr_chain *cur_seqr;

  seqr_base->cur_seqr_p = seqr_base;
  
  for (cur_seqr = seqr_base; cur_seqr; cur_seqr = cur_seqr->next) {
    cur_seqr->cur_seq_cnt = cur_seqr->max_seq_cnt;
  }
}

struct seq_record *
next_sequence_p(struct mseq_record **cur_mseq_p, struct seqr_chain *seqr_base) {
  struct seqr_chain *cur_seqr_p;
  struct seq_record *next_seq_p;
  int seq_offset;

  cur_seqr_p = seqr_base->cur_seqr_p;

  /* used the last seq_record in the current chain? */
  if (cur_seqr_p->cur_seq_cnt <= 0) {
    seqr_base->cur_seqr_p = cur_seqr_p = cur_seqr_p->next;
    if (cur_seqr_p == NULL) return NULL;
  }

  seq_offset = cur_seqr_p->max_seq_cnt - cur_seqr_p->cur_seq_cnt;

  *cur_mseq_p = cur_seqr_p->mseqr_base + seq_offset;
  next_seq_p = cur_seqr_p->seqr_base + seq_offset;

#ifdef DEBUG
  if (next_seq_p->n1 == 0 || next_seq_p->aa1b == NULL) {
    fprintf(stderr,"*** ERROR: empty seq_record (next_sequence_p)\n");
  }
#endif

  cur_seqr_p->cur_seq_cnt--;
  return next_seq_p;
};
