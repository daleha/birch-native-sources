/* copyright (c) 1996, 1997, 1998, 1999, 2002  William R. Pearson and the
   U. of Virginia */

/*  $Id: comp_lib5e.c 852 2011-10-27 15:14:53Z wrp $ */
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

char *mp_verstr="";

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

int closelib(struct lmf_str *m_fptr, int force);

void *my_srand();
unsigned int my_nrand(int, void *);

void
getlib_buf_work(unsigned char **aa0, unsigned char *aa0s,
		struct buf_head *lib_bhead_p, struct lib_struct *lib_list_p,
		const struct mng_thr *m_bufi_p, struct mngmsg *m_msp, struct pstruct *ppst,
		void *pstat_void, struct db_str *db, struct db_str *ldb,
		struct hist_str *histp,	struct score_count_s *s_info,
#if !defined(COMP_THR) && !defined(PCOMPLIB)
		unsigned char *aa1shuff, void *f_str, void *qf_str,
#endif
		struct seq_record *best_seqs, struct mseq_record *best_mseqs, struct beststr *best,
		FILE *fdata);

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

void
preserve_seq(struct buf2_data_s *, struct seq_record *, struct mseq_record *, struct beststr *);

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

extern void 	/* pre-alignment */
pre_load_best(unsigned char *aa1, int maxn,struct beststr **bbp_arr, int nbest,
	      struct mngmsg *m_msp);

extern char *	/* run link file */
build_link_data(char **, struct mngmsg *, struct beststr **, int);

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

struct lib_struct * 
lib_select(char *lname, char *ltitle, const char *flstr, int ldnaseq); /* lib_sel.c */

void query_parm(struct mngmsg *, struct pstruct *); /* initfa/sw.c */

/* doinit.c */
void markx_to_m_msp(struct mngmsg *m_msp, struct markx_str *markx);
void m_msp_to_markx(struct markx_str *markx, struct mngmsg *m_msp);

/* compacc.c */
void print_header1(FILE *fd, const char *argv_line,
		   const struct mngmsg *m_msp, const struct pstruct *ppst);
void print_header2(FILE *fd, int qlib, char *info_qlabel, unsigned char **aa0,
		   const struct mngmsg *m_msp, const struct pstruct *ppst, char *info_lib_range_p);
void print_header3(FILE *fd, int qlib, struct mngmsg *m_msp, struct pstruct *ppst);

void print_header4(FILE *fd, char *info_qlabel, char *argv_line, char *info_gstring3,
		   char *info_hstring_p[2], struct mngmsg *m_msp, struct pstruct *ppst);
void print_header4a(FILE *fd, struct mngmsg *m_msp);

void print_header5(FILE *fd, int qlib, struct db_str *qtt,
		   struct mngmsg *m_msp, struct pstruct *ppst, int in_mem, long tot_memK);

void prhist(FILE *, const struct mngmsg *, struct pstruct *, struct hist_str hist,
	    int nstats, int sstats, struct db_str, char *, char *, char **, char **, long);
void print_sum(FILE *, struct db_str db);
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
struct beststr **bestp_arr;	/* array of pointers */
int nbest;	/* number of best scores */

struct stat_str *stats; /* array of scores for statistics from real
			   (or shuffled) sequences*/
struct stat_str *qstats;	/* array of scores for shuffled query stats */
struct stat_str *rstats;	/* array of scores from shuffled library */

  /* these variables are global so they can be set both by the main()
     program and save_best() in threaded mode.
  */
int nstats, nqstats, nrstats, pre_nstats, kstats, shuff_tot, sstats;
double zbestcut;		/* cut off for best z-score */
int bestfull;		/* index for selectbest() */
int stats_done=0;	/* flag for z-value processing */
void *rand_state;

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

  char link_title[MAX_LSTR];
  char *link_lib_str;
  char *link_lib_file;

  int j;
  int shuff_mult, n1lib_req;
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
  struct lib_struct *lib_list_p, *cur_lib_p, *link_lib_list_p, *next_lib_p;

  struct db_str link_db={0,0,0}, link_ldb={0,0,0};
  struct hist_str link_hist;
  struct score_count_s link_s_info={0,0,0,0};

  double db_tt;
  int utmp;		/* user input tmp */

  struct pstruct pst;
  void *f_str[6], *qf_str;	/* different f_str[]'s for forward,reverse */
  int have_f_str=0;

  /* these variables track buffers of library sequences */
  struct buf_head *lib_bhead_p, *t_lib_bhead_p;
  struct buf2_data_s *lib_buf2_dp;

  struct mng_thr m_bufi;	/* has max_work_buf, max_buf2_res */
  int ave_seq_len;
  /*   int empty_reader_bufs; */
#ifdef COMP_THR
  /*   int t_reader_buf_readp; */
  struct thr_str *work_info;
#endif
#ifdef MPI_SRC
  int mpi_tid;
#endif
  /* end of library sequence buffers */

  struct mngmsg m_msg;		/* saves most non-param information
				   about library alignment */
  struct markx_str markx_save;	/* saves m_msg values for markx */
  struct markx_str *cur_markx;	/* follow m_msg.markx_list */

  struct hist_str hist2;	/* hist str for zsflag > 2 */
  int zsflag_save;		/* save zsflag > 20 */
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
  int igncnt=0;			/* count for ignoring sequences warning */
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
  /* document command line */
  argv_line[0]='\0';
  for (i=0; i<argc; i++) {
    SAFE_STRNCAT(argv_line," ",sizeof(argv_line));
    if (strchr(argv[i],' ')) {
      SAFE_STRNCAT(argv_line,"\"",sizeof(argv_line));
      SAFE_STRNCAT(argv_line,argv[i],sizeof(argv_line));
      SAFE_STRNCAT(argv_line,"\"",sizeof(argv_line));
    }
    else {
      SAFE_STRNCAT(argv_line,argv[i],sizeof(argv_line));
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

  /* print argv_line, program, version */
  print_header1(stdout, argv_line, &m_msg, &pst);

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
  lib_list_p = lib_select(m_msg.lname, m_msg.ltitle, m_msg.flstr, m_msg.ldb_info.ldnaseq);

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
     This section defines m_bufi.max_buf2_res, the average number of entries
     per buffer, and m_bufi.max_work_buf, the total number of buffers

     Use a 2 Mbyte (DEF_WORKER_BUF) buffer for each worker.  For
     proteins, that means 5,000 sequences of length 400 (average).
     For DNA, that means 2,000 sequences of length 1000.

     To accommodate larger libraries in memory, use more buffers, not
     bigger buffers.

     Once m_bufi.max_buf2_res/max_work_buf and ave_seq_len are set,
     allocate all the communication buffers with alloc_comp_bufs();
  */

  if (m_msg.ldb_info.ldnaseq== SEQT_DNA) {
    ave_seq_len = AVE_NT_LEN;
    m_bufi.max_buf2_res = DEF_WORKER_BUF/AVE_NT_LEN;
  }
  else {
    ave_seq_len = AVE_AA_LEN;
    m_bufi.max_buf2_res = DEF_WORKER_BUF/AVE_AA_LEN;
  }

  /* however - buffer sizes should be a function of the number of
     workers so that all the workers are kept busy.  Assuming a 10,000
     entry library is the smallest we want to schedule, then

  if (max_buf2_res > 10000/fa_max_workers) 
    max_buf2_res = 10000/(2*fa_max_workers);
  */

  m_bufi.max_buf2_res /= m_msg.thr_fact;

  /* finally, max_buf2_res should be mod 6 for tfasta/s/f */
  m_bufi.max_buf2_res -= (m_bufi.max_buf2_res % 6);

  /* max_work_buf is the number of buffers - if the worker buffers are
     small, then make lots more buffers */

#ifdef PCOMPLIB	/* PCOMPLIB -- one buffer per worker */
  m_bufi.max_work_buf = fa_max_workers;
#else	/* !PCOMPLIB */
  m_bufi.max_work_buf = (DEF_WORKER_BUF * 2 * fa_max_workers)/(ave_seq_len * m_bufi.max_buf2_res);
  if (m_bufi.max_work_buf < 2*fa_max_workers) m_bufi.max_work_buf = 2*fa_max_workers;

  m_bufi.max_work_buf -= (m_bufi.max_work_buf%fa_max_workers);
#ifndef COMP_THR
  /* if not threaded, only one (larger) buffer */
  m_bufi.max_buf2_res *= m_bufi.max_work_buf;
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

  m_msg.db.length = 0l;
  m_msg.db.entries = m_msg.db.carry = 0;

  outfd = stdout;  

  /* main loop for doing a search, getting the next query */
  while(1) {

    /* Initialize bestp_arr */
    for (nbest = 0; nbest < MAX_BEST; nbest++)
      bestp_arr[nbest] = &best[nbest];
    nbest = 0;

    qlib++;
    stats_done = 0;

    /* also sets ldb_info.l_overlap */
    /* happens inside while (1) {query} loop for m_msg.n0 */
    m_msg.ldb_info.maxn = maxn = reset_maxn(&m_msg,m_msg.n0, m_msg.max_tot);
    pst.maxlen = maxn;

    zbestcut = -FLT_MAX;
    nstats = nrstats = pre_nstats = shuff_tot = sstats = 0;

    /* ensure that link_list is NULL for no result */
    link_lib_list_p = NULL;	

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

    if (!(m_msg.markx & MX_MBLAST2)) {
	sprintf (&m_msg.qtitle[leng], " - %d %s", m_msg.n0, m_msg.sqnam);
    }

    tprev = s_time();
  
    if (fdata) fprintf(fdata,">>>%ld %3d\t%-50s\n",qtt.entries,m_msg.n0,m_msg.qtitle);

    qtt.length += m_msg.n0;
    qtt.entries++;

#if !defined(COMP_THR) && !defined(PCOMPLIB)
    have_f_str=1;

    /* allocates aa0 arrays, calls init_work, only for non-threaded */
    init_aa0(aa0, m_msg.n0, m_msg.nm0, &aa0s, &aa1shuff,
	     m_msg.qframe, m_msg.qshuffle, m_msg.max_tot,
	     &pst,  &f_str[0], &qf_str, rand_state);
    aa1shuff_b = aa1shuff-1;

    /* label library size limits -- must be called after init_aa0(),
       which calls init_work(), which can reset the n1_high/n1_low
       limits */
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

    getlib_buf_work(aa0, aa0s, lib_bhead_p, lib_list_p, &m_bufi, &m_msg, &pst,
		    &m_msg.pstat_void, &m_msg.db, &m_msg.ldb, &m_msg.hist, &m_msg.s_info,
#if !defined(COMP_THR) && !defined(PCOMPLIB)
		    aa1shuff, f_str, qf_str,
#endif
		    best_seqs, best_mseqs, best, fdata);

#if defined(COMP_THR) && !defined(PCOMPLIB)
    info_lib_range_p = work_info[0].info_lib_range;
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
      pst.zsflag_f = process_hist(rstats,nrstats,&m_msg, &pst, &hist2,
				  &m_msg.pstat_void2,&m_msg.ss_info, 0);
      pst.zs_off = zs_off_save;

      for (i=0; i<nbest; i++) {
	bestp_arr[i]->zscore2 =
	  find_z(bestp_arr[i]->rst.score[pst.score_ix],
		 bestp_arr[i]->rst.escore, bestp_arr[i]->seq->n1,
		 bestp_arr[i]->rst.comp,m_msg.pstat_void2);
      }
    }

    get_param(&pst, info_gstring2p,info_gstring3, &m_msg.s_info);

    /* **************************************************************** */
    /* label Library: output                                            */
    /* **************************************************************** */

    tscan = s_time();

    /* header2 print Query:, Annotation:, Library: */
    print_header2(stdout, qlib, info_qlabel, aa0, &m_msg, &pst, info_lib_range_p);

    if (m_msg.std_output) {
      prhist (stdout, &m_msg, &pst, m_msg.hist, nstats, sstats, m_msg.ldb,
	      (pst.zsflag > 20? hist2.stat_info:NULL),info_lib_range_p,
	      info_gstring2p, info_hstring_p, tscan-tprev);
    }
#ifdef COMP_MLIB
    ttscan += tscan-tprev;
#endif

    /* fprintf(stderr, " *** shuff_tot: %d ***\n",shuff_tot); */

  l3:
    if (!m_msg.quiet) {
      printf("Enter filename for results [%s]: ", m_msg.outfile);
      fflush(stdout);
      rline[0]='\0';
      if (fgets(rline,sizeof(rline),stdin)==NULL) goto end_l;
      else {	/* parse rline input */
	if ((bp=strchr(rline,'\n'))!=NULL) *bp = '\0';
	if (rline[0]!='\0') strncpy(m_msg.outfile,rline,sizeof(m_msg.outfile));
      }
    }

    if (m_msg.outfile[0]!='\0') {	/* have an output file name */
      if ((outfd=fopen(m_msg.outfile,"w"))==NULL) {
	fprintf(stderr," could not open %s\n",m_msg.outfile);
	if (!m_msg.quiet) goto l3;
	else goto l4;		/* skip output file */
      }

      if (m_msg.markx_list==NULL) {	/* no -m 9 options, need one */
	if ((m_msg.markx_list = (struct markx_str *)calloc(1,sizeof(struct markx_str)))==NULL) {
	  fprintf(stderr,"cannot allocate m_msg.markx_list\n");
	  goto l4;
	}
	else {
	  m_msp_to_markx(m_msg.markx_list, &m_msg);
	}
      }
      m_msg.markx_list->out_file = m_msg.outfile;
      m_msg.markx_list->out_fd = m_msg.outfd = outfd;
      /* save permanent m_msg marks info */
    }
      
  l4:   
    m_msp_to_markx(&markx_save, &m_msg);
    /* for each entry in markx_list */
    for (cur_markx = m_msg.markx_list; cur_markx;
	 cur_markx = cur_markx->next) {

      if (cur_markx->out_file && cur_markx->out_file[0] && cur_markx->out_fd == NULL) {
	if ((cur_markx->out_fd=fopen(cur_markx->out_file,"w"))==NULL) {
	  fprintf(stderr," could not open %s\n",cur_markx->out_file);
	}
      }
      if (cur_markx->out_fd == NULL) continue;

      markx_to_m_msp(&m_msg, cur_markx);

      if (qlib==1) print_header1(cur_markx->out_fd, argv_line, &m_msg, &pst);
      print_header2(cur_markx->out_fd, qlib, info_qlabel, aa0, &m_msg, &pst, info_lib_range_p);

      if (m_msg.std_output) {
	prhist(cur_markx->out_fd, &m_msg, &pst,m_msg.hist, nstats, sstats, m_msg.db,
	       (pst.zsflag > 20? hist2.stat_info:NULL), info_lib_range_p,
	       info_gstring2p, info_hstring_p, tscan-tprev);
      }

      /* done with this output, restore m_msg */
    }
    markx_to_m_msp(&m_msg,&markx_save);

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
    }
    else {
      m_msg.nskip = 0;
    }

    if (m_msg.quiet || m_msg.markx & MX_M9SUMM) {

      /* to determine how many sequences to re-align (either for
	 do_opt() or calc_id() we need to modify m_msg.mshow to get
	 the correct number of alignments */

      if (pst.zsflag >= 0) {	/* do we have e_values? */
	for (i=m_msg.nskip; i<nbest && bestp_arr[i]->rst.escore < m_msg.e_cut; i++) {}

	if (m_msg.mshow_set != 1) {
	  m_msg.nshow = min(i - m_msg.nskip, nbest-m_msg.nskip);
	}
	else {	/* mshow_set>0 */
	  if (m_msg.mshow_min) {	/* must show at least m_msg.mshow results */
	    m_msg.nshow = max(m_msg.mshow, i-m_msg.nskip);
	  }
	  else {	/* limit by e_cut */
	    m_msg.nshow = min(m_msg.mshow, i-m_msg.nskip);
	  }
	}
      }

      if (m_msg.nshow <= 0) { /* no results to display */
	if (m_msg.std_output) fprintf(outfd,"!! No sequences with E() < %0.5g\n",m_msg.e_cut);
	m_msg.nshow = 0;
	goto end_l;
      }
    }


    /* (a) check to see if there are additional sequences in the
       link-file

       (b) if so, read those sequences, score them, and add them to
       best_str if appropriate
     */

    m_msg.pre_load_done = 0;
    if (m_msg.link_lname[0]) {
      /* guarantee that we have the bline's for the best sequences */
      pre_load_best(aa1save, maxn, &bestp_arr[m_msg.nskip], m_msg.nshow, &m_msg);

      if ((link_lib_str = build_link_data(&link_lib_file, &m_msg, bestp_arr,pst.debug_lib))==NULL) {
	goto no_links;
      };

      /* get a list of files */
      link_lib_list_p = lib_select(link_lib_str, link_title, m_msg.flstr, m_msg.ldb_info.ldnaseq);
      if (link_lib_str != NULL) free(link_lib_str);

      /* read in the sequences */
      /* calculate scores for the sequences */

#if defined(COMP_THR) || defined(PCOMPLIB)
      get_rbuf(&lib_bhead_p,m_bufi.max_work_buf);
#else
      lib_bhead_p = lib_buf2_list;
#endif

      getlib_buf_work(aa0, aa0s, lib_bhead_p, link_lib_list_p, &m_bufi, &m_msg, &pst, &m_msg.pstat_void, 
		      &link_db, &link_ldb, NULL, &link_s_info, 
#if !defined(COMP_THR) && !defined(PCOMPLIB)
		      aa1shuff, f_str, qf_str,
#endif
		      best_seqs, best_mseqs, best, NULL);

      m_msg.pre_load_done = 0;	/* pre-load has not been done for link library sequences */

      /* need to resort results, and re-check how many should be
	 displayed */
      scale_scores(bestp_arr,nbest,m_msg.db, &pst,m_msg.pstat_void);

      if (pst.zsflag >= 0) {
	/* skip entries if -F e_low specified */
	for (i=0; i<nbest && bestp_arr[i]->rst.escore < m_msg.e_low; i++) {};
	m_msg.nskip = i;
      }
      else {
	/* no statistics, just use the same score */
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

	  if (m_msg.mshow_set != 1) {
	    m_msg.nshow = min(i - m_msg.nskip, nbest-m_msg.nskip);
	  }
	  else {	/* mshow_set>0 */
	    if (m_msg.mshow_min) {	/* must show at least m_msg.mshow results */
	      m_msg.nshow = max(m_msg.mshow, i-m_msg.nskip);
	    }
	    else {	/* limit by e_cut */
	      m_msg.nshow = min(m_msg.mshow, i-m_msg.nskip);
	    }
	  }
	}
      }
    }
    /* done with -e link_lname */

  no_links:
    if (m_msg.quiet && ((m_msg.stages > 1) || (m_msg.markx & MX_M9SUMM))) {
      /* pre-load sequence data for  alignments for showbest, showalign */
      pre_load_best(aa1save, maxn, &bestp_arr[m_msg.nskip], m_msg.nshow, &m_msg);

      /* must calculate repeat_thresh before buf_align_seq */
      if (pst.do_rep) {
	for (i=m_msg.nskip; i < m_msg.nskip+m_msg.nshow; i++) {
	  bestp_arr[i]->repeat_thresh = 
	    min(E1_to_s(pst.e_cut_r, m_msg.n0, bestp_arr[i]->seq->n1,
			pst.zdb_size, m_msg.pstat_void),bestp_arr[i]->rst.score[pst.score_ix]);
	}
      }

      buf_align_seq(aa0, m_msg.n0, &bestp_arr[m_msg.nskip], m_msg.nshow,
		    &pst, &m_msg, &m_bufi
#if !defined(COMP_THR) && !defined(PCOMPLIB)
		    , f_str
#endif
		    );
    }

    print_header3(stdout, qlib, &m_msg, &pst);

    if (m_msg.do_showbest) {
      showbest(stdout, aa0, aa1save, maxn, &bestp_arr[m_msg.nskip], nbest-m_msg.nskip,
	       qtt.entries, &m_msg, &pst,m_msg.db, info_gstring2p, f_str);

      m_msp_to_markx(&markx_save, &m_msg);
      t_quiet = m_msg.quiet;
      m_msg.quiet = -1;	/* should guarantee 1..nbest shown */
      for (cur_markx = m_msg.markx_list; cur_markx; cur_markx = cur_markx->next) {
	if (cur_markx->out_fd == NULL) continue;
	markx_to_m_msp(&m_msg, cur_markx);
	print_header3(cur_markx->out_fd, qlib, &m_msg, &pst);
	showbest(cur_markx->out_fd, aa0, aa1save, maxn, &bestp_arr[m_msg.nskip], nbest-m_msg.nskip,
		 qtt.entries, &m_msg, &pst, m_msg.db, info_gstring2p, f_str);
      }
      m_msg.quiet = t_quiet;
      markx_to_m_msp(&m_msg, &markx_save);
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
      }

      /* alignment transition
	 >>>query vs library for -m 9, -m 10
	 "\n" for MX_MBLAST
	 ; pg_name, other info for -m 10
	 
	 not sent to stdout if outfile specified
       */
      print_header4(outfd, info_qlabel, argv_line, info_gstring3, info_hstring_p, &m_msg, &pst);
      showalign (outfd, aa0, aa1save, maxn,
		 &bestp_arr[m_msg.nskip], nbest-m_msg.nskip,
		 qtt.entries, &m_msg, &pst, info_gstring2p, f_str, &m_bufi);

      fflush(outfd);
    }

    m_msp_to_markx(&markx_save, &m_msg);
    for (cur_markx = m_msg.markx_list; cur_markx; cur_markx=cur_markx->next) {
      if (cur_markx->out_fd == NULL) continue;
      if (cur_markx->out_fd == outfd) continue;
      markx_to_m_msp(&m_msg, cur_markx);
      print_header4(cur_markx->out_fd, info_qlabel, argv_line, info_gstring3, info_hstring_p, &m_msg, &pst);
      showalign (cur_markx->out_fd, aa0, aa1save, maxn,
		 &bestp_arr[m_msg.nskip], nbest-m_msg.nskip,
		 qtt.entries, &m_msg, &pst, info_gstring2p, f_str, &m_bufi);
      fflush(cur_markx->out_fd);
    }
    markx_to_m_msp(&m_msg, &markx_save);

  end_l:

    print_header4a(outfd, &m_msg);
    for (cur_markx = m_msg.markx_list; cur_markx; cur_markx=cur_markx->next) {
      if (cur_markx->out_fd == NULL) continue;
      if (cur_markx->out_fd == outfd) continue;
      print_header4a(cur_markx->out_fd, &m_msg);
    }

    if (fdata) {
      fprintf(fdata,"/** Algorithm : %s  **/\n",info_gstring2p[0]);
      fprintf(fdata,"/** Parameters : %s  **/\n",info_gstring2p[1]);
      fprintf(fdata,"%3ld%-50s\n",qtt.entries-1,m_msg.qtitle);
      fflush(fdata);
    }
    
    /* **************************************************************** */
    /* all done, start resetting everything                             */
    /* **************************************************************** */

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
      memset(lib_bhead_p->buf2_data,0,(size_t)(m_bufi.max_buf2_res+1)*sizeof(struct buf2_data_s));
      /* replace it */
      lib_bhead_p->buf2_data[0].seq = lib_bhead_p->hdr.seq_b;
      lib_bhead_p->buf2_data[0].mseq = lib_bhead_p->hdr.mseq_b;
      lib_bhead_p->buf2_data[0].seq->aa1b = lib_bhead_p->hdr.aa1b_start;
      lib_bhead_p->hdr.seq_record_continuous = 1;  /* should be 1 */
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
	if (closelib(cur_lib_p->m_file_p,0)) {
	  cur_lib_p->m_file_p = NULL;
	}
      }
    }

    if (link_lib_list_p) {
      for (cur_lib_p = link_lib_list_p; cur_lib_p != NULL; cur_lib_p = next_lib_p) {
	next_lib_p = cur_lib_p->next;
	if (cur_lib_p->m_file_p !=NULL) {
	  closelib(cur_lib_p->m_file_p,1);
	}
	free(cur_lib_p->file_name);
	free(cur_lib_p->m_file_p);
	free(cur_lib_p);
      }
      /* delete the library file */
      if (!pst.debug_lib) {
	unlink(link_lib_file);
      }
      if (link_lib_file) free(link_lib_file);
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

  next_query:
    m_msg.q_offset = next_q_offset;

    m_msg.n0 = 
      QGETLIB (aa0[0], MAXTST, m_msg.qtitle, sizeof(m_msg.qtitle),
	       &qseek, &qlcont,q_file_p,&m_msg.q_off);
    if (m_msg.n0 < 0) break;
    if (m_msg.n0 == 0) { 
      next_q_offset = 0;
      goto next_query;
    }
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

    if (qlcont) {
      next_q_offset = m_msg.q_offset + m_msg.n0 - m_msg.q_overlap;
    }
    else {
      next_q_offset = 0l;
    }
  }	/* end of while(1) for multiple queries */

#ifdef PCOMPLIB
  /* tell workers to quit */
  init_thr(fa_max_workers, NULL, NULL, NULL, NULL, NULL);
#endif

  tdone = s_time();

  /* 
     218 residues in 1 query   sequences
     5190103 residues in 13351 library sequences
     Tcomplib [36.3.5 Apr, 2011(preload7)] (2 proc)
     ...

     goes to stdout regardless of outfd
  */
  print_header5(stdout, qlib, &qtt, &m_msg, &pst, 0, 0);

  m_msp_to_markx(&markx_save, &m_msg);
  for (cur_markx = m_msg.markx_list; cur_markx; cur_markx=cur_markx->next) {
    if (cur_markx->out_fd == NULL) continue;
    markx_to_m_msp(&m_msg, cur_markx);
    print_header5(cur_markx->out_fd, qlib, &qtt, &m_msg, &pst, 0, 0);
    fflush(cur_markx->out_fd);
    fclose(cur_markx->out_fd);
  }
  markx_to_m_msp(&m_msg, &markx_save);

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
  fprintf(fd," Scomplib [%s%s]\n start: %s done: %s\n",verstr,mp_verstr,tstr1,tstr2);
#endif
#if defined(COMP_THR)
  fprintf(fd," Tcomplib [%s%s] (%d proc)\n start: %s done: %s\n", verstr,mp_verstr,
    fa_max_workers,tstr1,tstr2);
#endif
#if defined(PCOMPLIB)
  fprintf(fd," Pcomplib [%s%s] (%d proc)\n start: %s done: %s\n", verstr,mp_verstr,
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
  fprintf (fd, "\nFunction used was %s [%s%s]\n", prog_func,verstr,mp_verstr);
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

void preserve_seq(struct buf2_data_s *lib_buf2_dp,
		  struct seq_record *best_seqs,
		  struct mseq_record *best_mseqs,
		  struct beststr *best) {
  struct seq_record *dest_seq_p, *saved_seq_p;
  struct mseq_record *dest_mseq_p, *saved_mseq_p;
  struct beststr *next_bbp;

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

/* ****************************************************************
getlib_buf_work() is the comp_lib5e.c version of seqr_chain_work().
seqr_chain_work() works with a linked list of sequence buffers in memory
getlib_buf_work() (re-)reads the database for each query.

getlib_buf_work() distinguishes m_msg.db -- determined from reading
the database -- from m_msg.ldb -- set in save_best().  So both are
required here.  In addition, getlib_buf_work() uses *best_seqs/
*best_mseqs and best[] to save high scoring sequences so they need not
be re-read.

getlib_buf_work() and seqr_chain_work() allow additional sequences to
be compared and added to beststr during the expansion process, and
ensure that exactly the same process is followed for both libraries.
**************************************************************** */
void
getlib_buf_work(unsigned char **aa0, unsigned char *aa0s,
		struct buf_head *lib_bhead_p, struct lib_struct *lib_list_p,
		const struct mng_thr *m_bufi_p, struct mngmsg *m_msp, struct pstruct *ppst,
		void *pstat_void, struct db_str *db, struct db_str *ldb, struct hist_str *histp,
		struct score_count_s *s_info,
#if !defined(COMP_THR) && !defined(PCOMPLIB)
		unsigned char *aa1shuff, void *f_str, void *qf_str,
#endif
		struct seq_record *best_seqs, struct mseq_record *best_mseqs, struct beststr *best, FILE *fdata) {

  unsigned char *aa1ptr, *aa1, *aa1save;
  char *bp;
  int n1;
  int seq_index;
  struct lib_struct *cur_lib_p;
  struct lmf_str *m_file_p;
  struct buf2_data_s *lib_buf2_dp;
  struct seq_record *current_seq_p;
  struct mseq_record *current_mseq_p;
  long cur_aa1b_buf_size;			/* current length (in residues/bytes) of buffer */
  int i, itt, jstats;
  int seqr_chain_cnt, seqr_cnt;	/* count through sequences in seqr_chain */
  int buf2_shuff_mask;		/* pre-calculate when shuffles requested */
  int lcont, ocont, maxt;	/* continued sequence */

  int *n1tot_ptr=NULL, *n1tot_cur;
  int n1tot_cnt=0;
  int n1tot_v, aa1_loff;
  long loffset, l_off;
  fseek_t lseek;		/* seek into library of current sequence */

  char libstr[MAX_UID];		/* required because getlib() does not return lcont > 0 */
  int n_libstr=MAX_UID;			/* length of libstr */

  /* set shuffle mask for -z > 10 */
  buf2_shuff_mask = 0;
  if (ppst->zsflag >= 10 && ppst->zsflag < 20) {
    buf2_shuff_mask = BUF2_DOSHUF;
  }
  seq_index = 0;

  if ((aa1save = (unsigned char *)calloc((m_msp->ldb_info.maxn+1),sizeof (char))) == NULL) {
    s_abort ("Unable to allocate library overlap", "");
  }
  *aa1save= '\0';
  aa1save++;

  lib_bhead_p->hdr.buf2_cnt = 0;
  lib_bhead_p->hdr.have_results = 0;
  lib_bhead_p->hdr.aa1b_used = 0;
  lib_bhead_p->hdr.stop_work = 0;

  lib_bhead_p->hdr.buf2_type = (BUF2_DOWORK | buf2_shuff_mask);

  lib_buf2_dp = lib_bhead_p->buf2_data;
  cur_aa1b_buf_size = 0;

  /* open the library - start the search */

  /* step through each of the library files on the list */
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

    /* read sequence directly into buffer */

    current_seq_p = lib_buf2_dp->seq;
    current_mseq_p = lib_buf2_dp->mseq;
    aa1ptr = aa1 = current_seq_p->aa1b;

    while ((n1=GETLIB(aa1ptr, maxt, libstr, n_libstr,
		      &lseek,	&lcont, m_file_p,
		      &(current_seq_p->l_off)))>=0) {

#ifdef DEBUG
      /* check for out of range sequence */
      for (i=0; i<n1; i++) {
	if (aa1[i] > ppst->nsq_e) {
	  fprintf(stderr,"[%s/comp_lib5.c] %s[%d] = %d > %d out of range\n",
		  prog_func,libstr, i, aa1[i], ppst->nsq_e);
	  aa1[i] = 1;
	}
      }
#endif

      current_seq_p->index = seq_index;
      current_mseq_p->index = seq_index++;
      current_seq_p->n1 = n1;
      current_mseq_p->m_file_p = (void *)m_file_p;
      current_mseq_p->cont = ocont+1;
      current_seq_p->l_offset = loffset;
      current_mseq_p->lseek = lseek;

      if ((bp=strchr(libstr,' '))!=NULL) *bp='\0';
      strncpy(current_mseq_p->libstr,libstr,MAX_UID);	/* get old libstr for lcont>0 */

      /* add termination code for FASTX/FASTY if necessary */
      if (m_msp->ldb_info.term_code && !lcont &&
	  m_msp->ldb_info.ldnaseq==SEQT_PROT && aa1ptr[n1-1]!=m_msp->ldb_info.term_code) {
	aa1ptr[n1++]=m_msp->ldb_info.term_code;
	aa1ptr[n1]=0;
      }

#ifdef DEBUG
      if (n_libstr <= MAX_UID) {if ((bp=strchr(current_mseq_p->libstr,' '))!=NULL) *bp='\0';}
      if (aa1[-1]!='\0' || aa1ptr[n1]!='\0') {
	fprintf(stderr,"%s: aa1[%d] at %ld:%lld  missing NULL boundaries: %d %d\n",
		current_mseq_p->libstr,n1, ldb->entries+1,current_mseq_p->lseek,
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

      db->entries++;
      db->length += n1;
      if (db->length > LONG_MAX) {
	db->length -= LONG_MAX; db->carry++;
      }

      /* don't count long sequences more than once */
      if (aa1!=aa1ptr) {	/* this is a continuation */
	current_seq_p->n1 = n1 += m_msp->ldb_info.l_overlap;	/* corrected 28-June-2008 */
	ldb->entries--;
      }

#ifdef DEBUG
      /* This finds most reasons for core dumps */
      if (ppst->debug_lib)
	for (i=0; i<n1; i++) {
	  if (aa1[i]>ppst->nsq || aa1[i] <= 0) {
	    fprintf(stderr, 
		    "[%s/comp_lib5.c] %s residue[%d/%d] %d range (%d) lcont/ocont: %d/%d\n%s\n",
		    prog_func,current_mseq_p->libstr,i,current_seq_p->n1,aa1[i],ppst->nsq,lcont,ocont,aa1ptr+i);
	    aa1[i]=0;
	    n1=i-1;
	    break;
	  }
	}

      current_seq_p->adler32_crc =
	current_mseq_p->adler32_crc = adler32(1L,aa1,current_seq_p->n1);
#endif

      if ( n1 <= 1) {
	/*	if (igncnt++ <10)
		fprintf(stderr,"Ignoring: %s\n",current_seq_p->libstr);
	*/
	goto loop2;
      }

      /* check to see whether this score (or a shuff score) should
	 be included in statistics */
      jstats = samp_stats_idx(&pre_nstats, nstats, rand_state);

      cur_aa1b_buf_size += n1+1;
      lib_bhead_p->hdr.aa1b_used += n1+1;

#ifdef PCOMPLIB
      lib_buf2_dp->seq_dup = 0;	/* original, not duplicate */
#endif
      for (itt=m_msp->revcomp; itt<=m_msp->nitt1; itt++) {
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

      if (lcont) {
	memcpy(aa1save,&aa1[n1-m_msp->ldb_info.l_overlap],m_msp->ldb_info.l_overlap);
      }

      /* if the buffer is filled */
      if (lib_bhead_p->hdr.buf2_cnt >= m_bufi_p->max_buf2_res || cur_aa1b_buf_size >= m_bufi_p->seq_buf_size - m_msp->ldb_info.maxn) {

	/* fprintf(stderr," new empty buffer at: %lld\n", current_seq_p->lseek); */

#if defined(COMP_THR) || defined(PCOMPLIB)		/* if COMP_THR - fill and empty buffers */
	/* provide filled buffer to workers */
	/*
	  fprintf(stderr,"to:[%d] aa1b_used[1]: %d; aa1b_used[2] %d\n", lib_bhead_p->hdr.worker_idx,
	  lib_bhead_p->hdr.aa1b_used,  aa1 - lib_bhead_p->hdr.aa1b_start + n1+1);
	*/
	lib_bhead_p->hdr.have_data = 1;
	put_rbuf(lib_bhead_p,m_bufi_p->max_work_buf);
	/* empty_reader_bufs--; */
	/* get an empty buffer to fill */
	get_rbuf(&lib_bhead_p,m_bufi_p->max_work_buf);
	/* empty_reader_bufs++; */
#else			/* just do the searches */
	if (lib_bhead_p->hdr.buf2_type & BUF2_DOWORK) {
	  buf_do_work(aa0, m_msp->n0, lib_bhead_p, ppst, f_str);
	  if (m_msp->qshuffle)
	    buf_qshuf_work(aa0s,m_msp->n0, lib_bhead_p, ppst, qf_str, ppst->score_ix);
	}
	if (lib_bhead_p->hdr.buf2_type & BUF2_DOSHUF) {
	  buf_shuf_work(aa0,m_msp->n0, aa1shuff, lib_bhead_p, ppst,
			f_str, ppst->score_ix,rand_state);
	}
#endif

	/* "empty" buffers have results that must be processed */
	if (lib_bhead_p->hdr.buf2_cnt && lib_bhead_p->hdr.have_results) {
	  save_best(lib_bhead_p,m_msp, ppst, ldb, fdata,
		    histp, pstat_void, s_info);

	  /* this section of code is only used for re-cycled buffers */
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
	}

	/* now the buffer is truly empty, fill it up */
	lib_bhead_p->hdr.buf2_cnt = 0;
	lib_bhead_p->hdr.aa1b_used = 0;
	lib_bhead_p->hdr.have_results = 0;
	lib_bhead_p->hdr.stop_work = 0;
	lib_bhead_p->hdr.buf2_type = (BUF2_DOWORK | buf2_shuff_mask);
	lib_buf2_dp = lib_bhead_p->buf2_data;
	current_seq_p = lib_buf2_dp->seq = lib_bhead_p->hdr.seq_b;
	current_mseq_p = lib_buf2_dp->mseq = lib_bhead_p->hdr.mseq_b;
	aa1 = current_seq_p->aa1b = lib_bhead_p->hdr.aa1b_start;
	lib_bhead_p->hdr.aa1b_used = 0;
	cur_aa1b_buf_size = 0;
      }
      else {	/* room left in current buffer, increment ptrs */
	current_seq_p++;
	current_mseq_p++;
	lib_buf2_dp->seq = current_seq_p;
	lib_buf2_dp->mseq = current_mseq_p;
	aa1 = current_seq_p->aa1b = current_seq_p[-1].aa1b+n1+1;
      }

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
    /* all done reading data */

  free(--aa1save);

#if defined(COMP_THR) || defined(PCOMPLIB)
    /*
    fprintf(stderr,"to:[%d] aa1b_used[1]: %d; aa1b_used[2] %d\n", lib_bhead_p->hdr.worker_idx,
	    lib_bhead_p->hdr.aa1b_used,  aa1 - lib_bhead_p->hdr.aa1b_start + n1+1);
    */
    lib_bhead_p->hdr.have_data = 1;	/* ignored if buf2_cnt <= 0 */
    put_rbuf(lib_bhead_p,m_bufi_p->max_work_buf);

    /* wait for the threads to finish */
    wait_rbuf(m_bufi_p->max_work_buf);
    /* wait_rbuf(m_bufi_p->max_work_buf - empty_reader_bufs); */

    /* save the final results */
    for (i=0; i < num_reader_bufs; i++) {
      save_best(RESULTS_BUF[i],m_msp, ppst, ldb, fdata,
		histp, pstat_void, s_info);
      RESULTS_BUF[i]->hdr.buf2_cnt = RESULTS_BUF[i]->hdr.have_results = 0;
    }
#else
    if (lib_bhead_p->hdr.buf2_type & BUF2_DOWORK) {
      buf_do_work(aa0, m_msp->n0, lib_bhead_p, ppst, f_str);
      if (m_msp->qshuffle)
	buf_qshuf_work(aa0s,m_msp->n0, lib_bhead_p, ppst, qf_str, ppst->score_ix);
    }

    if (lib_bhead_p->hdr.buf2_type & BUF2_DOSHUF) {
      buf_shuf_work(aa0, m_msp->n0, aa1shuff, lib_bhead_p,
		    ppst, f_str, ppst->score_ix, rand_state);
    }

    save_best(lib_bhead_p,m_msp, ppst, ldb, fdata,
	      histp, pstat_void, s_info);
    lib_bhead_p->hdr.buf2_cnt = lib_bhead_p->hdr.have_results = 0;
#endif
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

  buf_siz = max(m_bufi_p->max_buf2_res*ave_seq_len, m_msp->max_tot * 4);
  if (buf_siz < m_msp->max_tot) buf_siz = m_msp->max_tot;
  m_bufi_p->seq_buf_size = buf_siz;

  for (i=0; i<m_bufi_p->max_work_buf; i++) {
    /* allocate max_buf2_res buf2_str's into each buf2 */

    if ((lib_buf2_list[i].buf2_ares =
	 (struct buf2_ares_s *)calloc((size_t)(m_bufi_p->max_buf2_res+1),
				   sizeof(struct buf2_ares_s)))
        ==NULL) {
      fprintf(stderr," cannot allocate buffer struct %d %d\n",
	      i,m_bufi_p->max_buf2_res+1);
      exit(1);
    }

    if ((lib_buf2_list[i].buf2_res =
	 (struct buf2_res_s *)calloc((size_t)(m_bufi_p->max_buf2_res+1),
				   sizeof(struct buf2_res_s)))
        ==NULL) {
      fprintf(stderr," cannot allocate buffer struct %d %d\n",
	      i,m_bufi_p->max_buf2_res+1);
      exit(1);
    }

    if ((lib_buf2_list[i].buf2_data =
	 (struct buf2_data_s *)calloc((size_t)(m_bufi_p->max_buf2_res+1),
				   sizeof(struct buf2_data_s)))
        ==NULL) {
      fprintf(stderr," cannot allocate buffer struct %d %d\n",
	      i,m_bufi_p->max_buf2_res+1);
      exit(1);
    }

    if ((lib_buf2_list[i].buf2_data[0].mseq = lib_buf2_list[i].hdr.mseq_b = 
	 (struct mseq_record *)calloc((size_t)(m_bufi_p->max_buf2_res+1),
				     sizeof(struct mseq_record)))
        ==NULL) {
      fprintf(stderr," cannot allocate mseq_record buffer %d %d\n",
	      i,m_bufi_p->max_buf2_res+1);
      exit(1);
    }

    if ((lib_buf2_list[i].buf2_data[0].seq = lib_buf2_list[i].hdr.seq_b = 
	 (struct seq_record *)calloc((size_t)(m_bufi_p->max_buf2_res+1),
				     sizeof(struct seq_record)))
        ==NULL) {
      fprintf(stderr," cannot allocate seq_record buffer %d %d\n",
	      i,m_bufi_p->max_buf2_res+1);
      exit(1);
    }
    lib_buf2_list[i].hdr.seq_record_continuous = 1;	/* should be 1 */

    if ((lib_buf2_list[i].hdr.aa1b_start =
         (unsigned char *)calloc((size_t)(buf_siz+1),sizeof(unsigned char)))
        ==NULL) {
      fprintf(stderr," cannot allocate buffer %d: %d\n",i, buf_siz);
      exit(1);
    }

    /* make certain there is a '\0' at the beginning */
    lib_buf2_list[i].hdr.aa1b_start++;
    lib_buf2_list[i].buf2_data[0].seq->aa1b = lib_buf2_list[i].hdr.aa1b_start;
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
