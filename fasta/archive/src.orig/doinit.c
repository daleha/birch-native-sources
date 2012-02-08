/*	doinit.c	general and function-specific initializations */

/* copyright (c) 1996, 1997, 1998  William R. Pearson and the U. of Virginia */


/*  $Id: doinit.c 768 2011-05-25 19:58:29Z wrp $ */
/* $Revision: 768 $  */

/* this file performs general initializations of search parameters

   In addition, it calls several functions in init??.c that provide
   program-specific initializations:

   f_initenv()	- called from initenv()
   f_getopt()	- called from initenv() during a getopt() scan
   f_getarg()	- called from initenv() after the getopt() scan

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(UNIX) || defined(_MACH)
#include <unistd.h>
#endif
#ifndef PCOMPLIB
#ifdef IRIX
#include <sys/sysmp.h>
#endif
#else
#include "msg.h"	/* need for FIRSTNODE */
#ifdef MPI_SRC
#include "mpi.h"
#endif
#endif

#include "defs.h"
#include "param.h"
#include "upam.h"	/* required for 'U' option change of nascii */

#include "structs.h"

#define XTERNAL
#include "uascii.h"
#undef XTERNAL

#ifdef UNIX
#include <getopt.h>
#else
extern int optind;		/* used by getopt() */
extern char *optarg;
#endif

char prog_name[MAX_FN];

extern void f_initenv(struct mngmsg *, struct pstruct *, unsigned char **);
extern void f_lastenv(struct mngmsg *, struct pstruct *);
extern void f_getopt(char, char *, struct mngmsg *, struct pstruct *);
extern void f_getarg(int, char **, int, struct mngmsg *, struct pstruct *);
extern void show_help(char *, int pgm_id);
extern void show_all_help(char *pgm_name, int pgm_id);
void g_init_opts(struct mngmsg *, struct pstruct *);

void add_ascii_ann(int *qascii, char *ann_arr);
static int set_markx(int markx, int val);
static void pre_parse_markx(char *opt_arg, struct mngmsg *m_msp);
static void parse_markx(char *opt_arg, struct markx_str *this_markx);
void markx_to_m_msp(struct mngmsg *m_msp, struct markx_str *this_markx);
void m_msp_to_markx(struct markx_str *this_markx, struct mngmsg *m_msp);

int optcnt;
int fa_max_workers=MAX_WORKERS;
#ifdef PCOMPLIB
int worker_1=0;
int worker_n=0;
#endif

extern struct opt_def_str f_options[];

void set_opt_disp_defs(char opt_char, struct opt_def_str *options, int type,
		       int i_param1, int i_param2,
		       double d_param1, double d_param2, char *s_param);

/* ****************************************************************
   The option/-help system has been substantially restructured to
   allow more consistent -h/-help messages.

   There are now two global arrays, opt_def_str g_options (global
   options, parsed in doinit.c), and opt_def_str f_options
   (function-specific options, parsed in initfa.c)

   struct opt_def_str {
     char opt_char;	# getopt single character option letter
     int has_arg;	# does it have an option?
     char *opt_str;	# getopt_long (future) long option name
     char *opt_descr_s;	# short description of option
     char *opt_descr_l; # long description of option (if NULL, use opt_descr_s)
     int opt_rank;	# rank of option (not used)
     int fmt_type;	# fmt type (for defaults): 1,2 ints, 3,4 doubles
     int i_param1;	# int default1
     int i_param2;	
     double d_param1;	# double default1
     double d_param2;
   };

   the g_opt_string and f_opt_string's parsed by getopt() are built
   from these structures, guaranteeing that the options and help
   messages are kept in sync.

   long options descriptions (opt_descr_l) are saved in static arrays
   (e.g. m_opt_descr[] in doinit.c, z_opt_descr[], s_opt_descr[] in
   initfa.c

   The default option values, which are displayed from i_param[1,2],
   d_param[1,2], are set by g_init_opts() and f_init_opts() using
   set_opt_disp_defs().  g_init_opts()/f_init_opts() should be called
   as late as possible in the program.

   **************************************************************** */

static char m_opt_descr[] ="Output/alignment format;\n      0 - standard \":. \" alignment; 1 - \" xX\"; 2 - \".MS..\"; 3 - separate >fasta entries;\n      4 - \"---\" alignment map; 5 - 0+4; 6 - <html>;\n      8 - BLAST tabular; 8C commented BLAST tabular;\n      B - BLAST Query/Sbjct alignments; BB - complete BLAST output;\n      9 - FASTA tabular; 9c - FASTA tabular encoded; 9C FASTA tabular CIGAR encoded;\n     10 - parseable key:value; 11 - lav for LALIGN;\n      F - 'F0,6,9c out_file' - alternate output formats to files;";

struct opt_def_str g_options[] = {
  {'C', 1, "aname_length", "length of the query/sbjct name in alignments", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'D', 0, "debug", "enable debugging output", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'e', 1, "expand", "expand_script to extend hits", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'F', 1, "evalue_min", "min E()-value displayed", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
#if defined(PCOMPLIB) || !defined(SHOW_HIST)
  {'H', 0, "histogram", "show histogram", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
#else
  {'H', 0, "nohist", "no histogram", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
#endif
  {'i', 0, "revcomp", "search with reverse-complement", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
#ifdef SHOW_HELP
  {'I', 0, "interact", "interactive mode", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
#endif
  {'l', 1, "fastlibs", "FASTLIBS abbreviation file", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'L', 0, "long_info", "long library descriptions", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'m', 1, "outfmt", "output format", &m_opt_descr[0], 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'N', 1, "lib_length", "max library length before overlapping", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'o', 1, "offsets", "offset coordinates of query/subject", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'O', 1, "out", "write results to file", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'q', 0, "quiet", "quiet -- do not prompt", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
#ifndef SHOW_HELP
  {'Q', 0, "\0", "quiet -- do not prompt", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
#else
  {'q', 0, "quiet", "[default] quiet -- do not prompt", NULL, 0},
  {'Q', 0, "\0", "[default] quiet -- do not prompt", NULL, 0},
#endif
  {'R', 1, "results_file", "raw score file", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'T', 1, "threads", "max threads/workers", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'v', 1, "shuffle_window", "shuffle window size", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'V', 1, "annotation", "annotation characters in query/library for aligments", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'w', 1, "aln_width", "width of alignment display", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'Z', 1, "db_size", "database size for E()-value", "[library entries] database size for E()-value", 0, 0, 0, 0, 0.0, 0.0, NULL},
  {'\0', 0, "", "", NULL, 0, 0, 0, 0, 0.0, 0.0, NULL}
};

/* set default option values for help */
void g_init_opts(struct mngmsg *m_msp, struct pstruct *ppst) {
  set_opt_disp_defs('C', g_options, 1, m_msp->nmlen, 0, 0.0, 0.0, NULL);
  set_opt_disp_defs('F', g_options, 3, 0, 0, m_msp->e_low, 0.0, NULL);
  set_opt_disp_defs('m', g_options, 1, m_msp->markx, 0, 0.0, 0.0, NULL);
  set_opt_disp_defs('o', g_options, 2, (int)m_msp->sq0off, (int)m_msp->sq1off, 0.0, 0.0,NULL);
  set_opt_disp_defs('T', g_options, 1, fa_max_workers, 0, 0.0, 0.0, NULL);
  set_opt_disp_defs('v', g_options, 1, ppst->zs_win, 0, 0.0, 0.0, NULL);
  set_opt_disp_defs('w', g_options, 1, m_msp->aln.llen, 0, 0.0, 0.0, NULL);
}

void
build_optstr(char *opt_str, int opt_len, struct opt_def_str *opt_defs);

static int long_info_set=0;
static int llen_set = 0;

/* initenv ()  initializes the environment */
void initenv (int argc, char **argv, struct mngmsg *m_msp, 
		 struct pstruct *ppst, unsigned char **aa0)
{
   char *cptr, ctmp;
   long l_tmp;
   int  copt, itmp, i;

   /* options for all search functions */
   /* char   *g_optstr = "b:BC:d:DE:F:HiK:l:Lm:N:O:QqR:T:v:V:w:W:X:Z:"; */

   char g_optstring[MAX_STR];
   char f_optstring[MAX_STR];
   char optstring[MAX_STR];

   build_optstr(g_optstring, sizeof(f_optstring), g_options);
   build_optstr(f_optstring, sizeof(f_optstring), f_options);

/*  these initializations will be used by all functions */

   /* prog_name[] is only used for error messages */
   strncpy(prog_name,argv[0],sizeof(prog_name));
   prog_name[sizeof(prog_name)-1]='\0';

#ifdef PCOMPLIB
#ifdef MPI_SRC
  MPI_Comm_size(MPI_COMM_WORLD,&fa_max_workers);
  if (fa_max_workers <= 1) {
    fprintf(stderr," nnodes = %d; no workers available\n",fa_max_workers);
    exit(1);
  }
  else {
    fa_max_workers -= FIRSTNODE;
    fprintf(stderr," have %d workers\n",fa_max_workers);
  }
#endif
#else
#if defined(IRIX)
   fa_max_workers = sysmp(MP_NPROCS);
#else
#if defined(UNIX) || defined(HAVE_SYSCONF)
   fa_max_workers = sysconf(_SC_NPROCESSORS_CONF);
#endif	/* UNIX || SYSCONF */
#endif  /* !IRIX */
#endif  /* !PCOMPLIB */

   m_msp->ltitle[0] = '\0';

   if ((cptr=getenv("FASTLIBS"))!=NULL) {
     strncpy(m_msp->flstr,cptr,MAX_FN);
     m_msp->flstr[MAX_FN-1] = '\0';
   }
   else m_msp->flstr[0]='\0';

   m_msp->std_output = 1;
   m_msp->hist.hist_a = NULL;
   m_msp->outfile[0] = '\0';
   m_msp->outfd = NULL;
   m_msp->ldb_info.ldnaseq = SEQT_PROT;	/* library is protein */
   m_msp->n1_low = ppst->n1_low = 0;
   m_msp->n1_high = ppst->n1_high = BIGNUM;
   m_msp->ql_start = 1;	/* start with first query sequence */
   m_msp->ql_stop = BIGNUM;	/* end with the last query sequence */
   m_msp->aa1save_buf_b = NULL;
   m_msp->bline_buf_b = NULL;

   m_msp->pamd1 = MAXSQ;
   m_msp->pamd2 = MAXSQ;

   m_msp->ldb_info.term_code = 0;
   ppst->tr_type = 0;
   ppst->debug_lib = 0;
   m_msp->nshow = 20;
   ppst->max_repeat = 50;
   m_msp->nohist = 1;
#if defined(PCOMPLIB)
   m_msp->mshow = 20;
#else
#ifdef SHOW_HIST
   m_msp->nohist = 0;
#endif
   m_msp->mshow = 50;
#endif
   m_msp->do_showbest = 1;
   m_msp->ashow = -1;
   m_msp->ashow_set = 0;
   m_msp->nmlen = DEF_NMLEN;
   m_msp->z_bits = 1;
   m_msp->mshow_set = 0;
   m_msp->mshow_min = 0;
   m_msp->aln.llen = 60;
   m_msp->aln.llcntx = 30;
   m_msp->aln.llcntx_set = 0;
   m_msp->e_low = 0.0;
   m_msp->e_cut_set = 0;
   m_msp->revcomp = 0;
   m_msp->long_info = 0;
   m_msp->ldb_info.maxn = 0;
   m_msp->ldb_info.dupn = SEQDUP;
   m_msp->dfile[0] = '\0';
   m_msp->tname[0] = '\0';
   m_msp->lname[0] = '\0';
   m_msp->link_lname[0] = '\0';
   m_msp->show_code = 0;
   m_msp->aln.showall = 0;
   m_msp->markx = 0;
   m_msp->markx_list = NULL;
   m_msp->sq0off = m_msp->sq1off = 1;
   strncpy(m_msp->sqnam,"aa",4);
   strncpy(m_msp->sqtype,"protein",10);
   m_msp->ann_flg = 0;
   m_msp->ann_arr[0] = '\0';
   m_msp->aa0a = NULL;
   
   ppst->LK_set = 0;
   ppst->e_cut = m_msp->e_cut = 10.0;
   ppst->e_cut_r = ppst->e_cut / 10.0;
   ppst->do_rep = 1;
   ppst->zs_win = 0;
   ppst->show_ident = 0;

   ppst->zdb_size = -1;
   ppst->zdb_size_set = 0;
   ppst->dnaseq = SEQT_PROT;	/* default is protein */
   ppst->nt_align = 0;

   ppst->other_info = NULL;

   g_init_opts(m_msp, ppst);

   f_initenv (m_msp, ppst, aa0);

   if (argc == 1) {
     show_help(m_msp->pgm_name, ppst->pgm_id);
   }
   if (strcmp(argv[1],"-help")==0) {
     show_all_help(m_msp->pgm_name, ppst->pgm_id);
   }

   strncpy (optstring, g_optstring, sizeof (optstring));
   strncat (optstring, f_optstring, sizeof (optstring));

   while ((copt = getopt (argc, argv, optstring)) != EOF)
   {
      if (strchr (g_optstring, copt) != NULL)
      {
	switch (copt) {  /* switches for all options */
	case 'C': 
	  sscanf(optarg,"%d",&m_msp->nmlen);
	  if (m_msp->nmlen > MAX_UID-1) m_msp->nmlen = MAX_UID-1;
	  break;
	case 'D': ppst->debug_lib = 1;
	  break;
	case 'e': 
	  strncpy(m_msp->link_lname, optarg, MAX_LSTR);
	  break;
	case 'F':
	  sscanf(optarg,"%lg",&m_msp->e_low);
	  m_msp->e_cut_set = 1;
	  break;
#if defined(PCOMPLIB) || !defined(SHOW_HIST)
	case 'H':
	  m_msp->nohist = 0; break;
#else
	case 'H':
	  m_msp->nohist = 1; break;
#endif
	case 'i':
	  m_msp->revcomp = 1; break;
	case 'I':
	  m_msp->quiet = 0; break;
	case 'l':
	  strncpy(m_msp->flstr,optarg,MAX_FN);
	  m_msp->flstr[MAX_FN-1]='\0';
	  break;
	case 'L':
	  m_msp->long_info = 1;
	  long_info_set = 1;
	  break;
	case 'm':
	  pre_parse_markx(optarg, m_msp);
	  break;
	case 'N':
	  sscanf(optarg,"%d",&m_msp->ldb_info.maxn);
	  break;
	case 'o':
	  sscanf (optarg,"%ld %ld",&m_msp->sq0off,&m_msp->sq1off); break;
	case 'O':
	  strncpy(m_msp->outfile,optarg,MAX_FN);
	  m_msp->outfile[MAX_FN-1]='\0';
	  break;
	case 'q':
	case 'Q':
	  m_msp->quiet = 1;
	  break;
	case 'R':
	  strncpy (m_msp->dfile, optarg, MAX_FN);
	  m_msp->dfile[MAX_FN-1]='\0';
	  break;
	case 'T':
#ifdef PCOMPLIB
	  if (strchr(optarg,'-') != NULL) {
	    sscanf(optarg,"%d-%d",&worker_1,&worker_n);
	    if (worker_1 > worker_n) {
	      worker_1 = worker_n = 0;
	    }
	  }
	  else 
#endif
	    sscanf (optarg, "%d", &fa_max_workers);
	  if (fa_max_workers < 0) fa_max_workers=1;
	  break;
	case 'v':
	  sscanf (optarg,"%d",&ppst->zs_win);
	  break;
	case 'V':
	  strncpy(m_msp->ann_arr+1,optarg,MAX_FN-2);
	  m_msp->ann_arr[0]='\0';
	  m_msp->ann_arr[MAX_FN-2]='\0';
	  m_msp->ann_flg = 1;
	  add_ascii_ann(qascii, m_msp->ann_arr);
	  break;
/*
	case 'V':
	  fprintf(stderr," -V option not currently supported in parallel\n");
	  break;
*/
	case 'w':
	  sscanf (optarg,"%d",&m_msp->aln.llen);
	  if (m_msp->aln.llen < 10) m_msp->aln.llen = 10;
	  if (m_msp->aln.llen > 200) m_msp->aln.llen = 200;
	  if (!m_msp->aln.llcntx_set) m_msp->aln.llcntx = m_msp->aln.llen/2;
	  llen_set = 1;
	  break;
	case 'Z':
	  sscanf(optarg,"%ld",&ppst->zdb_size);
	  ppst->zdb_size_set = 1;
	  break;
	}
      }
      else if (strchr (f_optstring, copt))
	 f_getopt (copt, optarg, m_msp, ppst);
   }
   optind--;

   f_lastenv (m_msp, ppst);

   if (argc - optind < 3) return;
   m_msp->tnamesize = sizeof (m_msp->tname);
   if (argc - optind > 1) {strncpy (m_msp->tname, argv[optind + 1],MAX_FN);}
   if (argc - optind > 2) {strncpy(m_msp->lname, argv[optind + 2],MAX_LSTR);}
   f_getarg (argc, argv, optind, m_msp, ppst);
}

/* ann_scan scans an aa0 query sequence if -V ann_chars, and returns
   an edited query sequence and allocates aa0a[n_n0+2] space for the
   annotation */

int
ann_scan(unsigned char *aa0, int n0, unsigned char **aa0a_p, int seqtype)
{
  unsigned char *aa0p, *aa0d, *aa0ad;
  int n_n0;

  /* count how many "real" residues */

  if (seqtype==SEQT_UNK) {
    /* with SEQT_UNK, annotation characters are all < @, 
       while sequence chars are all > @ */
    for (n_n0=0, aa0p = aa0; aa0p < aa0+n0; aa0p++) {
      if (*aa0p > '@' || *aa0p == ESS ) n_n0++;
    }
  }
  else {
    /* if the sequence type is known, then annotation chars are > NANN */
    for (n_n0=0, aa0p = aa0; aa0p < aa0+n0; aa0p++) {
      if (*aa0p < NANN ) n_n0++;
    }
  }

  if (n_n0 == n0) {
    *aa0a_p = NULL;
    return n_n0;
  }

  aa0d = aa0;
  /* n_n0 has the real sequence length */
  if ((*aa0a_p = calloc(n_n0+2, sizeof(char)))==NULL) {
    fprintf(stderr," cannot allocate annotation sequence: %d\n",n_n0);
    if (seqtype==SEQT_UNK) {
      for (aa0p = aa0; aa0p < aa0+n0; aa0p++) {
	if (*aa0p > '@' || *aa0p == ESS) {*aa0d++ = *aa0p;}
      }
    }
    else {
      for (aa0p = aa0; aa0p < aa0+n0; aa0p++) {
	if (*aa0p < NANN) {*aa0d++ = *aa0p;}
      }
    }
    *aa0d = '\0';
    return n_n0;
  }

  aa0ad = *aa0a_p;
  if (seqtype==SEQT_UNK) {
    for (aa0p = aa0; aa0p<aa0+n0; aa0p++) {
      if (*aa0p > '@' || *aa0p == ESS) {*aa0d++ = *aa0p; *aa0ad++='\0';}
      else if (aa0ad > *aa0a_p) { aa0ad[-1] = *aa0p - NANN;}
    }
  }
  else {
    for (aa0p = aa0; aa0p<aa0+n0; aa0p++) {
      if (*aa0p < NANN) {*aa0d++ = *aa0p; *aa0ad++='\0';}
      else if (aa0ad > *aa0a_p) { aa0ad[-1] = *aa0p - NANN;}
    }
  }
  *aa0ad = *aa0d = '\0';
  return n_n0;
}

/* renamed from ann_ascii() Feb, 2008 to allow ann_ascii[] */
void
add_ascii_ann(int *qascii, char *ann_arr)
{
  char *ann_p;
  int ann_ix = NANN+1;

  ann_arr[0] = ' ';
  if (strchr(ann_arr+1,'*')) {qascii['*'] = NA;}

  for (ann_p = ann_arr+1; *ann_p; ann_p++) {
    if (qascii[*ann_p] == NA) { qascii[*ann_p] = ann_ix++;}
  }
}

int 
set_markx(int markx, int val) {

  if (val < 3) {
    return markx | (MX_ATYPE & val);
  }
  else if (val == 3) {
    markx |= (MX_ATYPE + MX_ASEP);
  }
  else if (val == 4) {
    markx |= (MX_ATYPE + MX_AMAP);
  }
  else if (val == 5) {
    markx |= MX_AMAP;
  }
  else if (val == 6) {
    markx |= (MX_HTML) ;
  }
  else if (val == 8) {
    markx |= MX_M9SUMM+MX_M8OUT;
  }
  else if (val == 9) {
    markx |= MX_M9SUMM;
  }
  else if (val == 10) {
    markx |= MX_M10FORM;
  }
  else if (val == 11) {
    markx |= MX_M11OUT;
  }

  return markx;
}

void
pre_parse_markx(char *opt_arg, struct mngmsg *m_msp) {
  char *bp, *last_bp;
  struct markx_str *tmp_markx, *cur_markx, *last_markx;

  if (opt_arg[0] != 'F' && m_msp->markx_list != NULL) {
    tmp_markx = m_msp->markx_list;
  }
  else {
    if ((tmp_markx = (struct markx_str *)calloc(1,sizeof(struct markx_str)))==NULL) {
      fprintf(stderr,"[error] Cannot allocate markx_list\n");
      return;
    }

    /* initialize markx to m_msg defaults -- we do not use m_msp
       directly, because it might have been changed by an earlier -m
       out_fmt */

    tmp_markx->nohist = 1;
    if (m_msp->ashow_set) {tmp_markx->ashow = m_msp->ashow;}
    else {tmp_markx->ashow = -1;}

    tmp_markx->show_code = 0;
    if (long_info_set) tmp_markx->long_info = 1;
    else tmp_markx->long_info = 0;
    if (llen_set) {
      tmp_markx->aln_llen = m_msp->aln.llen;
      tmp_markx->aln_llcntx = m_msp->aln.llcntx;
      tmp_markx->aln_llcntx_set = m_msp->aln.llcntx_set;
    }
    else {
      tmp_markx->aln_llen = 60;
      if (m_msp->aln.llcntx_set) {
	tmp_markx->aln_llcntx = m_msp->aln.llcntx;
	tmp_markx->aln_llcntx_set = m_msp->aln.llcntx_set;
      }
      else {
	tmp_markx->aln_llcntx = 30;
	tmp_markx->aln_llcntx_set = 0;
      }
    }
    tmp_markx->std_output = 1;
  }

  /* first check for -m "F file" format */
  if (optarg[0] == 'F') {
    if ((bp=strchr(optarg+1,' '))==NULL) {
      fprintf(stderr,"-m F missing file name: %s\n",optarg);
      return;
    }
    /* allocate space for file name */
    if ((tmp_markx->out_file = calloc(strlen(bp+1)+1,sizeof(char)))==NULL) {
      fprintf(stderr,"[error] Cannot allocate markx->out_file\n");
      return;
    }
    strncpy(tmp_markx->out_file, bp+1, strlen(bp+1));
    *bp = '\0';

    last_bp = optarg+1;
  }
  else {
    last_bp = optarg;
  }

  if (opt_arg[0] != 'F') {
    m_msp_to_markx(tmp_markx, m_msp);
  }

  while ((bp=strchr(last_bp,','))!=NULL) {
    *bp = '\0';
    parse_markx(last_bp, tmp_markx);
    *bp = ',';
    last_bp = bp+1;
  }

  if (*last_bp) parse_markx(last_bp, tmp_markx);

  if (m_msp->markx_list!=NULL) {
    if (opt_arg[0] == 'F') {
      /* if file name, add this to the end of the list */
      last_markx = m_msp->markx_list;
      for (cur_markx=m_msp->markx_list->next; cur_markx; cur_markx = cur_markx->next) {
	last_markx = cur_markx;
      }
      last_markx->next = tmp_markx;
    }
    else if (tmp_markx != m_msp->markx_list) {
      /* if no file name, then make this the first in the list,
	 unless it is already there */
      cur_markx = m_msp->markx_list;
      m_msp->markx_list = tmp_markx;
      tmp_markx->next = cur_markx;
    }
  }
  else {
    m_msp->markx_list = tmp_markx;
  }

  /* if no -m F, save options into m_msp */
  if (optarg[0] != 'F') {
    markx_to_m_msp(m_msp, tmp_markx);
  }

  return;
}

void
parse_markx(char *optarg, struct markx_str *this) {
  int itmp;
  char ctmp;

  itmp = 0;
  ctmp = '\0';

  if (optarg[0] == 'B') {	/* BLAST alignment output */
    this->markx = MX_MBLAST;
    this->aln_llcntx = 0;
    this->aln_llcntx_set = 1;
    this->long_info=1;
    this->ashow = -1;
    if (optarg[1] == 'B') {	/* complete BLAST output */
      this->markx += MX_MBLAST2;
      this->nohist = 1;
      this->aln_llen = 65;
      this->std_output = 0;
      return;
    }
    else if (optarg[1] == '8') {
      sscanf(optarg,"%d%c",&itmp,&ctmp);
    }
    else {return;}		/* done with BLAST aligment output */
  }
  else {
    sscanf(optarg,"%d%c",&itmp,&ctmp);
  }
  if (itmp==9) {
    if (ctmp=='c') {this->show_code = SHOW_CODE_ALIGN;}
    else if (ctmp=='C') {this->show_code = SHOW_CODE_CIGAR;}
    else if (ctmp=='i') {this->show_code = SHOW_CODE_ID;}
  }
  if (itmp > 6 && itmp != 11 && itmp != 10 && itmp != 9 && itmp != 8) itmp = 0;
  this->markx = set_markx(this->markx,itmp);
  if (itmp == 11 ) { this->std_output = 0;}
  if (itmp == 8) {
    this->std_output = 0;
    this->ashow = 0;
    if (ctmp=='C') { this->markx += MX_M8COMMENT;}
  }
}

/* transfer markx values for m_msp to m_msp */
void
markx_to_m_msp(struct mngmsg *m_msp, struct markx_str *this) {

  m_msp->markx = this->markx;
  m_msp->nohist = this->nohist;
  m_msp->ashow = this->ashow;
  m_msp->show_code = this->show_code;
  m_msp->long_info = this->long_info;
  m_msp->aln.llen = this->aln_llen;
  m_msp->aln.llcntx = this->aln_llcntx;
  m_msp->aln.llcntx_set = this->aln_llcntx_set;
  m_msp->std_output = this->std_output;
}

/* save current m_msp values used with markx */
void
m_msp_to_markx(struct markx_str *this, struct mngmsg *m_msp) {

  this->markx = m_msp->markx ;
  this->nohist = m_msp->nohist ;
  this->ashow = m_msp->ashow ;
  this->show_code = m_msp->show_code ;
  this->long_info = m_msp->long_info ;
  this->aln_llen = m_msp->aln.llen ;
  this->aln_llcntx = m_msp->aln.llcntx ;
  this->aln_llcntx_set = m_msp->aln.llcntx_set ;
  this->std_output = m_msp->std_output ;
}

/* put options from option table [struct opt_def_str *opt_defs] into
   char *opt_str for getopt() */

void
build_optstr(char *opt_str, int max_len, struct opt_def_str *opt_defs) {
  int i, opt_len = 0;
  char *opt_pos;

  opt_pos = opt_str;
  for (i=0; opt_defs[i].opt_char != '\0'; i++) {
    if (opt_len + 2 > max_len) {
      fprintf(stderr," *** error -- options too long %d >= %d\n", opt_len, max_len);
      break;
    }
    *opt_pos++ = opt_defs[i].opt_char;
    opt_len++;
    if (opt_defs[i].has_arg) {
      *opt_pos++ = ':';
      opt_len++;
    }
  }
  *opt_pos = '\0';
}

/* set_opt_disp_defs associates parameter addresses with options */
void
set_opt_disp_defs(char opt_char, struct opt_def_str *options, int type,
		  int i_param1, int i_param2,
		  double d_param1, double d_param2,
		  char *s_param) {
  struct opt_def_str *this_opt;

  this_opt = options;
  while (this_opt->opt_char != '\0') {
    if (this_opt->opt_char == opt_char) {
      this_opt->fmt_type = type;
      switch (type) {
      case 1:
	this_opt->i_param1 = i_param1;
	break;
      case 2:
	this_opt->i_param1 = i_param1;
	this_opt->i_param2 = i_param2;
	break;
      case 3:
	this_opt->d_param1 = d_param1;
	break;
      case 4:
	this_opt->d_param1 = d_param1;
	this_opt->d_param2 = d_param2;
	break;
      case 5:
	if (s_param != NULL) {
	  this_opt->s_param = (char *)calloc(strlen(s_param)+1,sizeof(char));
	  strncpy(this_opt->s_param,s_param,strlen(s_param));
	}
	else this_opt->s_param = NULL;
	break;
      }
    }
    this_opt++;
  }
}
