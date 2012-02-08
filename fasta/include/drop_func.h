/* drop_func.h */

/* $Id: drop_func.h 625 2011-03-23 17:21:38Z wrp $ */
/* $Revision: 625 $  */

/* functions provided by each of the drop files */

/* Copyright (c) 2005 William R. Pearson and the University of Virginia */


#ifdef DEBUG
unsigned long adler32(unsigned long, const unsigned char *, unsigned int);
#endif

void	/* initializes f_struct **f_arg */
init_work (unsigned char *aa0, int n0,
	   struct pstruct *ppst,
#ifndef DROP_INTERN
	   void **f_arg
#else
	   struct f_struct **f_arg
#endif
);


void	/* frees memory allocated in f_struct */
close_work (const unsigned char *aa0, int n0,
	    struct pstruct *ppst,
#ifndef DROP_INTERN
	   void **f_arg
#else
	   struct f_struct **f_arg
#endif
);

void	/* documents search function, parameters */
get_param (const struct pstruct *pstr, 
	   char **pstring1, char *pstring2,
	   struct score_count_s *);

void	/* calculates alignment score(s), returns them in rst */
do_work (const unsigned char *aa0, int n0,
	 const unsigned char *aa1, int n1,
	 int frame,
	 const struct pstruct *ppst,
#ifndef DROP_INTERN
	 void *f_arg,
#else
	 struct f_struct *f_arg,
#endif
	 int qr_flg, int shuff_flg, struct rstruct *rst,
	 struct score_count_s *);

void	/* calculates optimal alignment score */
do_opt (const unsigned char *aa0, int n0,
	const unsigned char *aa1, int n1,
	int frame,
	struct pstruct *ppst,
#ifndef DROP_INTERN
	void *f_arg,
#else
	struct f_struct *f_arg,
#endif
	struct rstruct *rst
	);

struct a_res_str *	/* produces encoding of alignment */
do_walign (const unsigned char *aa0, int n0,
	   const unsigned char *aa1, int n1,
	   int frame, int repeat_thresh,
	   struct pstruct *ppst, 
#ifndef DROP_INTERN
	   void *f_arg,
#else
	   struct f_struct *f_arg,
#endif
	   int *have_ares);

void
pre_cons(const unsigned char *aa, int n, int frame, 
#ifndef DROP_INTERN
	   void *f_arg
#else
	   struct f_struct *f_arg
#endif
	);

void 
aln_func_vals(int frame, struct a_struct *aln);

/* calc_cons_a - takes aa0, aa1, a_res, and produces seqc0, seqc1, 
 *             and seqc0a, seqc1a - the annotated sequences 
 */
int
calc_cons_a(const unsigned char *aa0, int n0,
	  const unsigned char *aa1, int n1,
	  int *nc,
	  struct a_struct *aln,
	  struct a_res_str *a_res,
	  const struct pstruct *ppst,
	  char *seqc0, char *seqc1, char *seqca,
	  const char *ann_arr,
	  const unsigned char *aa0a, char *seqc0a, 
	  const unsigned char *aa1a, char *seqc1a, 
#ifndef DROP_INTERN
	  void *f_arg
#else
	  struct f_struct *f_arg
#endif
	  );

int	/* returns lenc - length of aligment */
calc_code(const unsigned char *aa0, int n0,
	  const unsigned char *aa1, int n1,
	  struct a_struct *aln,
	  struct a_res_str *a_res,
	  struct pstruct *ppst,
	  char *al_str, int al_str_n,
	  const char *ann_arr,
	  const unsigned char *aa0a, const unsigned char *aa1a, 
	  char *ann_str, int ann_str_n,
#ifndef DROP_INTERN
	  void *f_arg,
#else
	  struct f_struct *f_arg,
#endif
	  int code_fmt
	  );

int 	/* returns lenc - length of alignment */
calc_id(const unsigned char *aa0, int n0,
	unsigned char *aa1, int n1,
	struct a_struct *aln, 
	struct a_res_str *a_res,
	struct pstruct *ppst,
#ifndef DROP_INTERN
	void *f_arg
#else
	struct f_struct *f_arg
#endif
	);
