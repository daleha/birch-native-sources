
-- init libraries... it is really messy that they only differ by defines...

	project "init_sw"
		language    "C"
		kind        "StaticLib"
		targetdir PROJDIR.."lib"

		defines
		{
			"SSEARCH"
		}

		files
		{
			SRCDIR.."initfa.c"
			,SRCDIR.."doinit.c"
		}
		--init_sw_sse.o : initfa.c -DSW_SSE2 -DSSEARCH 
		--init_sw_alt.o : initfa.c -DSW_ALTIVEC -DSSEARCH 

	project "init_lal"

		language    "C"
		kind        "StaticLib"
		targetdir PROJDIR.."lib"

		defines
		{
			"LALIGN"
		}

		files
		{
			SRCDIR.."initfa.c"
			,SRCDIR.."doinit.c"
		}

	project "init_lnw"

		language    "C"
		kind        "StaticLib"
		targetdir PROJDIR.."lib"

		defines
		{
			"GLSEARCH"
		}

		files
		{
			SRCDIR.."initfa.c"
			,SRCDIR.."doinit.c"
		}

		--init_lnw_sse.o : initfa.c -DSW_SSE2 -DGLSEARCH 

	project "init_gnw"

		language    "C"
		kind        "StaticLib"
		targetdir PROJDIR.."lib"

		defines
		{
			"GGSEARCH"
		}

		files
		{
			SRCDIR.."initfa.c"
			,SRCDIR.."doinit.c"
		}
		--init_gnw_sse.o : initfa.c  -DSW_SSE2 -DGGSEARCH 


	project "init_gnw"

		language    "C"
		kind        "StaticLib"
		targetdir PROJDIR.."lib"

		defines
		{
			"GGSEARCH"
		}

		files
		{
			SRCDIR.."initfa.c"
			,SRCDIR.."doinit.c"
		}

	project "init_rss"

		language    "C"
		kind        "StaticLib"
		targetdir PROJDIR.."lib"

		defines
		{
			"PRSS"
		}

		files
		{
			SRCDIR.."initfa.c"
			,SRCDIR.."doinit.c"
		}


	project "init_rfx"

		language    "C"
		kind        "StaticLib"
		targetdir PROJDIR.."lib"

		defines
		{
			"PRSS"
			,"FASTX"
		}

		files
		{
			SRCDIR.."initfa.c"
			,SRCDIR.."doinit.c"
		}

	project "init_fa"

		language    "C"
		kind        "StaticLib"
		targetdir PROJDIR.."lib"

		defines
		{
			"FASTA"
		}

		files
		{
			SRCDIR.."initfa.c"
			,SRCDIR.."doinit.c"
		}


	project "init_ff"

		language    "C"
		kind        "StaticLib"
		targetdir PROJDIR.."lib"

		defines
		{
			"FASTF"
		}

		files
		{
			SRCDIR.."initfa.c"
			,SRCDIR.."doinit.c"
		}


	project "init_tf"

		language    "C"
		kind        "StaticLib"
		targetdir PROJDIR.."lib"

		defines
		{
			"FASTF"
			,"TFAST"
		}

		files
		{
			SRCDIR.."initfa.c"
			,SRCDIR.."doinit.c"
		}

	project "init_fs"

		language    "C"
		kind        "StaticLib"
		targetdir PROJDIR.."lib"

		defines
		{
			"FASTS"
		}

		files
		{
			SRCDIR.."initfa.c"
			,SRCDIR.."doinit.c"
		}

	project "init_fm"

		language    "C"
		kind        "StaticLib"
		targetdir PROJDIR.."lib"

		defines
		{
			"FASTM"
		}

		files
		{
			SRCDIR.."initfa.c"
			,SRCDIR.."doinit.c"
		}

	project "init_tfs"

		language    "C"
		kind        "StaticLib"
		targetdir PROJDIR.."lib"

		defines
		{
			"FASTS"
			,"TFAST"
		}

		files
		{
			SRCDIR.."initfa.c"
			,SRCDIR.."doinit.c"
		}

	project "init_tfm"

		language    "C"
		kind        "StaticLib"
		targetdir PROJDIR.."lib"

		defines
		{
			"FASTM"
			,"TFAST"
		}

		files
		{
			SRCDIR.."initfa.c"
			,SRCDIR.."doinit.c"
		}

	project "init_tfa"

		language    "C"
		kind        "StaticLib"
		targetdir PROJDIR.."lib"

		defines
		{
			"FASTA"
			,"TFAST"
		}

		files
		{
			SRCDIR.."initfa.c"
			,SRCDIR.."doinit.c"
		}

	project "init_fx"

		language    "C"
		kind        "StaticLib"
		targetdir PROJDIR.."lib"

		defines
		{
			"FASTX"
		}

		files
		{
			SRCDIR.."initfa.c"
			,SRCDIR.."doinit.c"
		}


	project "init_tfx"

		language    "C"
		kind        "StaticLib"
		targetdir PROJDIR.."lib"

		defines
		{
			"FASTX"
			,"TFAST"
		}

		files
		{
			SRCDIR.."initfa.c"
			,SRCDIR.."doinit.c"
		}



	project "init_fy"

		language    "C"
		kind        "StaticLib"
		targetdir PROJDIR.."lib"

		defines
		{
			"FASTY"
		}

		files
		{
			SRCDIR.."initfa.c"
			,SRCDIR.."doinit.c"
		}

	project "init_tfy"

		language    "C"
		kind        "StaticLib"
		targetdir PROJDIR.."lib"

		defines
		{
			"FASTY"
			,"TFAST"
		}

		files
		{
			SRCDIR.."initfa.c"
			,SRCDIR.."doinit.c"
		}




#================ miscellaneous

htime.o : htime.c

compacc2_t.o : compacc2.c  -DCOMP_THR -DCOMP_MLIB 

compacc2_s.o : compacc2.c 

compacc2_p.o : compacc2.c  -DMPI_SRC 

compacc.o : compacc.c -DCOMP_THR 

pssm_asn_subs.o : pssm_asn_subs.c 


#================ display list of best hits / alignments

showbest.o : $(SHOWBESTC) 

build_ares.o : build_ares.c 

$(SHOWALIGN_T).o : $(SHOWALIGN).c  -DCOMP_THR 

$(SHOWALIGN_P).o : $(SHOWALIGN).c -DMPI_SRC 

$(SHOWALIGN_S).o : $(SHOWALIGN).c 

$(LSHOWALIGN).o : $(SHOWALIGN).c  -DLALIGN 

re_getlib.o : re_getlib.c 

lib_sel.o : lib_sel.c

c_dispn.o : c_dispn.c 

#================ statistical functions

karlin.o : karlin.c 

scale_se.o : scaleswn.c -DLOCAL_SCORE 

scale_sn.o : scaleswn.c -DNORMAL_DIST 

scaleswtf.o : scaleswt.c  -DFASTF 

scaleswts.o : scaleswt.c 

tatstats_fs.o : tatstats.c  -DFASTS

tatstats_ff.o : tatstats.c  -DFASTF 

tatstats_fm.o : tatstats.c  -DFASTM 

last_tat.o : last_tat.c 

last_thresh.o : last_thresh.c 


#================ drop functions - actual scores/alignments

drop_nfa.o : dropnfa.c 

dropsbd.o : dropnfa.c 

# drop_ff, _fs, _fm must define FASTF, FASTS, and FASTM to ensure
# that tatstats.h is built appropriately

drop_ff2.o : dropff2.c  -DFASTF  

drop_tff.o : dropff2.c  -DFASTF -DTFAST 

drop_fs2.o : dropfs2.c  -DFASTS 

drop_tfs.o : dropfs2.c  -DTFAST -DFASTS 

drop_fm.o : dropfs2.c  -DFASTM 

drop_tfm.o : dropfs2.c  -DTFAST -DFASTM

drop_tfa.o : dropnfa.c  -DTFASTA 

drop_fx.o : dropfx.c 

drop_tfx.o : dropfx.c  -DTFAST 

drop_fz.o : dropfz2.c 

drop_tfz.o : dropfz2.c  -DTFAST 

dropnsw.o : dropnsw.c 

dropgsw2.o : dropgsw2.c 

dropgsw2_sse.o : dropgsw2.c -DSW_SSE2 

dropgsw2_alt.o : dropgsw2.c  -DSW_ALTIVEC 

droplal2.o : dropgsw2.c -DLALIGN 

droplal2_sse.o : dropgsw2.c -DLALIGN -DSW_SSE2 

droplal2_alt.o : -DLALIGN  -DSW_ALTIVEC 

lsim4.o : lsim4.c 

smith_waterman_altivec.o : smith_waterman_altivec.c -DSW_ALTIVEC 

smith_waterman_sse2.o : smith_waterman_sse2.c -DSW_SSE2 

global_sse2.o : global_sse2.c -DSW_SSE2 

glocal_sse2.o : glocal_sse2.c  -DSW_SSE2 

droplnw.o : dropnnw2.c 

droplnw_sse.o : -DSW_SSE2 

dropgnw.o : dropnnw2.c -DGLOBAL_GLOBAL 

dropgnw_sse.o : dropnnw2.c -DGLOBAL_GLOBAL -DSW_SSE2 

wm_align.o : wm_align.c 

calcons_fa.o : cal_cons.c  -DFASTA 

calcons_tfa.o : cal_cons.c  -DTFASTA 

calcons_sw.o : cal_cons.c  -DSSEARCH 

calcons_la.o : cal_cons.c -DLALIGN -DLCAL_CONS 

calcons_ff.o : cal_consf.c -DFASTF 

calcons_fs.o : cal_consf.c -DFASTS 

calcons_fm.o : cal_consf.c -DFASTM 

calcons_tff.o : cal_consf.c  -DTFAST -DFASTF 

calcons_tfs.o : cal_consf.c  -DTFAST -DFASTS 

calcons_tfm.o : cal_consf.c  -DTFAST -DFASTM 

#================ reading query, libraries

getseq.o : getseq.c 

llgetaa.o : llgetaa.c  -DNOLIB 

lgetlib.o : $(NGETLIB).c 

lgetaa_m.o : mmgetaa.c 

ncbl_lib.o : ncbl_lib.c 

ncbl2_mlib.o : ncbl2_mlib.c 

mysql_lib.o : mysql_lib.c 

pgsql_lib.o : pgsql_lib.c 

#================ threading functions

pthr_subs2.o : pthr_subs2.c 

uthr_subs.o : uthr_subs.c 

#================ MPI worker function

mpi_subs2.o : pcomp_subs2.c -DPCOMPLIB -DMPI_SRC 

#================ translation

faatran.o : faatran.c 

url_subs.o : url_subs.c 

#================ lav plotting functions

lav2plt.o : lav2plt.c 

lavplt_ps.o : lavplt_ps.c 

lavplt_svg.o : lavplt_svg.c 


#====== more libs

comp_mlib5e.o : comp_lib5e.c  -DCOMP_MLIB 

comp_mthr5e.o : comp_lib5e.c  -DCOMP_THR -DCOMP_MLIB 

comp_mlib7e.o : comp_lib7e.c  -DCOMP_MLIB 

comp_mthr7e.o : comp_lib7e.c  -DCOMP_THR -DCOMP_MLIB 

comp_mlib8.o : comp_lib8.c  -DCOMP_MLIB 

comp_mthr8.o : comp_lib8.c  -DCOMP_THR -DCOMP_MLIB 

work_thr2.o : work_thr2.c 

