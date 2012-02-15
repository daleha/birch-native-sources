

PROJDIR=ROOTDIR
SRCDIR=PROJDIR.."source/"
INCDIR=PROJDIR.."include/"


-- common configuration options
	configuration "Debug"
		buildoptions
		{
		--	"-O3"

		}

		defines
		{
			"SHOW_HELP" 
			,"SHOWSIM" 
			,"HZ=100" 
			,"MAX_WORKERS=4"
			,"THR_EXIT=pthread_exit" 
			,"PROGRESS" 
			,"FASTA_HOST='\"your_fasta_host_here\"'" 
			,"_REENTRANT" 
			,"HAS_INTTYPES" 
			,"_LARGEFILE_SOURCE" 
			,"_LARGEFILE64_SOURCE" 
			,"_FILE_OFFSET_BITS=64" 
			,"SAMP_STATS"
			,"PGM_DOC" 
			,"BIG_LIB64"

		}

		includedirs
		{
			INCDIR --add the header folder to include search path
		}

		links
		{
			"m" --math libraries
			,"z"
		}

		--dofile "./phylip4_lib.lua"

	configuration "not windows"
		defines
		{
			"UNIX" 
			,"TIMES" 
			,"USE_MMAP" 
			,"USE_FSEEKO"

		}

	configuration "windows"
		defines
		{
			"WIN32"
		}

		libdirs
		{
			"/usr/local/lib" -- used for zlib, which must be build separately
		}



	configuration "linux"

		buildoptions
		{
		--	"-fomit-frame-pointer" --optimizations
		}


	configuration "solaris"

		buildoptions
		{
	--		"-fomit-frame-pointer" --optimizations
		}


	configuration "macosx"

		buildoptions
		{
	--		"-Wall -mmacosx-version-min=10.1"	
		}
		linkoptions
		{
	--		"-Wl"	
		}




--# combinations of files for "composite" drop* functions
--#
--DROPNSW_O = dropnsw.o  wm_align.o calcons_sw.o
--DROPNFA_O = drop_nfa.o wm_align.o calcons_fa.o
--DROPBD_O = dropsbd.o wm_align.o calcons_fa.o
--DROPTFA_O = drop_tfa.o wm_align.o calcons_tfa.o
--DROPFF_O = drop_ff2.o calcons_ff.o
--DROPFS_O = drop_fs2.o calcons_fs.o
--DROPFM_O = drop_fm.o calcons_fm.o
--DROPTFF_O = drop_tff.o calcons_tff.o
--DROPTFS_O = drop_tfs.o calcons_tfs.o
--DROPTFM_O = drop_tfm.o calcons_tfm.o
--
--#COMPACC_TO = compacc.o	# used with comp_lib5.c/comp_lib7.c
--#COMPACC_SO = compacc.o
--COMPACC_TO = compacc2_t.o  # used with comp_lib5e.c/comp_lib7e.c/comp_lib8.c
--COMPACC_SO = compacc2_s.o
--
--SHOWBESTC = mshowbest.c
--SHOWBESTO = showbest.o build_ares.o
--SHOWALIGN = mshowalign2
--SHOWALIGN_T = mshowalign2_t
--SHOWALIGN_S = mshowalign2_s
--LSHOWALIGN = lshowalign
--MWH = mw.h 
--MWHP = mw.h
--
--TPROGS = ssearch36_t fasta36_t  fasts36_t fastx36_t tfastx36_t fasty36_t tfasty36_t tfasts36_t fastm36_t fastf36_t tfastf36_t glsearch36_t ggsearch36_t
--
--SPROGS = fasta36 ssearch36 lalign36 fasts36 fastx36 tfastx36 fasty36 tfasty36 tfasts36 fastm36 tfastm36 fastf36 tfastf36 glsearch36 ggsearch36 lav2ps lav2svg
--
--APROGS = map_db
--
--XTPROGS = fastx36_t tfastx36_t fasty36_t tfasty36_t
--XPROGS = fastx36 tfastx36  fasty36 tfasty36
--
--PROGS = $(SPROGS) $(TPROGS) $(APROGS)
--
--xall: $(XTPROGS) $(XPROGS) $(ZTPROGS) $(ZPROGS) 


project "fasta36" 
	language    "C"
	kind        "ConsoleApp"

	files
	{

		SRCDIR.."re_getlib.c"
		,SRCDIR.."htime.c"
		,SRCDIR.."apam.c"
		,SRCDIR.."initfa.c"
		,SRCDIR.."doinit.c"
		,SRCDIR.."scaleswn.c" -- scalese DLOCAL_SCORE --
		,SRCDIR.."karlin.c"
		,SRCDIR.."dropnfa.c"
		,SRCDIR.."wm_align.c"
		,SRCDIR.."cal_cons.c"
		,SRCDIR.."c_dispn.c"
		,SRCDIR.."lib_sel.c"
		,SRCDIR.."mrandom.c"
		,SRCDIR.."url_subs.c"
		,SRCDIR.."pssm_asn_subs.c"

		,SRCDIR.."comp_lib8.c"
		,SRCDIR.."compacc2.c"
		,SRCDIR.."mshowbest.c"
		,SRCDIR.."build_ares.c"
		,SRCDIR.."mshowalign2.c"
		,SRCDIR.."dropnfa.c"
		,SRCDIR.."nmgetlib.c"
		,SRCDIR.."ncbl2_mlib.c"


	}

	defines
	{
		"FASTA"
		,"LOCAL_SCORE"
		,"COMP_MLIB"
		--,"LALIGN"
		--,"LCAL_CONS"
	}


	configuration "not windows"
		files
		{
			SRCDIR.."mmgetaa.c"
			,SRCDIR.."getseq.c"
		}

		defines
		{
		}

	configuration "windows"

		files
		{
			SRCDIR.."getopt.c"
		}

	 
project "fastx36" 
       language    "C"
       kind        "ConsoleApp"

       files
       {

               SRCDIR.."re_getlib.c"
               ,SRCDIR.."htime.c"
               ,SRCDIR.."apam.c"
               ,SRCDIR.."initfa.c"
               ,SRCDIR.."doinit.c"
               ,SRCDIR.."scaleswn.c" -- scalese DLOCAL_SCORE
               ,SRCDIR.."karlin.c"
               ,SRCDIR.."dropfx.c"
               ,SRCDIR.."c_dispn.c"
               ,SRCDIR.."lib_sel.c"
               ,SRCDIR.."mrandom.c"
               ,SRCDIR.."url_subs.c"
               ,SRCDIR.."pssm_asn_subs.c"
               ,SRCDIR.."faatran.c"
        ,SRCDIR.."comp_lib8.c"
		,SRCDIR.."compacc2.c"
		,SRCDIR.."mshowbest.c"
		,SRCDIR.."build_ares.c"
		,SRCDIR.."mshowalign2.c"
		,SRCDIR.."nmgetlib.c"
		,SRCDIR.."ncbl2_mlib.c"



       }

       defines
       {
               "FASTX"
               ,"LOCAL_SCORE"
               ,"COMP_MLIB"
               --,"LALIGN"
               --,"LCAL_CONS"
       }


	configuration "not windows"
		files
		{
			SRCDIR.."mmgetaa.c"
			,SRCDIR.."getseq.c"
		}


	configuration "windows"

		files
		{
			SRCDIR.."getopt.c"
		}



        




 

 

project "fasty36" 
	language    "C"
	kind        "ConsoleApp"

	files
	{

		SRCDIR.."re_getlib.c"
		,SRCDIR.."htime.c"
		,SRCDIR.."apam.c"
		,SRCDIR.."initfa.c"
		,SRCDIR.."doinit.c"
		,SRCDIR.."scaleswn.c" -- scalese DLOCAL_SCORE --
		,SRCDIR.."karlin.c"
		,SRCDIR.."dropfz2.c"
		,SRCDIR.."c_dispn.c"
		,SRCDIR.."lib_sel.c"
		,SRCDIR.."mrandom.c"
		,SRCDIR.."url_subs.c"
		,SRCDIR.."pssm_asn_subs.c"
		,SRCDIR.."faatran.c"
,SRCDIR.."comp_lib8.c"
		,SRCDIR.."compacc2.c"
		,SRCDIR.."mshowbest.c"
		,SRCDIR.."build_ares.c"
		,SRCDIR.."mshowalign2.c"
		,SRCDIR.."nmgetlib.c"
		,SRCDIR.."ncbl2_mlib.c"


	}

	defines
	{
		"FASTY"
		,"LOCAL_SCORE"
        ,"COMP_MLIB"
		--,"LALIGN"
		--,"LCAL_CONS"
	}


	configuration "not windows"
		files
		{
			SRCDIR.."mmgetaa.c"
			,SRCDIR.."getseq.c"
		}


	configuration "windows"

		files
		{
			SRCDIR.."getopt.c"
		}





project "scaleswts"	 
	language "C"
	kind "StaticLib"

	files
	{

		SRCDIR.."scaleswt.c" 
	}

project "fastf36" 
	language    "C"
	kind        "ConsoleApp"

	files
	{

		SRCDIR.."re_getlib.c"
		,SRCDIR.."htime.c"
		,SRCDIR.."apam.c"
		,SRCDIR.."initfa.c"
		,SRCDIR.."doinit.c"
		,SRCDIR.."karlin.c"
		,SRCDIR.."dropff2.c"
		,SRCDIR.."last_tat.c"
		,SRCDIR.."tatstats.c"
		,SRCDIR.."c_dispn.c"
		,SRCDIR.."lib_sel.c"
		,SRCDIR.."mmgetaa.c"
		,SRCDIR.."cal_consf.c"
		,SRCDIR.."mrandom.c"
		,SRCDIR.."url_subs.c"
		,SRCDIR.."pssm_asn_subs.c"
,SRCDIR.."comp_lib8.c"
		,SRCDIR.."compacc2.c"
		,SRCDIR.."mshowbest.c"
		,SRCDIR.."build_ares.c"
		,SRCDIR.."mshowalign2.c"
		,SRCDIR.."nmgetlib.c"
		,SRCDIR.."ncbl2_mlib.c"


	}
	links {"scaleswts" }

	defines
	{
		"FASTF"
        ,"COMP_MLIB"
		--,"LALIGN"
		--,"LCAL_CONS"
	}



	configuration "not windows"
		files
		{
			SRCDIR.."mmgetaa.c"
			,SRCDIR.."getseq.c"
		}


	configuration "windows"

		files
		{
			SRCDIR.."getopt.c"
		}



	
project "fasts36" 
	language    "C"
	kind        "ConsoleApp"

	files
	{

		SRCDIR.."re_getlib.c"
		,SRCDIR.."htime.c"
		,SRCDIR.."apam.c"
		,SRCDIR.."initfa.c"
		,SRCDIR.."doinit.c"
		,SRCDIR.."karlin.c"
		,SRCDIR.."dropfs2.c"
		,SRCDIR.."last_tat.c"
		,SRCDIR.."tatstats.c"
		,SRCDIR.."c_dispn.c"
		,SRCDIR.."lib_sel.c"
		,SRCDIR.."mmgetaa.c"
		,SRCDIR.."cal_consf.c"
		,SRCDIR.."mrandom.c"
		,SRCDIR.."url_subs.c"
		,SRCDIR.."pssm_asn_subs.c"
,SRCDIR.."comp_lib8.c"
		,SRCDIR.."compacc2.c"
		,SRCDIR.."mshowbest.c"
		,SRCDIR.."build_ares.c"
		,SRCDIR.."mshowalign2.c"
		,SRCDIR.."nmgetlib.c"
		,SRCDIR.."ncbl2_mlib.c"


	}
	links {"scaleswts" }

	defines
	{
		"FASTS"
        ,"COMP_MLIB"
		--,"LALIGN"
		--,"LCAL_CONS"
	}

	configuration "not windows"
		files
		{
			SRCDIR.."mmgetaa.c"
			,SRCDIR.."getseq.c"
		}


	configuration "windows"

		files
		{
			SRCDIR.."getopt.c"
		}


	
project "fastm36" 
	language    "C"
	kind        "ConsoleApp"

	files
	{

		SRCDIR.."re_getlib.c"
		,SRCDIR.."htime.c"
		,SRCDIR.."apam.c"
		,SRCDIR.."initfa.c"
		,SRCDIR.."doinit.c"
		,SRCDIR.."karlin.c"
		,SRCDIR.."dropfs2.c"
		,SRCDIR.."last_tat.c"
		,SRCDIR.."tatstats.c"
		,SRCDIR.."c_dispn.c"
		,SRCDIR.."lib_sel.c"
		,SRCDIR.."mmgetaa.c"
		,SRCDIR.."cal_consf.c"
		,SRCDIR.."mrandom.c"
		,SRCDIR.."url_subs.c"
		,SRCDIR.."pssm_asn_subs.c"
,SRCDIR.."comp_lib8.c"
		,SRCDIR.."compacc2.c"
		,SRCDIR.."mshowbest.c"
		,SRCDIR.."build_ares.c"
		,SRCDIR.."mshowalign2.c"
		,SRCDIR.."nmgetlib.c"
		,SRCDIR.."ncbl2_mlib.c"


	}
	links {"scaleswts" }

	defines
	{
		"FASTM"
        ,"COMP_MLIB"
		--,"LALIGN"
		--,"LCAL_CONS"
	}

	configuration "not windows"
		files
		{
			SRCDIR.."mmgetaa.c"
			,SRCDIR.."getseq.c"
		}


	configuration "windows"

		files
		{
			SRCDIR.."getopt.c"
		}


project "tfastx36" 
	language    "C"
	kind        "ConsoleApp"

	files
	{

		SRCDIR.."re_getlib.c"
		,SRCDIR.."htime.c"
		,SRCDIR.."apam.c"
		,SRCDIR.."initfa.c"
		,SRCDIR.."doinit.c"
		,SRCDIR.."scaleswn.c" -- scalese DLOCAL_SCORE --
		,SRCDIR.."karlin.c"
		,SRCDIR.."dropfx.c"
		,SRCDIR.."c_dispn.c"
		,SRCDIR.."lib_sel.c"
		,SRCDIR.."mrandom.c"
		,SRCDIR.."url_subs.c"
		,SRCDIR.."pssm_asn_subs.c"
		,SRCDIR.."faatran.c"
,SRCDIR.."comp_lib8.c"
		,SRCDIR.."compacc2.c"
		,SRCDIR.."mshowbest.c"
		,SRCDIR.."build_ares.c"
		,SRCDIR.."mshowalign2.c"
		,SRCDIR.."nmgetlib.c"
		,SRCDIR.."ncbl2_mlib.c"


	}

	defines
	{
		"FASTX"
		,"TFAST"
		,"LOCAL_SCORE"
        ,"COMP_MLIB"
		--,"LALIGN"
		--,"LCAL_CONS"
	}

	configuration "not windows"
		files
		{
			SRCDIR.."mmgetaa.c"
			,SRCDIR.."getseq.c"
		}


	configuration "windows"

		files
		{
			SRCDIR.."getopt.c"
		}
	
project "tfasts36" 
	language    "C"
	kind        "ConsoleApp"

	files
	{

		SRCDIR.."re_getlib.c"
		,SRCDIR.."htime.c"
		,SRCDIR.."apam.c"
		,SRCDIR.."initfa.c"
		,SRCDIR.."doinit.c"
		,SRCDIR.."karlin.c"
		,SRCDIR.."dropfs2.c"
		,SRCDIR.."last_tat.c"
		,SRCDIR.."tatstats.c"
		,SRCDIR.."c_dispn.c"
		,SRCDIR.."lib_sel.c"
		,SRCDIR.."mmgetaa.c"
		,SRCDIR.."cal_consf.c"
		,SRCDIR.."mrandom.c"
		,SRCDIR.."url_subs.c"
		,SRCDIR.."faatran.c"
		,SRCDIR.."pssm_asn_subs.c"
,SRCDIR.."comp_lib8.c"
		,SRCDIR.."compacc2.c"
		,SRCDIR.."mshowbest.c"
		,SRCDIR.."build_ares.c"
		,SRCDIR.."mshowalign2.c"
		,SRCDIR.."nmgetlib.c"
		,SRCDIR.."ncbl2_mlib.c"


	}
	links {"scaleswts" }

	defines
	{
		"FASTS"
		,"TFAST"
        ,"COMP_MLIB"
		--,"LALIGN"
		--,"LCAL_CONS"
	}


	configuration "not windows"
		files
		{
			SRCDIR.."mmgetaa.c"
			,SRCDIR.."getseq.c"
		}


	configuration "windows"

		files
		{
			SRCDIR.."getopt.c"
		}

--tfasts36 : $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o init_tfs.o scaleswts.o tatstats_fs.o last_tat.o karlin.o $(DROPTFS_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o faatran.o mrandom.o url_subs.o
--	$(CC) $(HFLAGS) $(BIN)/tfasts36 $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o init_tfs.o $(DROPTFS_O) scaleswts.o tatstats_fs.o last_tat.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o faatran.o mrandom.o url_subs.o $(LIB_M)
--



--fasts36 : $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o init_fs.o scaleswts.o last_tat.o tatstats_fs.o karlin.o $(DROPFS_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o mrandom.o url_subs.o
--	$(CC) $(HFLAGS) $(BIN)/fasts36 $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o init_fs.o $(DROPFS_O) scaleswts.o last_tat.o tatstats_fs.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o mrandom.o url_subs.o $(LIB_M)


project "tfasty36" 
	language    "C"
	kind        "ConsoleApp"

	files
	{

		SRCDIR.."re_getlib.c"
		,SRCDIR.."htime.c"
		,SRCDIR.."apam.c"
		,SRCDIR.."initfa.c"
		,SRCDIR.."doinit.c"
		,SRCDIR.."scaleswn.c" -- scalese DLOCAL_SCORE --
		,SRCDIR.."karlin.c"
		,SRCDIR.."dropfz2.c"
		,SRCDIR.."c_dispn.c"
		,SRCDIR.."lib_sel.c"
		,SRCDIR.."mrandom.c"
		,SRCDIR.."url_subs.c"
		,SRCDIR.."pssm_asn_subs.c"
		,SRCDIR.."faatran.c"
,SRCDIR.."comp_lib8.c"
		,SRCDIR.."compacc2.c"
		,SRCDIR.."mshowbest.c"
		,SRCDIR.."build_ares.c"
		,SRCDIR.."mshowalign2.c"
		,SRCDIR.."nmgetlib.c"
		,SRCDIR.."ncbl2_mlib.c"


	}

	defines
	{
		"FASTY"
		,"TFAST"
		,"LOCAL_SCORE"
        ,"COMP_MLIB"
		--,"LALIGN"
		--,"LCAL_CONS"
	}

	configuration "not windows"
		files
		{
			SRCDIR.."mmgetaa.c"
			,SRCDIR.."getseq.c"
		}


	configuration "windows"

		files
		{
			SRCDIR.."getopt.c"
		}



project "tfastf36" 
	language    "C"
	kind        "ConsoleApp"

	files
	{

		SRCDIR.."re_getlib.c"
		,SRCDIR.."htime.c"
		,SRCDIR.."apam.c"
		,SRCDIR.."initfa.c"
		,SRCDIR.."doinit.c"
		,SRCDIR.."karlin.c"
		,SRCDIR.."dropff2.c"
		,SRCDIR.."last_tat.c"
		,SRCDIR.."tatstats.c"
		,SRCDIR.."c_dispn.c"
		,SRCDIR.."lib_sel.c"
		,SRCDIR.."scaleswt.c" 
		,SRCDIR.."mmgetaa.c"
		,SRCDIR.."cal_consf.c"
		,SRCDIR.."mrandom.c"
		,SRCDIR.."url_subs.c"
		,SRCDIR.."pssm_asn_subs.c"
		,SRCDIR.."faatran.c"
,SRCDIR.."comp_lib8.c"
		,SRCDIR.."compacc2.c"
		,SRCDIR.."mshowbest.c"
		,SRCDIR.."build_ares.c"
		,SRCDIR.."mshowalign2.c"
		,SRCDIR.."nmgetlib.c"
		,SRCDIR.."ncbl2_mlib.c"


	}

	defines
	{
		"FASTF"
		,"TFAST"
        ,"COMP_MLIB"
		--,"LALIGN"
		--,"LCAL_CONS"
	}

	configuration "not windows"
		files
		{
			SRCDIR.."mmgetaa.c"
			,SRCDIR.."getseq.c"
		}


	configuration "windows"

		files
		{
			SRCDIR.."getopt.c"
		}


	
project "tfastm36" 
	language    "C"
	kind        "ConsoleApp"

	files
	{

		SRCDIR.."re_getlib.c"
		,SRCDIR.."htime.c"
		,SRCDIR.."apam.c"
		,SRCDIR.."initfa.c"
		,SRCDIR.."doinit.c"
		,SRCDIR.."karlin.c"
		,SRCDIR.."dropfs2.c"
		,SRCDIR.."last_tat.c"
		,SRCDIR.."tatstats.c"
		,SRCDIR.."c_dispn.c"
		,SRCDIR.."lib_sel.c"
		,SRCDIR.."mmgetaa.c"
		,SRCDIR.."cal_consf.c"
		,SRCDIR.."mrandom.c"
		,SRCDIR.."url_subs.c"
		,SRCDIR.."pssm_asn_subs.c"
		,SRCDIR.."faatran.c"
,SRCDIR.."comp_lib8.c"
		,SRCDIR.."compacc2.c"
		,SRCDIR.."mshowbest.c"
		,SRCDIR.."build_ares.c"
		,SRCDIR.."mshowalign2.c"
		,SRCDIR.."nmgetlib.c"
		,SRCDIR.."ncbl2_mlib.c"


	}
	links {"scaleswts" }

	defines
	{
		"FASTM"
		,"TFAST"
        ,"COMP_MLIB"
		--,"LALIGN"
		--,"LCAL_CONS"
	}

	configuration "not windows"
		files
		{
			SRCDIR.."mmgetaa.c"
			,SRCDIR.."getseq.c"
		}


	configuration "windows"

		files
		{
			SRCDIR.."getopt.c"
		}


project "calcons_sw" 
	language    "C"
	kind        "StaticLib"

	files
	{
		SRCDIR.."cal_cons.c"
	}
	
	defines { "SSEARCH" }


project "ssearch36" 
	language    "C"
	kind        "ConsoleApp"

	files
	{

		SRCDIR.."re_getlib.c"
		,SRCDIR.."htime.c"
		,SRCDIR.."apam.c"
		,SRCDIR.."scaleswn.c"
		,SRCDIR.."doinit.c"
		,SRCDIR.."initfa.c"
		,SRCDIR.."scaleswn.c" -- scalese DLOCAL_SCORE --
		,SRCDIR.."dropgsw2.c"
		,SRCDIR.."karlin.c"
		,SRCDIR.."wm_align.c"
		,SRCDIR.."c_dispn.c"
		,SRCDIR.."lib_sel.c"
		,SRCDIR.."mrandom.c"
		,SRCDIR.."url_subs.c"
		,SRCDIR.."pssm_asn_subs.c"
,SRCDIR.."comp_lib8.c"
		,SRCDIR.."compacc2.c"
		,SRCDIR.."mshowbest.c"
		,SRCDIR.."build_ares.c"
		,SRCDIR.."mshowalign2.c"
		,SRCDIR.."nmgetlib.c"
		,SRCDIR.."ncbl2_mlib.c"


	}

	defines
	{
		"SSEARCH"
		,"LOCAL_SCORE"
        ,"COMP_MLIB"
		--,"LALIGN"
		--,"LCAL_CONS"
	}
	
	links { "calcons_sw" }

	configuration "not windows"
		files
		{
			SRCDIR.."mmgetaa.c"
			,SRCDIR.."getseq.c"
		}


	configuration "windows"

		files
		{
			SRCDIR.."getopt.c"
		}

 	 
project "glsearch36" 
	language    "C"
	kind        "ConsoleApp"

	files
	{

		SRCDIR.."re_getlib.c"
		,SRCDIR.."htime.c"
		,SRCDIR.."apam.c"
		,SRCDIR.."scaleswn.c"
		,SRCDIR.."doinit.c"
		,SRCDIR.."initfa.c"
		,SRCDIR.."scaleswn.c" -- scalese DLOCAL_SCORE --
		,SRCDIR.."dropnnw2.c"
		,SRCDIR.."karlin.c"
		,SRCDIR.."wm_align.c"
		,SRCDIR.."c_dispn.c"
		,SRCDIR.."lib_sel.c"
		,SRCDIR.."mrandom.c"
		,SRCDIR.."url_subs.c"
		,SRCDIR.."pssm_asn_subs.c"
,SRCDIR.."comp_lib8.c"
		,SRCDIR.."compacc2.c"
		,SRCDIR.."mshowbest.c"
		,SRCDIR.."build_ares.c"
		,SRCDIR.."mshowalign2.c"
		,SRCDIR.."nmgetlib.c"
		,SRCDIR.."ncbl2_mlib.c"


	}

	links { "calcons_sw" }

	defines
	{
		"GLSEARCH"
		,"NORMAL_DIST"
        ,"COMP_MLIB"
		--,"LALIGN"
		--,"LCAL_CONS"
	}

	configuration "not windows"
		files
		{
			SRCDIR.."mmgetaa.c"
			,SRCDIR.."getseq.c"
		}


	configuration "windows"

		files
		{
			SRCDIR.."getopt.c"
		}

project "lalign36" 
	language    "C"
	kind        "ConsoleApp"

	files
	{

		SRCDIR.."re_getlib.c"
		,SRCDIR.."htime.c"
		,SRCDIR.."apam.c"
		,SRCDIR.."initfa.c"
		,SRCDIR.."doinit.c"
		,SRCDIR.."scaleswn.c" -- scalese DLOCAL_SCORE --
		,SRCDIR.."karlin.c"
		,SRCDIR.."dropgsw2.c"
		,SRCDIR.."wm_align.c"
		,SRCDIR.."cal_cons.c"
		,SRCDIR.."lsim4.c"
		,SRCDIR.."c_dispn.c"
		,SRCDIR.."lib_sel.c"
		,SRCDIR.."scaleswn.c"
		,SRCDIR.."mrandom.c"
		,SRCDIR.."url_subs.c"
		,SRCDIR.."pssm_asn_subs.c"
		,SRCDIR.."last_thresh.c"
,SRCDIR.."comp_lib8.c"
		,SRCDIR.."compacc2.c"
		,SRCDIR.."mshowbest.c"
		,SRCDIR.."build_ares.c"
		,SRCDIR.."mshowalign2.c"
		,SRCDIR.."nmgetlib.c"
		,SRCDIR.."ncbl2_mlib.c"


	}

	defines
	{
		"LALIGN"
		,"LOCAL_SCORE"
		,"LCAL_CONS"
        ,"COMP_MLIB"
	}

	configuration "not windows"
		files
		{
			SRCDIR.."mmgetaa.c"
			,SRCDIR.."getseq.c"
		}


	configuration "windows"

		files
		{
			SRCDIR.."getopt.c"
		}




project "map_db" 
	language    "C"
	kind        "ConsoleApp"

	files
	{
		SRCDIR.."map_db.c"
	}

	
project "lav2ps" 
	language    "C"
	kind        "ConsoleApp"

	files
	{
		SRCDIR.."lav2plt.c"
		,SRCDIR.."lavplt_ps.c"
	}

project "lav2svg" 
	language    "C"
	kind        "ConsoleApp"

	files
	{
		SRCDIR.."lav2plt.c"
		,SRCDIR.."lavplt_svg.c"
	}

	 

 
--lav2ps : lav2plt.o lavplt_ps.o	-DUNIX 
--
--lav2svg : lav2plt.o lavplt_svg.o -DUNIX 
--





--lalign36 : $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(LSHOWALIGN).o htime.o apam.o doinit.o init_lal.o scale_se.o karlin.o last_thresh.o $(DROPLAL_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o pssm_asn_subs.o
--	$(CC) $(HFLAGS) $(BIN)/lalign36 $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(LSHOWALIGN).o htime.o apam.o doinit.o init_lal.o $(DROPLAL_O) scale_se.o karlin.o last_thresh.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o pssm_asn_subs.o $(LIB_M)
--
	
--ggsearch36 : $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o scale_sn.o karlin.o $(DROPGNW_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o pssm_asn_subs.o
--	$(CC) $(HFLAGS) $(BIN)/ggsearch36 $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o $(DROPGNW_O) scale_sn.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o pssm_asn_subs.o $(LIB_M)
--


--glsearch36 : $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o scale_sn.o karlin.o $(DROPLNW_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o pssm_asn_subs.o
--	$(CC) $(HFLAGS) $(BIN)/glsearch36 $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o $(DROPLNW_O) scale_sn.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o pssm_asn_subs.o $(LIB_M)
--

--ssearch36 : $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o scale_se.o karlin.o $(DROPGSW_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o pssm_asn_subs.o
--	$(CC) $(HFLAGS) $(BIN)/ssearch36 $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o $(DROPGSW_O) scale_se.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o pssm_asn_subs.o $(LIB_M)
--







--fasta36 : $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o init_fa.o scale_se.o karlin.o $(DROPNFA_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o mrandom.o url_subs.o
--	$(CC) $(HFLAGS) $(BIN)/fasta36 $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o init_fa.o $(DROPNFA_O) scale_se.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o $(LIB_M)
--

--fasta36_t : $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o init_fa.o scale_se.o karlin.o $(DROPNFA_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o
--	$(CC) $(HFLAGS) $(BIN)/fasta36_t $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o init_fa.o $(DROPNFA_O) scale_se.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o $(LIB_M) $(THR_LIBS)
--
--fasta36sum_t : $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) showsum.o re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o init_fa.o scale_se.o karlin.o $(DROPNFA_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o
--	$(CC) $(HFLAGS) $(BIN)/fasta36sum_t $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) showsum.o re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o init_fa.o $(DROPNFA_O) scale_se.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o $(LIB_M) $(THR_LIBS)
--
--fasta36u_t : $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) showun.o $(SHOWALIGN_T).o htime.o apam.o doinit.o init_fa.o scale_se.o karlin.o $(DROPNFA_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o mrandom.o url_subs.o
--	$(CC) $(HFLAGS) $(BIN)/fasta36u_t $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) showun.o $(SHOWALIGN_T).o htime.o apam.o doinit.o init_fa.o $(DROPNFA_O) scale_se.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o mrandom.o url_subs.o $(LIB_M) $(THR_LIBS)
--
--fasta36r_t : $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) showrel.o $(SHOWALIGN_T).o htime.o apam.o doinit.o init_fa.o scale_se.o karlin.o $(DROPNFA_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o mrandom.o url_subs.o
--	$(CC) $(HFLAGS) $(BIN)/fasta36r_t $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) showrel.o $(SHOWALIGN_T).o htime.o apam.o doinit.o init_fa.o $(DROPNFA_O) scale_se.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o mrandom.o url_subs.o $(LIB_M) $(THR_LIBS)
--
--fastf36_t : $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o init_ff.o scaleswtf.o last_tat.o tatstats_ff.o karlin.o $(DROPFF_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o mrandom.o url_subs.o
--	$(CC) $(HFLAGS) $(BIN)/fastf36_t $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o init_ff.o $(DROPFF_O) scaleswtf.o last_tat.o tatstats_ff.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o mrandom.o url_subs.o $(LIB_M) $(THR_LIBS)
--
--fastf36s_t : $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) showsum.o $(SHOWALIGN_T).o htime.o apam.o doinit.o init_ff.o scaleswtf.o karlin.o $(DROPFF_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o mrandom.o url_subs.o
--	$(CC) $(HFLAGS) $(BIN)/fastf36s_t $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) showsum.o $(SHOWALIGN_T).o htime.o apam.o doinit.o init_ff.o $(DROPFF_O) scaleswtf.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o mrandom.o url_subs.o $(LIB_M) $(THR_LIBS)
--
--fasts36_t : $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o init_fs.o scaleswts.o last_tat.o tatstats_fs.o karlin.o $(DROPFS_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o mrandom.o url_subs.o
--	$(CC) $(HFLAGS) $(BIN)/fasts36_t $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o init_fs.o $(DROPFS_O) scaleswts.o last_tat.o tatstats_fs.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o mrandom.o url_subs.o $(LIB_M) $(THR_LIBS)
--
--fastm36_t : $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o init_fm.o scaleswts.o last_tat.o tatstats_fm.o karlin.o $(DROPFM_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o mrandom.o url_subs.o
--	$(CC) $(HFLAGS) $(BIN)/fastm36_t $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o init_fm.o $(DROPFM_O) scaleswts.o last_tat.o tatstats_fm.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o mrandom.o url_subs.o $(LIB_M) $(THR_LIBS)
--
--fastx36_t : $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o c_dispn.o htime.o apam.o doinit.o init_fx.o faatran.o scale_se.o karlin.o drop_fx.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o
--	$(CC) $(HFLAGS) $(BIN)/fastx36_t $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o init_fx.o drop_fx.o faatran.o scale_se.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o $(LIB_M) $(THR_LIBS)
--
--fasty36_t : $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o c_dispn.o htime.o apam.o doinit.o init_fy.o faatran.o scale_se.o karlin.o drop_fz.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o
--	$(CC) $(HFLAGS) $(BIN)/fasty36_t $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o init_fy.o drop_fz.o faatran.o scale_se.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o $(LIB_M) $(THR_LIBS)
--
--tfasta36 : $(COMP_LIBO) compacc.o $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o init_tfa.o scale_se.o karlin.o $(DROPTFA_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o faatran.o url_subs.o mrandom.o
--	$(CC) $(HFLAGS) $(BIN)/tfasta36 $(COMP_LIBO) compacc.o $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o init_tfa.o $(DROPTFA_O) scale_se.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o faatran.o url_subs.o mrandom.o $(LIB_M)
--
--tfasta36_t : $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o c_dispn.o htime.o apam.o doinit.o init_tfa.o scale_se.o karlin.o $(DROPTFA_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o faatran.o url_subs.o mrandom.o
--	$(CC) $(HFLAGS) $(BIN)/tfasta36_t $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o init_tfa.o $(DROPTFA_O) scale_se.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o faatran.o url_subs.o mrandom.o $(LIB_M) $(THR_LIBS)
--
--tfastf36_t : $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o c_dispn.o htime.o apam.o doinit.o init_tf.o  scaleswtf.o last_tat.o tatstats_ff.o karlin.o $(DROPTFF_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o faatran.o mrandom.o url_subs.o
--	$(CC) $(HFLAGS) $(BIN)/tfastf36_t $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o init_tf.o $(DROPTFF_O) scaleswtf.o last_tat.o tatstats_ff.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o faatran.o mrandom.o url_subs.o $(LIB_M) $(THR_LIBS)
--
--tfasts36_t : $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o c_dispn.o htime.o apam.o doinit.o init_tfs.o scaleswts.o last_tat.o tatstats_fs.o karlin.o $(DROPTFS_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o faatran.o mrandom.o url_subs.o
--	$(CC) $(HFLAGS) $(BIN)/tfasts36_t $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o init_tfs.o $(DROPTFS_O) scaleswts.o last_tat.o tatstats_fs.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o faatran.o mrandom.o url_subs.o $(LIB_M) $(THR_LIBS)
--
--tfastx36_t : $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o init_tfx.o scale_se.o karlin.o drop_tfx.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o faatran.o url_subs.o mrandom.o
--	$(CC) $(HFLAGS) $(BIN)/tfastx36_t $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o init_tfx.o drop_tfx.o scale_se.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o faatran.o url_subs.o mrandom.o $(LIB_M) $(THR_LIBS)
--
--tfasty36_t : $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o init_tfy.o scale_se.o karlin.o drop_tfz.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o faatran.o url_subs.o mrandom.o
--	$(CC) $(HFLAGS) $(BIN)/tfasty36_t $(COMP_THRO) $(WORK_THRO) $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o init_tfy.o drop_tfz.o scale_se.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o faatran.o url_subs.o mrandom.o $(LIB_M) $(THR_LIBS)



--# do not use accelerated Smith-Waterman
--ssearch36s : $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o scale_se.o karlin.o $(DROPGSW_NA_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o pssm_asn_subs.o
--	$(CC) $(HFLAGS) $(BIN)/ssearch36s $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o $(DROPGSW_NA_O) scale_se.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o pssm_asn_subs.o $(LIB_M)
--
--osearch36 : $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o init_ssw.o scale_se.o karlin.o $(DROPNSW_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o
--	$(CC) $(HFLAGS) $(BIN)/osearch36 $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o init_ssw.o $(DROPNSW_O) scale_se.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o $(LIB_M)
--

--prss36 : ssearch36
--	ln -sf ssearch36 prss36
--
--ssearch36_t : $(COMP_THRO) ${WORK_THRO}  $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o scale_se.o karlin.o $(DROPGSW_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o pssm_asn_subs.o
--	$(CC) $(HFLAGS) $(BIN)/ssearch36_t $(COMP_THRO) ${WORK_THRO}  $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o $(DROPGSW_O) scale_se.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o pssm_asn_subs.o $(LIB_M) $(THR_LIBS)
--
--ssearch36s_t : $(COMP_THRO) ${WORK_THRO}  $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o scale_se.o karlin.o $(DROPGSW_NA_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o pssm_asn_subs.o
--	$(CC) $(HFLAGS) $(BIN)/ssearch36s_t $(COMP_THRO) ${WORK_THRO}  $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o $(DROPGSW_NA_O) scale_se.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o pssm_asn_subs.o $(LIB_M) $(THR_LIBS)
--
--glsearch36_t : $(COMP_THRO) ${WORK_THRO}  $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o scale_sn.o karlin.o $(DROPLNW_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o pssm_asn_subs.o
--	$(CC) $(HFLAGS) $(BIN)/glsearch36_t $(COMP_THRO) ${WORK_THRO}  $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o $(DROPLNW_O) scale_sn.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o pssm_asn_subs.o $(LIB_M) $(THR_LIBS)
--
--glsearch36s_t : $(COMP_THRO) ${WORK_THRO}  $(THR_SUBS).o $(COMPACC_TO) showsum.o re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o scale_sn.o karlin.o $(DROPLNW_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o pssm_asn_subs.o
--	$(CC) $(HFLAGS) $(BIN)/glsearch36s_t $(COMP_THRO) ${WORK_THRO}  $(THR_SUBS).o $(COMPACC_TO) showsum.o re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o $(DROPLNW_O) scale_sn.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o pssm_asn_subs.o $(LIB_M) $(THR_LIBS)
--
--ggsearch36_t : $(COMP_THRO) ${WORK_THRO}  $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o scale_sn.o karlin.o $(DROPGNW_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o pssm_asn_subs.o
--	$(CC) $(HFLAGS) $(BIN)/ggsearch36_t $(COMP_THRO) ${WORK_THRO}  $(THR_SUBS).o $(COMPACC_TO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o $(DROPGNW_O) scale_sn.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o pssm_asn_subs.o $(LIB_M) $(THR_LIBS)
--
--ggsearch36s_t : $(COMP_THRO) ${WORK_THRO}  $(THR_SUBS).o $(COMPACC_TO) showsum.o re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o scale_sn.o karlin.o $(DROPGNW_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o pssm_asn_subs.o
--	$(CC) $(HFLAGS) $(BIN)/ggsearch36s_t $(COMP_THRO) ${WORK_THRO}  $(THR_SUBS).o $(COMPACC_TO) showsum.o re_getlib.o $(SHOWALIGN_T).o htime.o apam.o doinit.o $(DROPGNW_O) scale_sn.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o url_subs.o mrandom.o pssm_asn_subs.o $(LIB_M) $(THR_LIBS)
--

--tfastm36 : $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o init_tfm.o scaleswts.o tatstats_fm.o last_tat.o karlin.o $(DROPTFM_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o faatran.o mrandom.o url_subs.o
--	$(CC) $(HFLAGS) $(BIN)/tfastm36 $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o init_tfm.o $(DROPTFM_O) scaleswts.o tatstats_fm.o last_tat.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o faatran.o mrandom.o url_subs.o $(LIB_M)
--



--fastm36 : $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o init_fm.o scaleswts.o last_tat.o tatstats_fm.o karlin.o $(DROPFM_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o mrandom.o url_subs.o
--	$(CC) $(HFLAGS) $(BIN)/fastm36 $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o init_fm.o $(DROPFM_O) scaleswts.o last_tat.o tatstats_fm.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o mrandom.o url_subs.o $(LIB_M)
--


--tfastf36 : $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o init_tf.o scaleswtf.o last_tat.o tatstats_ff.o karlin.o $(DROPTFF_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o faatran.o mrandom.o url_subs.o
--	$(CC) $(HFLAGS) $(BIN)/tfastf36 $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o init_tf.o $(DROPTFF_O) scaleswtf.o last_tat.o tatstats_ff.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o faatran.o mrandom.o url_subs.o $(LIB_M)
--


--fastf36 : $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o init_ff.o scaleswts.o last_tat.o tatstats_ff.o karlin.o $(DROPFF_O) $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o mrandom.o url_subs.o
--	$(CC) $(HFLAGS) $(BIN)/fastf36 $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o init_ff.o $(DROPFF_O) scaleswts.o last_tat.o tatstats_ff.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o mrandom.o url_subs.o $(LIB_M)

--fasty36 : $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o init_fy.o scale_se.o karlin.o drop_fz.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o faatran.o url_subs.o mrandom.o
--	$(CC) $(HFLAGS) $(BIN)/fasty36 $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o init_fy.o drop_fz.o scale_se.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o faatran.o url_subs.o mrandom.o $(LIB_M)
--
--tfasty36 : $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o init_tfy.o scale_se.o karlin.o drop_tfz.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o faatran.o url_subs.o mrandom.o
--	$(CC) $(HFLAGS) $(BIN)/tfasty36 $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o init_tfy.o drop_tfz.o scale_se.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o faatran.o url_subs.o mrandom.o $(LIB_M)
--

	
--tfastx36 : $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o init_tfx.o scale_se.o karlin.o drop_tfx.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o faatran.o url_subs.o mrandom.o
--	$(CC) $(HFLAGS) $(BIN)/tfastx36 $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o init_tfx.o drop_tfx.o scale_se.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o faatran.o url_subs.o mrandom.o $(LIB_M)
--



--fastx36 : $(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o init_fx.o scale_se.o karlin.o drop_fx.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o faatran.o url_subs.o mrandom.o

--	$(COMP_LIBO) $(COMPACC_SO) $(SHOWBESTO) re_getlib.o $(SHOWALIGN_S).o htime.o apam.o doinit.o init_fx.o drop_fx.o scale_se.o karlin.o $(LGETLIB) c_dispn.o $(NCBL_LIB) lib_sel.o faatran.o url_subs.o mrandom.o $(LIB_M)
--


--print_pssm : print_pssm.c getseq.c karlin.c apam.c $(LIB_M)
--
--
--list_db : list_db.c
--



	-- Executables section --
--
--	project "clique"
--		language    "C"
--		kind        "ConsoleApp"
--
--		files
--		{
--			SRCDIR.."clique.c"
--			,INCDIR.."disc.h"
--		}
--
--		
--		links
--
--		{
--			"disc"
--		}

--	project "drawgram"
--		language    "C"
--		kind        "ConsoleApp"
--
--drawgram.o:     drawgram.c draw.h phylip.h
--	$(CC) $(DFLAGS) -c drawgram.c
--
--drawgram:     drawgram.o draw.o draw2.o phylip.o
--	$(CC) $(DFLAGS) draw.o draw2.o drawgram.o phylip.o $(DLIBS) -o drawgram
--
--
--
--	project "drawtree"
--		language    "C"
--		kind        "ConsoleApp"
--
--
--drawtree.o:     drawtree.c draw.h phylip.h
--	$(CC) $(DFLAGS) -c drawtree.c
--
--drawtree:     drawtree.o draw.o draw2.o phylip.o
--	$(CC) $(DFLAGS) draw.o draw2.o drawtree.o phylip.o $(DLIBS) -o drawtree
--







newaction {
   trigger     = "build",
   description = "Run a command line build",
   execute = function ()
	os.chdir(PROJDIR.."build")
	local ver = os.getversion()
	
	print("Running a build")
   	if ver.description == "solaris" then
		print ("Detected platform as solaris")
		os.execute("gmake CC=gcc")
	else
		os.execute("make")
		
	end	
   end
}


newaction {
   trigger     = "clean",
   description = "Run a command line clean",
   execute = function ()
	os.chdir(PROJDIR.."build")
	local ver = os.getversion()
			
	print("Performing a clean") 
   	if ver.description == "solaris" then
		print ("Detected platform as solaris")
		os.execute("gmake clean")
	else
		os.execute("make clean")
	end	


	os.chdir(PROJDIR)
	os.execute("rm -rf build bin lib")

   end
}




--	files 
--	{
--		"*.txt", "**.lua", 
--		"src/**.h", "src/**.c",
--		"src/host/scripts.c"
--	}
--
--	excludes
--	{
--		"src/premake.lua",
--		"src/host/lua-5.1.4/src/lua.c",
--		"src/host/lua-5.1.4/src/luac.c",
--		"src/host/lua-5.1.4/src/print.c",
--		"src/host/lua-5.1.4/**.lua",
--		"src/host/lua-5.1.4/etc/*.c"
--	}
--		
--	configuration "Debug"
--		targetdir   "bin/debug"
--		defines     "_DEBUG"
--		flags       { "Symbols" }
--		
--	configuration "Release"
--		targetdir   "bin/release"
--		defines     "NDEBUG"
--		flags       { "OptimizeSize" }
--
--	configuration "vs*"
--		defines     { "_CRT_SECURE_NO_WARNINGS" }
--	
--	configuration "vs2005"
--		defines	{"_CRT_SECURE_NO_DEPRECATE" }
--
--	configuration "linux"
--		defines     { "LUA_USE_POSIX", "LUA_USE_DLOPEN" }
--		links       { "m", "dl" } 
--		
--	configuration "macosx"
--		defines     { "LUA_USE_MACOSX" }
--		links       { "CoreServices.framework" }
--		
--	configuration { "macosx", "gmake" }
--		buildoptions { "-mmacosx-version-min=10.4" }
--		linkoptions  { "-mmacosx-version-min=10.4" }
--
--	configuration { "linux", "bsd", "macosx" }
--		linkoptions { "-rdynamic" }
--		
--	configuration { "solaris" }
--		linkoptions { "-Wl,--export-dynamic" }
--
--
--
----
---- A more thorough cleanup.
----
--
--if _ACTION == "clean" then
--	os.rmdir("bin")
--	os.rmdir("build")
--end
--


