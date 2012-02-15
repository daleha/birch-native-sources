

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

	


--newaction {
--   trigger     = "build",
--   description = "Run a command line build",
--   execute = function ()
--	os.chdir(PROJDIR.."build")
--	local ver = os.getversion()
--	
--	print("Running a build")
--   	if ver.description == "solaris" then
--		print ("Detected platform as solaris")
--		os.execute("gmake CC=gcc")
--	else
--		os.execute("make")
--		
--	end	
--   end
--}
--
--
--newaction {
--   trigger     = "clean",
--   description = "Run a command line clean",
--   execute = function ()
--	os.chdir(PROJDIR.."build")
--	local ver = os.getversion()
--			
--	print("Performing a clean") 
--   	if ver.description == "solaris" then
--		print ("Detected platform as solaris")
--		os.execute("gmake clean")
--	else
--		os.execute("make clean")
--	end	
--
--
--	os.chdir(PROJDIR)
--	os.execute("rm -rf build bin lib")
--
--   end
--}
--



