

PROJDIR=ROOTDIR
SRCDIR=PROJDIR.."source/"
INCDIR=PROJDIR.."include/"


-- common configuration options
	configuration "Debug"
		buildoptions
		{
			"-O3"

		}


		includedirs
		{
			INCDIR --add the header folder to include search path
		}

		files
		{
			INCDIR.."phylip.h"
		}

		links
		{
			"m" --math libraries
			,"phylip_core"  -- everything links phylip core
		}

		location "build"
		targetdir "bin" -- and the bins into bin

		dofile "./phylip4_lib.lua"

	configuration "linux"

		buildoptions
		{
			"-fomit-frame-pointer" --optimizations
		}


	configuration "solaris"

		buildoptions
		{
			"-fomit-frame-pointer" --optimizations
		}


	configuration "macosx"

		buildoptions
		{
			"-Wall -mmacosx-version-min=10.1"	
		}
		linkoptions
		{
			"-Wl"	
		}



	-- Executables section --

	project "clique"
		language    "C"
		kind        "ConsoleApp"

		files
		{
			SRCDIR.."clique.c"
			,INCDIR.."disc.h"
		}

		
		links

		{
			"disc"
		}

	project "consense"
		language    "C"
		kind        "ConsoleApp"

		files
		{
			SRCDIR.."consense.c"
			,INCDIR.."cons.h"
		}


		links
		
		{
			"cons"
		}


	project "contml"
		language    "C"
		kind        "ConsoleApp"

		files
		{
			SRCDIR.."contml.c"
			,INCDIR.."cont.h"
		}

		links
		{
			"cont"
		}


	project "contrast"
		language    "C"
		kind        "ConsoleApp"

		files
		{

			SRCDIR.."contrast.c"
			,INCDIR.."cont.h"
		}


		links
		{
			"cont"
		}	




	project "dnacomp"
		language    "C"
		kind        "ConsoleApp"

		files
		{
			SRCDIR.."dnacomp.c"
			,INCDIR.."seq.h"
		}

		links
		{
			"seq"
		}





	project "dnadist"
		language    "C"
		kind        "ConsoleApp"


		files
		{
			SRCDIR.."dnadist.c"
			,INCDIR.."seq.h"
		}

		links
		{
			"seq"
		}

	


	project "dnainvar"
		language    "C"
		kind        "ConsoleApp"



		files
		{
			SRCDIR.."dnainvar.c"
			,INCDIR.."seq.h"
		}

		links
		{
			"seq"
		}

	


	project "dnaml"
		language    "C"
		kind        "ConsoleApp"

		files
		{
			SRCDIR.."dnaml.c"
			,INCDIR.."seq.h"


		}

		links
		{
			"seq"
		}

	


	project "dnamlk"
		language    "C"
		kind        "ConsoleApp"

		files
		{
			SRCDIR.."dnamlk.c"
			,INCDIR.."seq.h" 
			,INCDIR.."mlclock.h" 
			,INCDIR.."printree.h"
		}
		
		links
		{
			"seq"
			,"mlclock"
			,"printree"
		}	



	project "dnamove"
		language    "C"
		kind        "ConsoleApp"

		files
		{
			SRCDIR.."dnamove.c"
			,INCDIR.." seq.h" 
			,INCDIR.."moves.h" 
		}

		links
		{
			"seq"
			,"moves"
		}


	project "dnapars"
		language    "C"
		kind        "ConsoleApp"


		files
		{
			SRCDIR.."dnapars.c"
			,INCDIR.."seq.h"
		}
	
		links
		{
			"seq"
		}



	project "dnapenny"
		language    "C"
		kind        "ConsoleApp"


		files
		{
			SRCDIR.."dnapenny.c"
			,INCDIR.."seq.h"
		}

		links
		{
			"seq"
		}



	project "dolmove"
		language    "C"
		kind        "ConsoleApp"

		files
		{
			SRCDIR.."dolmove.c"
			,INCDIR.."disc.h"
			,INCDIR.."moves.h"
			,INCDIR.."dollo.h"
		}

		
		links
		{
			"disc"
			,"moves"
			,"dollo"
		}


	project "dollop"
		language    "C"
		kind        "ConsoleApp"


		files
		{
			SRCDIR.."dollop.c"
			,INCDIR.."disc.h"
			,INCDIR.."dollo.h"
			
		}
		
		links
		{
			"disc"
			,"dollo"	
		}


	project "dolpenny"
		language    "C"
		kind        "ConsoleApp"

		files
		{
			SRCDIR.."dolpenny.c"
			,INCDIR.."disc.h"
			,INCDIR.."dollo.h"
		}

		links
		{
			"disc"
			,"dollo"
		}


	project "factor"
		language    "C"
		kind        "ConsoleApp"

		files
		{
			SRCDIR.."factor.c"
		}


	project "fitch"
		language    "C"
		kind        "ConsoleApp"


		files
		{	
			SRCDIR.."fitch.c"
			,INCDIR.."dist.h"
		}

		links
		{
			"dist"
		}


	project "gendist"
		language    "C"
		kind        "ConsoleApp"

		files
		{
			SRCDIR.."gendist.c"
		}
		

	project "kitsch"
		language    "C"
		kind        "ConsoleApp"

		files
		{
			SRCDIR.."kitsch.c"
			,INCDIR.."dist.h"
		}

		links
		{
			"dist"
		}


	project "mix"
		language    "C"
		kind        "ConsoleApp"

		files
		{
			SRCDIR.."mix.c"
			,INCDIR.."disc.h"
			,INCDIR.."wagner.h"
		}
	
		links
		{
			"disc"
			,"wagner"
		}


	project "move"
		language    "C"
		kind        "ConsoleApp"

		files
		{
			SRCDIR.."move.c"
			,INCDIR.."disc.h"
			,INCDIR.."moves.h"
			,INCDIR.."wagner.h"
		}
		
		links
		{
			"disc"
			,"moves"
			,"wagner"
		}


	project "neighbor"
		language    "C"
		kind        "ConsoleApp"

		files
		{
			SRCDIR.."neighbor.c"
			,INCDIR.."dist.h"
			
		}

		links
		{
			"dist"
		}


	project "pars"
		language    "C"
		kind        "ConsoleApp"

		files
		{
			SRCDIR.."pars.c"
			,INCDIR.."discrete.h"
		}

		links
		{
			"discrete"
		}

	project "penny"
		language    "C"
		kind        "ConsoleApp"

		files
		{
			SRCDIR.."penny.c"
			,INCDIR.."disc.h"
			,INCDIR.."wagner.h"
			
		}

		links
		{
			"disc"
			,"wagner"
		}


	project "proml"
		language    "C"
		kind        "ConsoleApp"

		files
		{
			SRCDIR.."proml.c"
			,INCDIR.."seq.h"
		}

		links
		{
			"seq"
		}



	project "promlk"
		language    "C"
		kind        "ConsoleApp"


		files
		{
			SRCDIR.."promlk.c"
			,INCDIR.."seq.h"
			,INCDIR.."mlclock.h"
			,INCDIR.."printree.h"
		}

		links
		{
			"seq"
			,"mlclock"
			,"printree"
		}


	project "protdist"
		language    "C"
		kind        "ConsoleApp"

		
		files
		{
			SRCDIR.."protdist.c"
			,INCDIR.."seq.h"

		}


		configuration "macosx"
			files { SRCDIR.."seq.c" }
		

		configuration "not macosx"
			links
			{
				"seq"
			}




	project "protpars"
		language    "C"
		kind        "ConsoleApp"

		files
		{
			SRCDIR.."protpars.c"
			,INCDIR.."seq.h"
		}
		
		links
		{
			"seq"
		}


	project "restdist"
		language    "C"
		kind        "ConsoleApp"

		files
		{
			SRCDIR.."restdist.c"
			,INCDIR.."seq.h"
		}

		links
		{
			"seq"
		}


	project "restml"
		language    "C"
		kind        "ConsoleApp"

		files
		{
			SRCDIR.."restml.c"
			,INCDIR.."seq.h"	
		}

		links
		{
			"seq"
		}


	project "retree"
		language    "C"
		kind        "ConsoleApp"

		files
		{
			SRCDIR.."retree.c"
			,INCDIR.."moves.h"
		}	

		links
		{
			"moves"
		}


	project "seqboot"
		language    "C"
		kind        "ConsoleApp"

		files
		{
			SRCDIR.."seqboot.c"
		}

		configuration "macosx"
			files { SRCDIR.."seq.c" }
		

		configuration "not macosx"
			links
			{
				"seq"
			}



		links
		{
			"seq"
		}


	project "treedist"
		language    "C"
		kind        "ConsoleApp"

		files
		{
			SRCDIR.."treedist.c"
			,INCDIR.."cons.h"
		}
	
		links
		{
			"cons"
		}



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


