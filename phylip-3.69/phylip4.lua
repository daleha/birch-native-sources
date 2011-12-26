

PROJDIR=ROOTDIR
SRCDIR=PROJDIR.."source/"
INCDIR=PROJDIR.."include/"

-- common configuration options
	configuration "Debug"
		buildoptions
		{
			"-O3"
			,"-fomit-frame-pointer" --optimizations

		}

		links
		{
			"m" --math libraries
		}

		includedirs
		{
			INCDIR
		}

		files
		{
			INCDIR.."phylip.h"
		}

		links

		targetdir "bin"



	-- Executables section --

	project "clique "
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
			,"phylip_core"	
		}

	project "consense "
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
			,"phylip_core"
		}


	project "contml "
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
			,"phylip_core"
		}


	project "contrast "
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
			,"phylip_core"
		}	




	project "dnacomp "
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
			,"phylip_core"
		}





	project "dnadist "
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
			,"phylip_core"
		}

	


	project "dnainvar "
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
			,"phylip_core"
		}

	


	project "dnaml "
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
			,"phylip_core"
		}

	


	project "dnamlk "
		language    "C"
		kind        "ConsoleApp"


dnamlk.o:      dnamlk.c seq.h phylip.h mlclock.h printree.h

dnamlk:      dnamlk.o seq.o phylip.o mlclock.o printree.o
	$(CC) $(CFLAGS) dnamlk.o seq.o phylip.o mlclock.o printree.o $(LIBS) -o dnamlk



	project "dnamove "
		language    "C"
		kind        "ConsoleApp"

dnamove.o:      dnamove.c seq.h moves.h phylip.h

dnamove:      dnamove.o seq.o moves.o phylip.o
	$(CC) $(CFLAGS) dnamove.o seq.o moves.o phylip.o $(LIBS) -o dnamove



	project "dnapars "
		language    "C"
		kind        "ConsoleApp"

dnapars.o:      dnapars.c seq.h phylip.h

dnapars:      dnapars.o seq.o phylip.o
	$(CC) $(CFLAGS) dnapars.o seq.o phylip.o $(LIBS) -o dnapars



	project "dnapenny "
		language    "C"
		kind        "ConsoleApp"

dnapenny.o:      dnapenny.c seq.h phylip.h

dnapenny:      dnapenny.o seq.o phylip.o
	$(CC) $(CFLAGS) dnapenny.o seq.o phylip.o $(LIBS) -o dnapenny



	project "dolmove "
		language    "C"
		kind        "ConsoleApp"

dolmove.o:       dolmove.c disc.h moves.h dollo.h phylip.h

dolmove:       dolmove.o disc.o moves.o dollo.o phylip.o
	$(CC) $(CFLAGS) dolmove.o disc.o moves.o dollo.o phylip.o $(LIBS) -o dolmove


	project "dollop "
		language    "C"
		kind        "ConsoleApp"

dollop.o:       dollop.c disc.h dollo.h phylip.h

dollop:       dollop.o disc.o dollo.o phylip.o
	$(CC) $(CFLAGS) dollop.o disc.o dollo.o phylip.o $(LIBS) -o dollop


	project "dolpenny "
		language    "C"
		kind        "ConsoleApp"

dolpenny.o:       dolpenny.c disc.h dollo.h phylip.h

dolpenny:       dolpenny.o disc.o dollo.o phylip.o
	$(CC) $(CFLAGS) dolpenny.o disc.o dollo.o phylip.o $(LIBS) -o dolpenny


	project "factor "
		language    "C"
		kind        "ConsoleApp"

factor.o:       factor.c phylip.h

factor:       factor.o phylip.o
	$(CC) $(CFLAGS) factor.o phylip.o $(LIBS) -o factor


	project "fitch "
		language    "C"
		kind        "ConsoleApp"



fitch.o:        fitch.c dist.h phylip.h

fitch:        fitch.o dist.o phylip.o
	$(CC) $(CFLAGS) fitch.o dist.o phylip.o $(LIBS) -o fitch


	project "gendist "
		language    "C"
		kind        "ConsoleApp"

gendist.o:      gendist.c phylip.h

gendist:      gendist.o phylip.o
	$(CC) $(CFLAGS) gendist.o phylip.o $(LIBS) -o gendist


	project "kitsch "
		language    "C"
		kind        "ConsoleApp"

kitsch.o:        kitsch.c dist.h phylip.h

kitsch:        kitsch.o dist.o phylip.o
	$(CC) $(CFLAGS) kitsch.o dist.o phylip.o $(LIBS) -o kitsch


	project "mix "
		language    "C"
		kind        "ConsoleApp"

mix.o:        mix.c disc.h wagner.h phylip.h

mix:        mix.o disc.o wagner.o phylip.o
	$(CC) $(CFLAGS) mix.o disc.o wagner.o phylip.o $(LIBS) -o mix


	project "move "
		language    "C"
		kind        "ConsoleApp"


move.o:        move.c disc.h moves.h wagner.h phylip.h

move:        move.o disc.o moves.o wagner.o phylip.o
	$(CC) $(CFLAGS) move.o disc.o moves.o wagner.o phylip.o $(LIBS) -o move


	project "neighbor "
		language    "C"
		kind        "ConsoleApp"


neighbor.o:        neighbor.c dist.h phylip.h

neighbor:        neighbor.o dist.o phylip.o
	$(CC) $(CFLAGS) neighbor.o dist.o phylip.o $(LIBS) -o neighbor


	project "pars "
		language    "C"
		kind        "ConsoleApp"

pars.o:   pars.c discrete.h phylip.h

pars: pars.o discrete.o phylip.o
	$(CC) $(CFLAGS) pars.o discrete.o phylip.o $(LIBS) -o pars


	project "penny "
		language    "C"
		kind        "ConsoleApp"

penny.o:  penny.c disc.h wagner.h phylip.h

penny:  penny.o disc.o wagner.o phylip.o
	$(CC) $(CFLAGS) penny.o disc.o wagner.o  phylip.o $(LIBS) -o penny


	project "proml "
		language    "C"
		kind        "ConsoleApp"

proml.o:      proml.c seq.h phylip.h

proml:      proml.o seq.o phylip.o
	$(CC) $(CFLAGS) proml.o seq.o phylip.o $(LIBS) -o proml


	project "promlk "
		language    "C"
		kind        "ConsoleApp"


promlk.o:      promlk.c seq.h phylip.h mlclock.h printree.h

promlk:      promlk.o seq.o phylip.o mlclock.o printree.o
	$(CC) $(CFLAGS) promlk.o seq.o phylip.o mlclock.o printree.o $(LIBS) -o promlk


	project "protdist "
		language    "C"
		kind        "ConsoleApp"

protdist.o:      protdist.c seq.h phylip.h

protdist:      protdist.o seq.o phylip.o
	$(CC) $(CFLAGS) protdist.o seq.o phylip.o $(LIBS) -o protdist


	project "protpars "
		language    "C"
		kind        "ConsoleApp"

protpars.o: protpars.c seq.h phylip.h

protpars: protpars.o seq.o phylip.o
	$(CC) $(CFLAGS) protpars.o seq.o phylip.o $(LIBS) -o protpars


	project "restdist "
		language    "C"
		kind        "ConsoleApp"

restdist.o: restdist.c seq.h phylip.h

restdist: restdist.o seq.o phylip.o
	$(CC) $(CFLAGS) restdist.o seq.o phylip.o $(LIBS) -o restdist


	project "restml "
		language    "C"
		kind        "ConsoleApp"

restml.o: restml.c seq.h phylip.h

restml: restml.o seq.o phylip.o
	$(CC) $(CFLAGS) restml.o seq.o phylip.o $(LIBS) -o restml


	project "retree "
		language    "C"
		kind        "ConsoleApp"

retree.o:       retree.c moves.h phylip.h

retree:       retree.o moves.o phylip.o
	$(CC) $(CFLAGS) retree.o moves.o phylip.o $(LIBS) -o retree


	project "seqboot "
		language    "C"
		kind        "ConsoleApp"


seqboot.o:      seqboot.c phylip.h

seqboot:      seqboot.o seq.o phylip.o
	$(CC) $(CFLAGS) seqboot.o seq.o phylip.o $(LIBS) -o seqboot



	project "treedist "
		language    "C"
		kind        "ConsoleApp"


treedist.o:     treedist.c cons.h phylip.h

treedist:     treedist.o phylip.o cons.o
	$(CC) $(CFLAGS) treedist.o cons.o phylip.o $(LIBS) -o treedist



	project "drawgram "
		language    "C"
		kind        "ConsoleApp"

drawgram.o:     drawgram.c draw.h phylip.h
	$(CC) $(DFLAGS) -c drawgram.c

drawgram:     drawgram.o draw.o draw2.o phylip.o
	$(CC) $(DFLAGS) draw.o draw2.o drawgram.o phylip.o $(DLIBS) -o drawgram



	project "drawtree"
		language    "C"
		kind        "ConsoleApp"


drawtree.o:     drawtree.c draw.h phylip.h
	$(CC) $(DFLAGS) -c drawtree.c

drawtree:     drawtree.o draw.o draw2.o phylip.o
	$(CC) $(DFLAGS) draw.o draw2.o drawtree.o phylip.o $(DLIBS) -o drawtree














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


