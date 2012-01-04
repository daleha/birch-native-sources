
PROJDIR=ROOTDIR

SRCDIR=PROJDIR.."lib-src/"
GUIDIR=PROJDIR.."gui-src/"

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


		targetdir "lib"

	-- shared libs section --

	project "phylip_core"     
		language    "C"
		kind        "SharedLib"

		files 
		{ 
			SRCDIR.."**.c"
			,INCDIR.."*.h"
		}

		configuration "not macosx"
			excludes
			{
				SRCDIR.."interface.c"
				,SRCDIR.."macface.c"
				,INCDIR.."interface.h"
			}

		excludes
		{
			SRCDIR.."newmove.c"
		}



	-- Gui lib
	project "phylip_gui"
		language    "C"
		kind        "SharedLib"

		files
		{

			PROJDIR.."gui/**.c"
			,INCDIR.."/gui/**.h"
		}

