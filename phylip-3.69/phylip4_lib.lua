
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

		targetdir "lib"

	-- shared libs section --

	project "phylip_core"     
		language    "C"
		kind        "SharedLib"

		files 
		{ 
			SRCDIR.."phylip.c" 
		}

	project "seq"        
		language    "C"
		kind        "SharedLib"

		files 
		{
			SRCDIR.."seq.c"
			,INCDIR.."seq.h"
		}

	project "disc"       
		language    "C"
		kind        "SharedLib"

		files
		{
			SRCDIR.."disc.c"
			,INCDIR.."disc.h"
		}

	project "discrete"   
		language    "C"
		kind        "SharedLib"

		files
		{
			SRCDIR.."discrete.c"
			,INCDIR.."discrete.h"
		}

	project "dollo"      
		language    "C"
		kind        "SharedLib"

		files
		{
			SRCDIR.."dollo.c"
			,INCDIR.."dollo.h"

		}

	project "wagner"     
		language    "C"
		kind        "SharedLib"

		files
		{
			SRCDIR.."wagner.c"
			INCDIR.."wagner.h"
		}

	project "dist"       
		language    "C"
		kind        "SharedLib"

		files
		{
			SRCDIR.."dist.c"
			,INCDIR.."dist.h"
		}

	project "cont"       
		language    "C"
		kind        "SharedLib"

		files
		{
			SRCDIR.."cont.c"
			,INCDIR.."cont.h"
		}

	project "mlclock"    
		language    "C"
		kind        "SharedLib"

		files
		{
			SRCDIR.."mlclock.c"
			,INCDIR.."mlclock.h"
		}
	
	project "moves"      
		language    "C"
		kind        "SharedLib"

		files
		{
			SRCDIR.."moves.c"
			,INCDIR.."moves.h"
		}

	project "printree"   
		language    "C"
		kind        "SharedLib"

		files
		{
			SRCDIR.."printree.c"
			INCDIR.."printree.h"
		}

	project "cons"
		language    "C"
		kind        "SharedLib"


		files
		{

			SRCDIR.."cons.c" 
			,INCDIR.."cons.h"
		}		

	project "draw"
		language    "C"
		kind        "SharedLib"

		files
		{

			SRCDIR.."draw.c"
			,INCDIR.."draw.h"
		}

	-- Graphics libs --


	--convenience lib
	project "phylip_gui"
		language    "C"
		kind        "SharedLib"

		files
		{

			SRCDIR.."draw.c"
			,SRCDIR.."draw2.c"
			,INCDIR.."draw.h"
			,INCDIR.."draw2.h"
		}




	project "draw"
		language    "C"
		kind        "SharedLib"

		files
		{

			SRCDIR.."draw.c"
			,INCDIR.."draw.h"
		}



	project "draw2"

		language    "C"
		kind        "SharedLib"

		files
		{

			SRCDIR.."draw.c"
			,INCDIR.."draw.h"
		}


