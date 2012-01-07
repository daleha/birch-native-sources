
PROJDIR=ROOTDIR
SRCDIR=PROJDIR.."source/"
INCDIR=PROJDIR.."include/"

-- common configuration options

	-- shared libs section --

	project "phylip_core"     
		language    "C"
		kind        "StaticLib"
		targetdir "lib"

		files 
		{ 
			SRCDIR.."phylip.c" 
		}
		

	project "seq"        
		language    "C"
		kind        "StaticLib"
		targetdir "lib"

		files 
		{
			SRCDIR.."seq.c"
			,INCDIR.."seq.h"
		}

	project "disc"       
		language    "C"
		kind        "StaticLib"
		targetdir "lib"

		files
		{
			SRCDIR.."disc.c"
			,INCDIR.."disc.h"
		}

	project "discrete"   
		language    "C"
		kind        "StaticLib"
		targetdir "lib"

		files
		{
			SRCDIR.."discrete.c"
			,INCDIR.."discrete.h"
		}

	project "dollo"      
		language    "C"
		kind        "StaticLib"
		targetdir "lib"

		files
		{
			SRCDIR.."dollo.c"
			,INCDIR.."dollo.h"

		}

	project "wagner"     
		language    "C"
		kind        "StaticLib"
		targetdir "lib"

		files
		{
			SRCDIR.."wagner.c"
			,INCDIR.."wagner.h"
		}

	project "dist"       
		language    "C"
		kind        "StaticLib"
		targetdir "lib"

		files
		{
			SRCDIR.."dist.c"
			,INCDIR.."dist.h"
		}

	project "cont"       
		language    "C"
		kind        "StaticLib"
		targetdir "lib"

		files
		{
			SRCDIR.."cont.c"
			,INCDIR.."cont.h"
		}

	project "mlclock"    
		language    "C"
		kind        "StaticLib"
		targetdir "lib"

		files
		{
			SRCDIR.."mlclock.c"
			,INCDIR.."mlclock.h"
		}
	
	project "moves"      
		language    "C"
		kind        "StaticLib"
		targetdir "lib"

		files
		{
			SRCDIR.."moves.c"
			,INCDIR.."moves.h"
		}

	project "printree"   
		language    "C"
		kind        "StaticLib"
		targetdir "lib"

		files
		{
			SRCDIR.."printree.c"
			,INCDIR.."printree.h"
		}

	project "cons"
		language    "C"
		kind        "StaticLib"
		targetdir "lib"


		files
		{

			SRCDIR.."cons.c" 
			,INCDIR.."cons.h"
		}		

	project "draw"
		language    "C"
		kind        "StaticLib"
		targetdir "lib"

		files
		{

			SRCDIR.."draw.c"
			,INCDIR.."draw.h"
		}

	-- Graphics libs --


	--convenience lib
	project "phylip_gui"
		language    "C"
		kind        "StaticLib"
		targetdir "lib"

		files
		{

			SRCDIR.."draw.c"
			,SRCDIR.."draw2.c"
			,INCDIR.."draw.h"
			,INCDIR.."draw2.h"
		}




	project "draw"
		language    "C"
		kind        "StaticLib"
		targetdir "lib"

		files
		{

			SRCDIR.."draw.c"
			,INCDIR.."draw.h"
		}



	project "draw2"

		language    "C"
		kind        "StaticLib"
		targetdir "lib"

		files
		{

			SRCDIR.."draw.c"
			,INCDIR.."draw.h"
		}


