

PROJDIR=ROOTDIR.."fsap/"
SRCDIR=PROJDIR.."source/"
INCDIR=PROJDIR.."include/"



-- common configuration options
	configuration "Debug"


		includedirs
		{
			INCDIR --add the header folder to include search path
		}

			
		links 
		{ 
			"m"
		}	
		
		location "build"	
		targetdir "bin"


	project "p2c" 
		language    "C"
		kind        "StaticLib"

		files
		{
			SRCDIR.."p2clib.c" 
		}



	project "bachrest" 
		language    "C"
		kind        "ConsoleApp"

		links { "p2c" }
		files
		{

			SRCDIR.."bachrest.c" 

		}

	project "d3hom" 
		language    "C"
		kind        "ConsoleApp"

		links { "p2c" }
		files
		{

			SRCDIR.."d3hom.c" 

		}


	project "d4hom" 
		language    "C"
		kind        "ConsoleApp"

		links { "p2c" }
		files
		{

			SRCDIR.."d4hom.c" 

		}


	project "multidigest"
		language    "C"
		kind        "ConsoleApp"

		links { "p2c" }
		files
		{

			SRCDIR.."multidigest.c" 

		}


	project "funnel"
		language    "C"
		kind        "ConsoleApp"

		links { "p2c" }
		files
		{

			SRCDIR.."funnel.c" 

		}


	project "gel"
		language    "C"
		kind        "ConsoleApp"

		links { "p2c" }
		files
		{

			SRCDIR.."gel.c" 

		}


	project "intrest"
		language    "C"
		kind        "ConsoleApp"

		links { "p2c" }
		files
		{

			SRCDIR.."intrest.c" 

		}


	project "numseq"
		language    "C"
		kind        "ConsoleApp"

		links { "p2c" }
		files
		{

			SRCDIR.."numseq.c" 

		}

	project "p1hom"
		language    "C"
		kind        "ConsoleApp"

		links { "p2c" }
		files
		{

			SRCDIR.."p1hom.c" 

		}


	project "p2hom"
		language    "C"
		kind        "ConsoleApp"

		links { "p2c" }
		files
		{

			SRCDIR.."p2hom.c" 

		}


	project "prostat"
		language    "C"
		kind        "ConsoleApp"

		links { "p2c" }
		files
		{

			SRCDIR.."prostat.c" 

		}


	project "testcode"
		language    "C"
		kind        "ConsoleApp"

		links { "p2c" }
		files
		{

			SRCDIR.."testcode.c" 

		}






