

PROJDIR=ROOTDIR.."xylem/"
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



	project "dbstat" 
		language    "C"
		kind        "ConsoleApp"

		links { "p2c" }
		files
		{

			SRCDIR.."dbstat.c" 

		}



	project "getloc" 
		language    "C"
		kind        "ConsoleApp"

		links { "p2c" }
		files
		{

			SRCDIR.."getloc.c" 

		}


	project "getob" 
		language    "C"
		kind        "ConsoleApp"

		links { "p2c" }
		files
		{

			SRCDIR.."getob.c" 

		}


	project "identify"
		language    "C"
		kind        "ConsoleApp"

		links { "p2c" }
		files
		{

			SRCDIR.."identify.c" 

		}


	project "prot2nuc"
		language    "C"
		kind        "ConsoleApp"

		links { "p2c" }
		files
		{

			SRCDIR.."prot2nuc.c" 

		}


	project "reform"
		language    "C"
		kind        "ConsoleApp"

		links { "p2c" }
		files
		{

			SRCDIR.."reform.c" 

		}


	project "ribosome"
		language    "C"
		kind        "ConsoleApp"

		links { "p2c" }
		files
		{

			SRCDIR.."ribosome.c" 

		}


	project "splitdb"
		language    "C"
		kind        "ConsoleApp"

		links { "p2c" }
		files
		{

			SRCDIR.."splitdb.c" 

		}


	project "shuffle"
		language    "C"
		kind        "ConsoleApp"

		links { "p2c" }
		files
		{

			SRCDIR.."shuffle.c" 

		}


	project "flat2phyl"
		language    "C"
		kind        "ConsoleApp"

		links { "p2c" }
		files
		{

			SRCDIR.."flat2phyl.c" 

		}


	project "phyl2flat"
		language    "C"
		kind        "ConsoleApp"

		links { "p2c" }
		files
		{

			SRCDIR.."phyl2flat.c" 

		}







