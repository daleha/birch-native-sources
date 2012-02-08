

PROJDIR=ROOTDIR.."dialign/"
SRCDIR=PROJDIR.."source/"
INCDIR=PROJDIR.."include/"



-- common configuration options
	configuration "Debug"


		includedirs
		{
			INCDIR --add the header folder to include search path
		}
		
		defines
		{
			"CONS"
		}
			
		links 
		{ 
			"m"
		}	
		
		buildoptions
		{
			"-O"
		}		
	

	project "dialign" 
		language    "C"
		kind        "ConsoleApp"

		files
		{

			SRCDIR.."*.c" 

		}

