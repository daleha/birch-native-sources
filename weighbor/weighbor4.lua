

PROJDIR=ROOTDIR.."weighbor/"
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
		--	"USECONFIG_H"
		}
		
		buildoptions
		{
			"-O4"
		}	
		links 
		{ 
			"m"
		}	
		
	

	project "weighbor" 
		language    "C"
		kind        "ConsoleApp"

		files
		{

			SRCDIR.."*.c" 

		}

