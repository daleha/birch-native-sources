

PROJDIR=ROOTDIR.."dialigntx/"
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
			
		links 
		{ 
			"m"
		}	
	
	configuration "windows"
		
		defines
		{
			"__MINGW__"
		}	
	

	project "dialigntx" 
		language    "C"
		kind        "ConsoleApp"

		files
		{

			SRCDIR.."*.c" 

		}

