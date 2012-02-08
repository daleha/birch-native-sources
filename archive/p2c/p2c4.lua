

PROJDIR=ROOTDIR.."lib/p2c/"
SRCDIR=PROJDIR.."source/"
INCDIR=PROJDIR.."include/"

	project "p2c" 
		language    "C"
		kind        "StaticLib"

		includedirs
		{
			INCDIR
		}
	
		files
		{
			SRCDIR.."p2clib.c" 
			,SRCDIR.."loc.p2clib.c"
		}


