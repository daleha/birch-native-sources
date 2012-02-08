ROOTDIR=os.getcwd().."/../"


solution "xylem"
	configurations {  "Debug", "Release"}

	location "build"	
	targetdir "bin"


dofile "xylem4.lua"
