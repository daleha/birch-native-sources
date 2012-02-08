ROOTDIR=os.getcwd().."/../"


solution "dialign"
	configurations {  "Debug", "Release"}

	location "build"	
	targetdir "bin"


dofile "dialign4.lua"
