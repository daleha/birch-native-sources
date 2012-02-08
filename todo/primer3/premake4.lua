ROOTDIR=os.getcwd().."/../"


solution "primer3"
	configurations {  "Debug", "Release"}

	location "build"	
	targetdir "bin"


dofile "primer34.lua"
