ROOTDIR=os.getcwd().."/../"


solution "fsap"
	configurations {  "Debug", "Release"}

	location "build"	
	targetdir "bin"


dofile "fsap4.lua"
