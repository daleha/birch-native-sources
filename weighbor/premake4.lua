ROOTDIR=os.getcwd().."/../"


solution "weighbor"
	configurations {  "Debug", "Release"}

	location "build"	
	targetdir "bin"


dofile "weighbor4.lua"
