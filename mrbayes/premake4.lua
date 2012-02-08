ROOTDIR=os.getcwd().."/../"


solution "mrbayes"
	configurations {  "Debug", "Release"}

	location "build"	
	targetdir "bin"


dofile "mrbayes4.lua"
