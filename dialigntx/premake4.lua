ROOTDIR=os.getcwd().."/../"


solution "dialigntx"
	configurations {  "Debug", "Release"}

	location "build"	
	targetdir "bin"


dofile "dialigntx4.lua"
