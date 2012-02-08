
ROOTDIR=os.getcwd().."/"


solution "FASTA"
	configurations {  "Debug", "Release"}

	location "build"
	targetdir "bin" -- and the bins into bin


dofile "fasta4.lua"
