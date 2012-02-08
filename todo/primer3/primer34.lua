

PROJDIR=ROOTDIR.."primer3/"
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
			"CONS"
		}
			
		links 
		{ 
			"m"
		}	
		
		buildoptions
		{
			"-O"
		}		
	

	project "oligotmlib" 
		language    "C++"
		kind        "StaticLib"

		files
		{
			SRCDIR.."oligotm.c" 
		}

	project "dpal" 
		language    "C++"
		kind        "StaticLib"

		files
		{
			SRCDIR.."dpal.c" 
		}

	project "thal" 
		language    "C++"
		kind        "StaticLib"

		buildoptions
		{
			"-ffloat-store"
		}
		files
		{
			SRCDIR.."thal.c" 
		}

	project "primer3" 
		language    "C++"
		kind        "StaticLib"

		files
		{
			SRCDIR.."p3_seq_lib.c" 
			,SRCDIR.."libprimer3.c" 
		}



--LIBOLIGOTM      = liboligotm.a
--oligotm.o: oligotm.c 
--
--LIBDPAL         = libdpal.a
--dpal_primer.o: dpal.c 
--
--LIBTHAL         = libthal.a
--thal_primer.o: thal.c -ffloat-store 
--
--
--LIBPRIMER3      = libprimer3.a
--p3_seq_lib.o: p3_seq_lib.c 
--libprimer3.o: libprimer3.c -Wno-deprecated 
--
--# We use '-ffloat-store' on windows to prevent undesirable
--# precision which may lead to differences in floating point results.
--
--
--
--NTDPAL_EXE	= ntdpal
--dpal.o: dpal.c 
--ntdpal_main.o: ntdpal_main.c
--
--NTTHAL_EXE	= ntthal
--thal.o: thal.c -ffloat-store 
--thal_main.o: thal_main.c #no  -O2
--
--
--PRIMER_EXE      = primer3_core
--primer3_boulder_main.o: primer3_boulder_main.c 
--read_boulder.o: read_boulder.c 
--print_boulder.o: print_boulder.c 
--format_output.o: format_output.c 
--
--
--
--OLIGOTM_EXE	= oligotm
--oligotm_main.c $(LIBOLIGOTM)
--
--LONG_SEQ_EXE	= long_seq_tm_test
--long_seq_tm_test_main.c oligotm.o



