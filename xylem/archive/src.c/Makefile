#This Makefile works with Sun cc and gcc

#sun (bsd)
BIN = ../bin
INCLIB = .
LDLIB = .
LIBOBJS = p2clib.o
#CC= gcc
CC= cc
CFLAGS= -O  -L$(LDLIB) -lp2c -lm -I$(INCLIB)
LFLAGS=  -c

PROGS= dbstat getloc getob identify prot2nuc reform ribosome splitdb shuffle flat2phyl phyl2flat

all : libp2c.a $(PROGS)

clean : 
	rm libp2c.a *.o $(PROGS)

install : 
	mv $(PROGS) $(BIN)

# Making the p2c runtime library
libp2c.a: $(LIBOBJS)
	ar r libp2c.a $(LIBOBJS); ranlib libp2c.a

p2clib.o: p2clib.c
	$(CC) $(LFLAGS)  p2clib.c

#libp2c.a: libp2c.a
#	ranlib libp2c.a

# Making the programs
dbstat :   
	$(CC) dbstat.c $(CFLAGS)  -o dbstat

flat2phyl :   
	$(CC) flat2phyl.c $(CFLAGS)  -o flat2phyl

phyl2flat :   
	$(CC) phyl2flat.c $(CFLAGS)  -o phyl2flat

getloc :   
	$(CC) getloc.c $(CFLAGS) -o getloc

getob :   
	$(CC) getob.c $(CFLAGS)  -o getob

identify :   
	$(CC) identify.c $(CFLAGS) -o identify 

prot2nuc :   
	$(CC) prot2nuc.c $(CFLAGS) -o prot2nuc 

reform :   
	$(CC) reform.c $(CFLAGS) -o reform 

ribosome :   
	$(CC) ribosome.c $(CFLAGS) -o ribosome 

shuffle :   
	$(CC) shuffle.c $(CFLAGS) -o shuffle 

splitdb :   
	$(CC) splitdb.c $(CFLAGS) -o splitdb 

