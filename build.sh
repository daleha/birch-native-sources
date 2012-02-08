#!/bin/bash


# usage: build.sh [premake path] [ make path ] [cc type ]
rm -rf bin
mkdir bin

PREMAKE=$1
MAKE=$2 # for solaris, set this to gmake
MAKEOPTS="CC=$3" # for solaris, set this to CC=gcc
ROOTDIR=`pwd`

cd phylip
"$PREMAKE" gmake
cd build
"$MAKE" clean
"$MAKE" "$MAKEOPTS"
cd $ROOTDIR

cd fasta
"$PREMAKE" gmake
cd build
"$MAKE" clean
"$MAKE" "$MAKEOPTS"
cd $ROOTDIR

cd xylem
"$PREMAKE" gmake
cd build
"$MAKE" clean
"$MAKE" "$MAKEOPTS"
cd $ROOTDIR

cd mrbayes
"$PREMAKE" gmake
cd build
"$MAKE" clean
"$MAKE" "$MAKEOPTS"
cd $ROOTDIR



echo "Copying binaries..."
cp */bin/* bin
rm bin/*.a # delete un-needed static libs
#rm -rf */bin */build #delete intermediate build folders
