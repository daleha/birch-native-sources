#!/bin/bash

rm -rf bin
mkdir bin

MAKE="make" # for solaris, set this to gmake
MAKEOPTS="" # for solaris, set this to CC=gcc
PREMAKE=`which premake4`
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



echo "Copying binaries..."
cp */bin/* bin
rm bin/*.a # delete un-needed static libs
#rm -rf */bin */build #delete intermediate build folders
