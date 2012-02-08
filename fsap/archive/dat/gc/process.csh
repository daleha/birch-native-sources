#!/bin/csh
set file = $1
set name = $file:r
#echo $file
#cat $file  | tr 'i' ' ' | tr -s " " | cut -f2,5,8,11 -d" "

cat $file | cut -c5-6,20-21,35-36,50 > $name.gc
 