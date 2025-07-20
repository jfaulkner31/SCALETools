#!/bin/bash

file1=$1
file2=$2
time1=$3
time2=$4
interptime=$5
fileloc=$6

obiwan tag -idtags='celi=1' -interptags="t=$time1" $file1
obiwan tag -idtags='celi=1' -interptags="t=$time2" $file2
obiwan interp -idtags='celi=1' -interptags="t=$interptime" $file1 $file2 > f33_interp.log
mv InterpLib.bof $fileloc
