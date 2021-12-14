#!/bin/sh

dir_path="./Input_files"

files=`find $dir_path -maxdepth 1 -type f`

for file in $files;
do
	echo $file
	./calcPSF.o < $file
done

