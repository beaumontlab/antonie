#!/bin/bash

rm -f data.js loci.0.js
./antonie -1 sbw25/P1-1-35_S5_L001_R1_001.fastq -2 sbw25/P1-1-35_S5_L001_R2_001.fastq -r sbw25/NC_012660.fna -a sbw25/NC_012660.gbk 

if grep 895351 data.js  -q
then
	echo Found SNP 895151
else
	echo missed SNP 895151
	exit 1
fi

if grep YP_002874253.1 data.js -q
then
	echo Found YP_002874253.1
else
	echo missed YP_002874253.1
	exit 1
fi


exit 0
