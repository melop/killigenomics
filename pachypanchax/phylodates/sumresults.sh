#!/bin/bash
nCount=1

#for d in output/*_RELAX.txt ; do
 #   echo "$d"
  #  if [ "$nCount" = "1" ]; then
#	echo first
#	cat $d > ./ret_RELAX.txt
#	nCount=0
 #   else
  #     	tail -n+2 $d >> ./ret_RELAX.txt
   # fi
#done


cat output/*_cleanDNA.fasta > ret_improved_cleanDNA.fasta
cat output/*_AA.fasta > ret_AA.fasta

