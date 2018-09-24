#!/usr/bin/bash

samtools faidx contig.orig.fa
cut -f1-2 contig.orig.fa.fai > original_congitlen.txt
sort -k 2,2nr original_congitlen.txt > original_congitlen.sorted.txt
cut -f1 original_congitlen.sorted.txt > original_congitnames.txt
LN=$1

sed "${LN}q;d" original_congitnames.txt
