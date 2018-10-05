#!/bin/bash

arrInput=( CTO )
sourceFolder=./reads/
SCRIPT=./mapstarPE.sh
outFolder=./mappedPE/

mkdir $outFolder
for i in "${arrInput[@]}" 
do 
echo doing $i ...
bash $SCRIPT ${sourceFolder}/${i}_paired_1.fq.gz ${sourceFolder}/${i}_paired_2.fq.gz $outFolder > $outFolder/${i}_trim.log 2>&1 &
done
wait
