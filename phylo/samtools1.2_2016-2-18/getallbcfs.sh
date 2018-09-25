#!/bin/bash


arrInput=( Xmac Fht Alm PLP )

rm *.sorted.bam


#mkdir $outFolder
ln -s `realpath ../cactusmafs/*.sorted.bam*` ./;

for i in "${arrInput[@]}" 
do 
echo converting $i ...;

bash pipeline2.sh ${i} >> ${i}.log 2>&1 &

done
wait
