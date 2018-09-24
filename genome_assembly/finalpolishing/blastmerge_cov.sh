#!/bin/bash
#SBATCH -p blade,himem,hugemem
#SBATCH -o slurm-%a.out
#SBATCH -e slurm-%a.err
#SBATCH -a 0-39

#SET THIS TO THE CORRECT NUMBER
nTotalPart=40;
sIn="./denovo/release/NOR/1.0/NORv1.0.fa"; #The genome assembly to close overhangs
sTmpFolder="./tmp_aggresive/"; //temporary folder

nMinBlockIdentity=90; # Minimal block identity from 0-100
nMaxHitDistance=100000; #if hits are more than 100kb away,ignore them, probably just interspersed repeats.
nMinOverallIdentity=0.95; #the min identity of the whole dup block, 0-1
nMaxGapPerc=0.15; #the maximal amount of unaligned sequence within the dup block
nMaxInterHitDistance=3000; #if the sub hits are < this bp away, merge and treat as one big block
nMaxNonNDistanceBetweenDups=20000; # the two duplicated blocks need to be <3000bp
nMaxNonNDistanceBetweenDupsPerc=5; # the two duplicated blocks need to be <500% of duplicated block length

sBAMFile="./mapped/mergedPE.bam"; #the sequencing reads mapped using BWA-MEM
nAvgGenomeCov=64.96965; #what is the average coverage in this bam file?
sRepeatGFF="./repeatmasker/NOR_v1.0/scf.fa.out.gff"; #GFF file giving repeat coordinates
nMinLenAfterRepeatExclusion=20; #at least 20bp left after masking the repeats in the overlap

INFO="#############################################
### Job array example - templer@age.mpg.de
### date $(date)
### hostname $(hostname)
### array ID ${SLURM_ARRAY_ID}
### task ID  ${SLURM_ARRAY_TASK_ID}
#############################################"

echo -e "$INFO" 1>&2


hostname


module load php
#module load BLAST
module load R
php blastmergegap_cov.php  -N ${nTotalPart} -f ${SLURM_ARRAY_TASK_ID} \
-i "$sIn" -T "$sTmpFolder" -b $nMinBlockIdentity -h $nMaxHitDistance \
-O $nMinOverallIdentity -M $nMaxGapPerc -D $nMaxInterHitDistance \
-g $nMaxNonNDistanceBetweenDups -G $nMaxNonNDistanceBetweenDupsPerc \
-B $sBAMFile -C $nAvgGenomeCov -R $sRepeatGFF -r $nMinLenAfterRepeatExclusion

exit

