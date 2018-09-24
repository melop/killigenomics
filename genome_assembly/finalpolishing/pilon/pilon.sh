#!/bin/bash
#SBATCH -p hugemem
#SBATCH -n 1 
#SBATCH -c 64
#SBATCH --mem=1000G

THREADS=40
BAMDIR=`realpath mapped`
BAMDIRFIXED=`realpath fixedbam`
PILON=./pilon/pilon-1.22.jar #PATH TO THE PILON EXECUTABLE

export JAVA_HOME=/jre1.8.0_131 #WHERE IS JAVA?
export PATH=/jre1.8.0_131/bin/:$PATH; #WHERE IS JAVA BIN?
module load BWA/0.7.12
module load R


OUTDIR=pilon_out

mkdir -p $OUTDIR
cd $OUTDIR

ln -sf $BAMDIR/discovar.assembly.fa ./ref.fa

java -Xmx1000G -jar $PILON --genome ref.fa \
			--frags $BAMDIR/mergedPE.bam \
			--bam $BAMDIRFIXED/3kbMP.lib1.merged.fixflag.bam \
			--bam $BAMDIRFIXED/8kbMP.lib2.merged.fixflag.bam \
			--bam $BAMDIRFIXED/12kbMP.lib1.merged.fixflag.bam \
                        --bam $BAMDIRFIXED/12kbMP.lib2.merged.fixflag.bam \
			--output pilon --outdir output --changes --tracks --diploid --fix all,breaks --threads $THREADS --mindepth 0.25 

