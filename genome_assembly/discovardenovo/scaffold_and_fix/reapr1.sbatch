#!/bin/bash
#SBATCH -p blade
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem 50G


OUTDIR=Reapr_1 #Output folder of the broken scaffolds
ASSEMBLY=./Gapfilled_1/OUTPUT/OUTPUT.gapfilled.final.fa #gapfilled alignments to break
MPREADBASE=./reads/ #mate pair reads
THREADS=40

source ~/.bashrc

function MapMPReads {
	R1FILE=$1
	R1BASENAME=`basename $R1FILE`
	READDIR=`dirname  $R1FILE`
	R2BASENAME="${R1BASENAME/_R1/_R2}"

	OUTBAM=${R1BASENAME}.bam
	if [ -e  ${OUTBAM} ]; then
		echo ${OUTBAM}
	else
		reapr smaltmap -n $THREADS assembly.renamed.fa $READDIR/$R1BASENAME $READDIR/$R2BASENAME $OUTBAM >> smaltmap.log 2>&1
		echo ${OUTBAM}
	fi
}  


module load R

mkdir -p $OUTDIR;

cd $OUTDIR;

reapr facheck ${ASSEMBLY} assembly.renamed
arrMPs="$MPREADBASE/MP*lib2/*12kb*R1.fastq.revcomp.fq"; #BE CAREFUL! repear requires the mate pair reads to be in FR direction; and it only accepts 1 library, usually the longest one.

sAllBams="";

for f in $arrMPs; do
	sAllBams="$sAllBams $(MapMPReads $f)";
done

samtools1.2 merge alljmp.12kb.bam $sAllBams

reapr pipeline assembly.renamed.fa alljmp.12kb.bam output_12kb 



