#!/bin/bash
#SBATCH -p blade
#SBATCH -o slurmlog/fc-%a.out
#SBATCH -e slurmlog/fc-%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40

FEATURECOUNT=/beegfs/group_dv/software/source/subread-1.6.2-Linux-x86_64/bin/featureCounts
outdir=featurecount

mkdir -p lnmappedreads

for i in ./mapped/*/Aligned.sortedByCoord.out.bam; do
	sDIRname=`dirname $i`
	sSample=`basename $sDIRname`;

	ln -sf `realpath $i` lnmappedreads/${sSample}.bam

done

$FEATURECOUNT -T 40 -a genes.gff -O -o readcounts.txt -F GTF -t CDS -g Parent -Q 30 lnmappedreads/*.bam

