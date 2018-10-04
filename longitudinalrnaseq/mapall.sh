#!/bin/bash
#SBATCH -p blade
#SBATCH -o slurmlog/star-%a.out
#SBATCH -e slurmlog/star-%a.err
#SBATCH -a 0-101
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

STAR=/beegfs/group_dv/software/source/STAR-2.6.0c/bin/Linux_x86_64_static/STAR
outdir=mapped
tmpdir=tmp

mkdir -p $outdir
mkdir -p $tmpdir

outdir=`realpath $outdir`
tmpdir=`realpath $tmpdir`

CURRJOB=0

for i in `realpath ./download/data`/*.gz; do
	if (( CURRJOB == SLURM_ARRAY_TASK_ID )); then

	sFqFileName=`basename $i`;
	sBamFileName=${sFqFileName/.fq.gz/}

	echo $sBamFileName;

	rm -rf $tmpdir/${sBamFileName}
	mkdir -p $outdir/${sBamFileName}/
	sDoneFlag=$outdir/${sBamFileName}/all.done;

	if [ -e $sDoneFlag ]; then

		echo $sBamFileName already finished;
	else 
		$STAR	--sjdbGTFfeatureExon exon \
			--sjdbGTFtagExonParentTranscript Parent  \
			--sjdbGTFtagExonParentGene ID \
			--runThreadN 4 \
			--genomeDir /beegfs/group_dv/home/RCui/killifish_genomes/jenalongitudinalrnaseq/NORref \
			--readFilesIn $i \
			--readFilesCommand 'zcat' \
			--sjdbGTFfile genes.gff \
			--outFileNamePrefix $outdir/${sBamFileName}/ \
			--outSAMmapqUnique 200 \
			--outSAMtype BAM SortedByCoordinate \
			--outTmpDir $tmpdir/${sBamFileName}/ \
	 		--genomeLoad NoSharedMemory && touch $sDoneFlag;
		
	fi

		exit 0;
	fi

	CURRJOB=$(( CURRJOB + 1 ));
done
