#!/bin/bash
#SBATCH -p blade,himem,hugemem
#SBATCH -n 1 
#SBATCH -c 40

THREADS=40
OUTDIR=BEEST_RNA_1
ASSEMBLY=./Gapfilled_2/OUTPUT/OUTPUT.gapfilled.final.fa #the gap filled assembly to be scaffolded by RNAseq data
THREADS=40
BESST=./BESST_RNA/src/Main.py #full path to besstrna program
MAPSCRIPT=./reads_to_ctg_map.ray.py #mapping script. these changes have been incorporated in the latest version of BESST.
tmppathcmd="--tmp_path ./mapped/ --norebuildindex" #temporary folder
BWA=/software/bwa/0.6.2/bwa     #this bwa is used for the mate pairs due to short reads. the aln/sampe method is used
BWAMEM=/software/BWA/0.7.12/bwa #this is for bwa-mem 

EXCLUDE_LEN=1000 #EXCLUDE CONTIGS BELOW THIS LENGTH BEFORE SCAFFOLDING. THESE ARE ALSO REMOVED IN THE FINAL RESULTS!
MIN_LINK=4


source ~/.bashrc

module load R
module load STAR/2.5.1b


mkdir -p $OUTDIR;

cd $OUTDIR;

rm ./discovar.assembly.fa
ln -s `realpath $ASSEMBLY` ./discovar.assembly.fa


function MapMPReads {

	READ1=$1
	READ2=$2
	
	OUTBAM=`basename ${READ1}`
	if [ -e  ${OUTBAM}.ok ]; then
		echo ${OUTBAM}Aligned.sortedByCoord.out.bam
	else

		echo ${OUTBAM}Aligned.sortedByCoord.out.bam
	fi
}  

function MergeBam {
     MergeBamName=$1
     BamList=$2
     if [ -e  $MergeBamName ]; then
          echo $MergeBamName is present, skip merge...
     else
          echo merging $MergeBamName from $BamList ... 
          samtools merge $MergeBamName  $BamList 
     fi

     if [ -e  ${MergeBamName}.bai ]; then
          echo ${MergeBamName} already indexed ;
     else
          echo Indexing ${MergeBamName} ...
          samtools index ${MergeBamName}
     fi
}


arrR1=( `ls ../rnaseq_reads/*1.f*q*` )
arrR2=( `ls ../rnaseq_reads/*2.f*q*` )

sBamLIST=""
for i in "${!arrR1[@]}"; do 
	sR1="${arrR1[$i]}"
	sR2="${arrR2[$i]}"
	sBAMNAME=`basename $sR1`
	sBamLIST="$sBamLIST $sBAMNAME.bam"

     if [ -e  $sBAMNAME.bam.bai ]; then
	echo "$sR1 and $R2 already mapped skip..."
     else
	echo "Mapping $sR1 and $R2..."
	python $MAPSCRIPT  --bwa_path $BWAMEM --threads $THREADS --clear $tmppathcmd $sR1 $sR2 discovar.assembly.fa $sBAMNAME
     fi
done





MergeBam merged.bam "$sBamLIST" &

wait #wait for bam merge to complete

 
python $BESST 1 -c discovar.assembly.fa --mapq  30 -f merged.bam -o scaffold/  -e $MIN_LINK  -T 150000 -k $EXCLUDE_LEN -z 1000


