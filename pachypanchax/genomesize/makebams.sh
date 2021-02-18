#!/bin/bash
#SBATCH -p blade
#SBATCH -n 1 
#SBATCH -c 40


OUTDIR=$1
ASSEMBLY=$2
THREADS=$3
R1=`realpath $4`
R2=`realpath $5`
SPCODE=$6

MAPSCRIPT=/beegfs/group_dv/software/source/BESST/scripts/reads_to_ctg_map.ray.py
tmppathcmd="--tmp_path ./tmp/ --norebuildindex"
BWA=/software/bwa/0.6.2/bwa     #this bwa is used for the mate pairs due to short reads. the aln/sampe method is used
BWAMEM=/software/BWA/0.7.12/bwa #this is for mem 


source ~/.bashrc

module load R
module load bwa/0.6.2
module load BWA/0.7.12


mkdir -p $OUTDIR;

cd $OUTDIR;

ln -s `realpath $ASSEMBLY` ./ref.fa

function MapMPReads {
	LIBNAME=$1
	READLEN=$2
	INSERTSIZE=$3
	LIBTYPE=$4
	SPECIESCODE=$5
	BWAPATH=$6
	noMEM=""

	if [[ -z "${BWAPATH// }" ]]; then
		BWAPATH=$BWA;
	fi

        if [ "$BWAPATH" == "$BWA" ]; then
		noMEM="--nomem"
	fi

	OUTBAM=${READLEN}_${LIBNAME}_${INSERTSIZE}_${LIBTYPE}.bam
	if [ -e  ${OUTBAM} ]; then
		echo ${OUTBAM}
	else
		python $MAPSCRIPT --bwa_path $BWAPATH $noMEM --threads $THREADS --clear $tmppathcmd ../reads/MP_${READLEN}_${LIBNAME}/${SPECIESCODE}*${INSERTSIZE}*_${LIBTYPE}_R1.fastq ../reads/MP_${READLEN}_${LIBNAME}/${SPECIESCODE}*${INSERTSIZE}*_${LIBTYPE}_R2.fastq ref.fa ${READLEN}_${LIBNAME}_${INSERTSIZE}_${LIBTYPE} >> mapscriptlog.txt 2>&1
		echo ${OUTBAM}
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


python $MAPSCRIPT  --bwa_path $BWAMEM --threads $THREADS --clear $tmppathcmd "$R1" "$R2" ref.fa $SPCODE

