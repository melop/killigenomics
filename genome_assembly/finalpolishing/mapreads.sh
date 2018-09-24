#!/bin/bash
#SBATCH -p blade
#SBATCH -n 1 
#SBATCH -c 40

THREADS=40
OUTDIR=mapped
ASSEMBLY=./release/NOR/1.0/NORv1.0.fa #the assembly to be fixed
MPREADBASE=./reads	#where are the reads?
THREADS=40
MAPSCRIPT=/beegfs/group_dv/software/source/BESST/scripts/reads_to_ctg_map.ray.py #mapping script
tmppathcmd="--tmp_path ./mapped/ --norebuildindex" #tmp folder for mapping
BWA=/software/bwa/0.6.2/bwa     #this bwa is used for the mate pairs due to short reads. the aln/sampe method is used
BWAMEM=/software/BWA/0.7.12/bwa #this is for mem 
EXCLUDE_LEN=1000 #EXCLUDE CONTIGS BELOW THIS LENGTH BEFORE SCAFFOLDING. THESE ARE ALSO REMOVED IN THE FINAL RESULTS!


ln -sf `realpath $MPREADBASE` ./reads
source ~/.bashrc

module load R
module load bwa/0.6.2
module load BWA/0.7.12


mkdir -p $OUTDIR;

cd $OUTDIR;

ln -s `realpath $ASSEMBLY` ./discovar.assembly.fa

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
		python $MAPSCRIPT --bwa_path $BWAPATH $noMEM --threads $THREADS --clear $tmppathcmd ../reads/MP_${READLEN}_${LIBNAME}/${SPECIESCODE}*${INSERTSIZE}*_${LIBTYPE}_R1.fastq ../reads/MP_${READLEN}_${LIBNAME}/${SPECIESCODE}*${INSERTSIZE}*_${LIBTYPE}_R2.fastq discovar.assembly.fa ${READLEN}_${LIBNAME}_${INSERTSIZE}_${LIBTYPE} >> mapscriptlog.txt 2>&1
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


#Only map Pair-end reads

python $MAPSCRIPT  --bwa_path $BWAMEM --threads $THREADS --clear $tmppathcmd ../reads/250PE_run4_R1.fq.gz ../reads/250PE_run4_R2.fq.gz discovar.assembly.fa 250PE_run4
python $MAPSCRIPT  --bwa_path $BWAMEM --threads $THREADS --clear $tmppathcmd ../reads/250PE_run5_R1.fq.gz ../reads/250PE_run5_R2.fq.gz discovar.assembly.fa 250PE_run5
python $MAPSCRIPT  --bwa_path $BWAMEM --threads $THREADS --clear $tmppathcmd ../reads/250PE_run7_L1_R1.fq.gz ../reads/250PE_run7_L1_R2.fq.gz discovar.assembly.fa 250PE_run7_L1
python $MAPSCRIPT  --bwa_path $BWAMEM --threads $THREADS --clear $tmppathcmd ../reads/250PE_run7_L2_R1.fq.gz ../reads/250PE_run7_L2_R2.fq.gz discovar.assembly.fa 250PE_run7_L2


MergeBam mergedPE.bam "250PE_run4.bam 250PE_run5.bam  250PE_run7_L1.bam 250PE_run7_L2.bam" &
wait


