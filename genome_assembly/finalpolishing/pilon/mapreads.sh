#!/bin/bash
#SBATCH -p blade
#SBATCH -n 1 
#SBATCH -c 40

THREADS=40
OUTDIR=mapped
ASSEMBLY=`realpath ../summed.blastmerged.fa` # the assembly from the previous step
MPREADBASE=./reads # path to the reads
THREADS=40
MAPSCRIPT=./reads_to_ctg_map.ray.py #MAPPING SCRIPT FROM BESST, this is a custom version, the functionalties had been incoroprated into the latest version of BESST
tmppathcmd="--tmp_path ./mapped/ --norebuildindex" #temporary folder 
BWA=/software/bwa/0.6.2/bwa     #this bwa is used for the mate pairs due to short reads. the aln/sampe method is used
BWAMEM=/software/BWA/0.7.12/bwa #this is for mem , for mapping PE reads
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


python $MAPSCRIPT  --bwa_path $BWAMEM --threads $THREADS --clear $tmppathcmd ../reads/250PE_run4_R1.fq.gz ../reads/250PE_run4_R2.fq.gz discovar.assembly.fa 250PE_run4
python $MAPSCRIPT  --bwa_path $BWAMEM --threads $THREADS --clear $tmppathcmd ../reads/250PE_run5_R1.fq.gz ../reads/250PE_run5_R2.fq.gz discovar.assembly.fa 250PE_run5
python $MAPSCRIPT  --bwa_path $BWAMEM --threads $THREADS --clear $tmppathcmd ../reads/250PE_run7_L1_R1.fq.gz ../reads/250PE_run7_L1_R2.fq.gz discovar.assembly.fa 250PE_run7_L1
python $MAPSCRIPT  --bwa_path $BWAMEM --threads $THREADS --clear $tmppathcmd ../reads/250PE_run7_L2_R1.fq.gz ../reads/250PE_run7_L2_R2.fq.gz discovar.assembly.fa 250PE_run7_L2

#mate pairs:
s3kbMPBamsLib1=""
s8kbMPBamsLib1=""
s12kbMPBamsLib1=""




sSpCode=NOR

sLibName=lib1
sReadLen=75bp
sInsertSize=3kb

s3kbMPBamsLib1="$s3kbMPBamsLib1 $(MapMPReads $sLibName $sReadLen $sInsertSize A $sSpCode)"
s3kbMPBamsLib1="$s3kbMPBamsLib1 $(MapMPReads $sLibName $sReadLen $sInsertSize B $sSpCode)"
s3kbMPBamsLib1="$s3kbMPBamsLib1 $(MapMPReads $sLibName $sReadLen $sInsertSize C $sSpCode)"


sInsertSize=8kb

s8kbMPBamsLib1="$s8kbMPBamsLib1 $(MapMPReads $sLibName $sReadLen $sInsertSize A $sSpCode)"
s8kbMPBamsLib1="$s8kbMPBamsLib1 $(MapMPReads $sLibName $sReadLen $sInsertSize B $sSpCode)"
s8kbMPBamsLib1="$s8kbMPBamsLib1 $(MapMPReads $sLibName $sReadLen $sInsertSize C $sSpCode)"


sInsertSize=12kb

s12kbMPBamsLib1="$s12kbMPBamsLib1 $(MapMPReads $sLibName $sReadLen $sInsertSize A $sSpCode)"
s12kbMPBamsLib1="$s12kbMPBamsLib1 $(MapMPReads $sLibName $sReadLen $sInsertSize B $sSpCode)"
s12kbMPBamsLib1="$s12kbMPBamsLib1 $(MapMPReads $sLibName $sReadLen $sInsertSize C $sSpCode)"


sReadLen=150bp
sInsertSize=3kb

s3kbMPBamsLib1="$s3kbMPBamsLib1 $(MapMPReads $sLibName $sReadLen $sInsertSize A $sSpCode)"
s3kbMPBamsLib1="$s3kbMPBamsLib1 $(MapMPReads $sLibName $sReadLen $sInsertSize B $sSpCode)"
s3kbMPBamsLib1="$s3kbMPBamsLib1 $(MapMPReads $sLibName $sReadLen $sInsertSize C $sSpCode)"



sInsertSize=8kb

s8kbMPBamsLib1="$s8kbMPBamsLib1 $(MapMPReads $sLibName $sReadLen $sInsertSize A $sSpCode)"
s8kbMPBamsLib1="$s8kbMPBamsLib1 $(MapMPReads $sLibName $sReadLen $sInsertSize B $sSpCode)"
s8kbMPBamsLib1="$s8kbMPBamsLib1 $(MapMPReads $sLibName $sReadLen $sInsertSize C $sSpCode)"


sInsertSize=12kb

s12kbMPBamsLib1="$s12kbMPBamsLib1 $(MapMPReads $sLibName $sReadLen $sInsertSize A $sSpCode)"
s12kbMPBamsLib1="$s12kbMPBamsLib1 $(MapMPReads $sLibName $sReadLen $sInsertSize B $sSpCode)"
s12kbMPBamsLib1="$s12kbMPBamsLib1 $(MapMPReads $sLibName $sReadLen $sInsertSize C $sSpCode)"



sLibName=lib2
sReadLen=150bp

sInsertSize=8kb

s8kbMPBamsLib2="$s8kbMPBamsLib2 $(MapMPReads $sLibName $sReadLen $sInsertSize A $sSpCode)"
s8kbMPBamsLib2="$s8kbMPBamsLib2 $(MapMPReads $sLibName $sReadLen $sInsertSize B $sSpCode)"
s8kbMPBamsLib2="$s8kbMPBamsLib2 $(MapMPReads $sLibName $sReadLen $sInsertSize C $sSpCode)"


sInsertSize=12kb

s12kbMPBamsLib2="$s12kbMPBamsLib2 $(MapMPReads $sLibName $sReadLen $sInsertSize A $sSpCode)"
s12kbMPBamsLib2="$s12kbMPBamsLib2 $(MapMPReads $sLibName $sReadLen $sInsertSize B $sSpCode)"
s12kbMPBamsLib2="$s12kbMPBamsLib2 $(MapMPReads $sLibName $sReadLen $sInsertSize C $sSpCode)"

MergeBam 250PE.lib5.merged.bam "250PE_run5.bam  250PE_run7_L1.bam 250PE_run7_L2.bam" &
MergeBam 3kbMP.lib1.merged.bam "$s3kbMPBamsLib1" &
MergeBam 8kbMP.lib1.merged.bam "$s8kbMPBamsLib1" &
MergeBam 12kbMP.lib1.merged.bam "$s12kbMPBamsLib1" &
MergeBam 8kbMP.lib2.merged.bam "$s8kbMPBamsLib2" &
MergeBam 12kbMP.lib2.merged.bam "$s12kbMPBamsLib2" &


wait #wait for bam merge to complete



