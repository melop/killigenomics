###############
PICARDDIR="/beegfs/group_dv/software/source/picard-tools-1.119/"
SCRATCHDIR="./tmp/"
REF_GENOME="notfur.fa"
REF_GENOME_NAME=${REF_GENOME%.*}
SAMPLENAME=${1}
GATKPATH="/beegfs/group_dv/software/source/GATK/GenomeAnalysisTK-3.2-2/" #"/beegfs/group_dv/software/source/gatk.nightly.10.8.2015/" #  /beegfs/group_dv/software/source/gatk-source/target/
GATKPATHNEW="/beegfs/group_dv/software/source/gatk3.4.46/"
IGVTOOLSPATH="/beegfs/group_dv/software/source/IGVTools/igvtools"
KNOWNSITES="KNOWN_liftup.sorted.vcf.sorted.vcf"
FIXQUALSCORES=" "         #" -fixMisencodedQuals"
STARTFROMSTEP=1    #only start from step 7. previous steps not needed.
THREAD=3
nStep=0
SAMTOOLSCMD="samtools1.2"
BCFTOOLSCMD="bcftools1.2"



function error_exit
{
	echo "Error occured at step $1, pipeline stopped." 1>&2
	exit 1
}


#Step_1
let nStep="nStep + 1"
echo Step $nStep : Making mpileup  ...
[ $nStep -ge $STARTFROMSTEP ] && { $SAMTOOLSCMD mpileup -B --count-orphans -t DP,DV,DPR,INFO/DPR,DP4,SP --adjust-MQ 0 --min-MQ 0 -go ${SAMPLENAME}.raw.bcf -f ${REF_GENOME} ${SAMPLENAME}.sorted.bam ;

#[ $? -ne 0 ] &&  error_exit $nStep ;
 }

#Step_2
let nStep="nStep + 1"
echo Step $nStep : calling snps...
[ $nStep -ge $STARTFROMSTEP ] && { $BCFTOOLSCMD call -m -M -A  -O b -o ${SAMPLENAME}.called.bcf ${SAMPLENAME}.raw.bcf ;
[ $? -ne 0 ] &&  error_exit $nStep ;
 }

#Step_3
let nStep="nStep + 1"
echo Step $nStep index bcfs...
[ $nStep -ge $STARTFROMSTEP ] && { $BCFTOOLSCMD index ${SAMPLENAME}.called.bcf; 
[ $? -ne 0 ] &&  error_exit $nStep ;
 }

#Step_4
let nStep="nStep + 1"
echo Step $nStep reports...
[ $nStep -ge $STARTFROMSTEP ] && { $BCFTOOLSCMD stats -F ${REF_GENOME} -s - ${SAMPLENAME}.called.bcf > ${SAMPLENAME}.called.bcf.stats;
mkdir -p plots_${SAMPLENAME};
plot-vcfstats -p plots_${SAMPLENAME}/ ${SAMPLENAME}.called.bcf.stats;
 
[ $? -ne 0 ] &&  error_exit $nStep ;
 }
