VCFOUT=vcfs

mkdir -p $VCFOUT;

VCFOUT=`realpath $VCFOUT`
#REFSEQ=/beegfs/group_dv/home/RCui/killifish_genomes/denovo/discovardenovo/NOR/LG/v1.2/FTH.chr.fa
REFSEQ=/beegfs/group_dv/home/RCui/killifish_genomes/denovo/discovardenovo/NOR/LG/v1.2/Both_Crosses.fasta #use NOR as ref. use folded alleles instead
LNBAMS=`realpath lnbams/`

#http://www.htslib.org/workflow/#mapping_to_variant

SAMTOOLS=/beegfs/group_dv/software/source/samtools-1.6/bin/bin/samtools
BCFTOOLS=/beegfs/group_dv/software/source/samtools-1.6/bin/bin/bcftools
arrPops=( ORTWET ORTDRY RACWET RACDRY )
MINMAPQ=30


for i in "${arrPops[@]}";do
	#call variants using the FTH outgroup as the reference.
	BAMDIR=${LNBAMS}/${i}
	cd $BAMDIR
	$BCFTOOLS mpileup -q $MINMAPQ -Ou -f $REFSEQ *.bam | $BCFTOOLS call -vmO z -o ${VCFOUT}/${i}.var.vcf.gz  > ${VCFOUT}/${i}.call.log 2>&1 &

done
wait
