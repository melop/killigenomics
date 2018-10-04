mkdir vcfs
#http://www.htslib.org/workflow/#mapping_to_variant

SAMTOOLS=/beegfs/group_dv/software/source/samtools-1.6/bin/bin/samtools
BCFTOOLS=/beegfs/group_dv/software/source/samtools-1.6/bin/bin/bcftools
arrPops=( ORTWET ORTDRY )
MINMAPQ=60


for i in "${arrPops[@]}";do
	#call variants using the FTH outgroup as the reference.
	$BCFTOOLS mpileup -q $MINMAPQ -Ou -f /beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/SimNORLG1/ref.fa mapped/bam/${i}*.bam | $BCFTOOLS call -vmO z -o vcfs/${i}.var.vcf.gz  > vcfs/${i}.call.log 2>&1 &

done
wait
