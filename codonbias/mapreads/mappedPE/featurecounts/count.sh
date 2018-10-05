GTF1=../../CTOv1.2.gff
BAM="../CTO_paired_1.fq.gz_out/Aligned.sortedByCoord.out.bam"

FEATURECOUNT=/beegfs/group_dv/software/source/subread-1.5.0-source/bin/featureCounts

$FEATURECOUNT -g Parent -t CDS -a $GTF1 -p -o counts.txt $BAM

