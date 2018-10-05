STAR=/beegfs/group_dv/software/source/STAR-STAR_2.4.2a/bin/STAR
OUTDIR=./CTO
genome=CTOv1.2.fa
GTF=CTOv1.2.gff

mkdir $OUTDIR
$STAR \
--runThreadN 10 \
--runMode genomeGenerate \
--genomeDir $OUTDIR \
--genomeFastaFiles $genome \
--sjdbGTFfile $GTF \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbGTFtagExonParentGene Parent \
--sjdbGTFfeatureExon exon \
--genomeChrBinNbits 14.38532 


