STAR=/beegfs/group_dv/software/source/STAR-STAR_2.4.2a/bin/STAR
GENOMEDIR=./CTO
genome=CTOv1.2.fa

FILE1NAME=`basename $1`

mkdir -p ${3}/${FILE1NAME}_out/

$STAR \
--runThreadN 20 \
--genomeDir $GENOMEDIR \
--readFilesIn $1 $2 \
--readFilesCommand zcat \
--sjdbGTFfile CTOv1.2.gff \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbGTFtagExonParentGene Parent \
--sjdbGTFfeatureExon exon \
--outFileNamePrefix ${3}/${FILE1NAME}_out/ \
--outSAMmapqUnique 254 \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts

