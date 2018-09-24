#!/usr/bin/bash
#SBATCH -p blade
#SBATCH -c 40

module load BLAST/2.2.31+
THREADS=40

contamsource=FTH_DDNcontig.fa

convert2blastmask -in $contamsource -parse_seqids -masking_algorithm repeat \
 -masking_options "repeatmasker, default" -outfmt maskinfo_asn1_bin \
 -out $contamsource.asnb

makeblastdb -in $contamsource -dbtype nucl -out $contamsource -mask_data $contamsource.asnb -parse_seqids

blastn -task megablast -query contigs.fa -db $contamsource -outfmt 7 -evalue 1e-95 -num_threads $THREADS > contig.megablast.txt
