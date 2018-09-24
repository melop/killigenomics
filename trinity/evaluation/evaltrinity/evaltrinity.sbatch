#!/bin/bash
#SBATCH -p hugemem,himem,blade
#SBATCH -n 1
#SBATCH -c 40

export PERL5LIB=/perl5/share/perl5/:$PERL5LIB #Add lib5 lib
export PATH=./software/source/trinityrnaseq-2.1.1/
TRINITY_HOME=./software/source/trinityrnaseq-2.1.1/ #path to trinity executable
Assembly=Trinity.fasta #trinity assembly
READPFX=NOR_paired #read prefix
ESTMETHOD=RSEM
THREADS=40

INFO="#############################################
### Job array example - templer@age.mpg.de
### date $(date)
### hostname $(hostname)
#############################################"

echo -e "$INFO" 1>&2


hostname


module load Bowtie/1.1.2
module load Bowtie2/2.2.6
module load Java

perlbrew switch 5.22.1t
ln -s ../../trinityassem1/${Assembly} ./${Assembly}
bash lnk.sh
#$TRINITY_HOME/util/align_and_estimate_abundance.pl --output_dir $READPFX --transcripts $Assembly --est_method $ESTMETHOD --aln_method bowtie2 --trinity_mode --prep_reference --thread_count $THREADS
$TRINITY_HOME/util/align_and_estimate_abundance.pl --output_dir $READPFX --transcripts $Assembly --seqType fq --left ${READPFX}_1.fq.gz --right ${READPFX}_2.fq.gz --est_method $ESTMETHOD --aln_method bowtie --trinity_mode --prep_reference --thread_count $THREADS
ln -s RSEM.isoforms.results $READPFX/RSEM.isoforms.results2
$TRINITY_HOME/util/abundance_estimates_to_matrix.pl --est_method $ESTMETHOD --out_prefix trans_counts  $READPFX/RSEM.isoforms.results $READPFX/RSEM.isoforms.results2
$TRINITY_HOME/util/misc/contig_ExN50_statistic.pl trans_counts.TMM.EXPR.matrix $Assembly | tee ExN50.stats

exit

