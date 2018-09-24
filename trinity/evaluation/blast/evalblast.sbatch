#!/bin/bash
#SBATCH -p hugemem,himem,blade,dontuseme
#SBATCH -n 1
#SBATCH -c 40

export PERL5LIB=/perl5/share/perl5/:$PERL5LIB #PERL 5 lib
export PATH=./trinityrnaseq-2.1.1/ #path to trinity executables
TRINITY_HOME=./trinityrnaseq-2.1.1/
Assembly=Trinity.fasta #Trinity assembly
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
module load BLAST

perlbrew switch 5.22.1t

#makeblastdb -in jenaproteins.fa -dbtype prot
ln -s ../../trinityassem1/${Assembly} ./${Assembly}
ln -s ../../../PLP/evaluation/blast/jenaprotein* ./
blastx -query $Assembly -db jenaproteins.fa -out blastx.outfmt6 \
        -evalue 1e-20 -num_threads 40 -max_target_seqs 1 -outfmt 6

$TRINITY_HOME/util/analyze_blastPlus_topHit_coverage.pl blastx.outfmt6  $Assembly jenaproteins.fa > report1.txt

$TRINITY_HOME/util/misc/blast_outfmt6_group_segments.pl \
      blastx.outfmt6   $Assembly jenaproteins.fa > blastx.outfmt6.grouped

$TRINITY_HOME/util/misc/blast_outfmt6_group_segments.tophit_coverage.pl blastx.outfmt6.grouped > report2_grouped.txt

exit

