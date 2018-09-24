#!/bin/bash
#SBATCH -p himem,hugemem,blade
#SBATCH -n 1
#SBATCH -c 40

export PERL5LIB=~/perl5/share/perl5/:$PERL5LIB
export PATH=$PATH:/source/trinityrnaseq-2.1.1/ #Add the trinity executable path


hostname


module load Bowtie/1.1.2
module load Java

perlbrew switch 5.22.1t


Trinity --grid_conf cluster.config --seqType fq --max_memory 250G --left NOR_paired_1.fq.gz --right NOR_paired_2.fq.gz --CPU 40 --output trinityassem1 

exit

