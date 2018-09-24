#!/bin/bash
#SBATCH -p hugemem
#SBATCH -n 1
#SBATCH -c 64
#SBATCH --mem=1000G
source ~/.bashrc

export MALLOC_PER_THREAD=1


DiscovarDeNovo READS="sample:A :: run4*.fq.gz   +   sample:B :: run5*.fq.gz + sample:C :: run7*L001*.fq.gz + sample:D :: run7*L002*.fq.gz" OUT_DIR=RUN1 NUM_THREADS=64 MAX_MEM_GB=1000



