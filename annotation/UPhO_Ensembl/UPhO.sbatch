#!/usr/bin/bash
#SBATCH -n 1
#SBATCH -c 40
#SBATCH -p blade,himem,hugemem


module load BLAST
PROT=allproteins.fa

Blast_helper.sh -i $PROT



