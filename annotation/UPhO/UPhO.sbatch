#!/usr/bin/bash
#SBATCH -n 1
#SBATCH -c 80
#SBATCH -p blade,himem,hugemem


export PATH=/software/source/UPhO/:$PATH #add the path of UPHO executables

module load BLAST
PROT=allproteins.fa #File containing all protein sequences formatted by preprocessAA.php

Blast_helper.sh -i $PROT



