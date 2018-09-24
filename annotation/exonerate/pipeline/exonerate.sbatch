#!/bin/bash
#SBATCH -p blade
#SBATCH -o slurm-%a.out
#SBATCH -e slurm-%a.err
#SBATCH -a 0-1999

#SET THIS TO THE CORRECT NUMBER
TOTALPARTS=2000


module load BLAST
module load HMMER/2.3.2
module load php

#-i The transcripts to align
#-a the annotation of the transcript, marking the region of the CDS
#-r the genome assembly to align to 
#-w output folder
# check the blastandexonerate.php for other parameters
php blastandexonerate.php  -i ../Trinity.cdhitout.filtered.fa -a ../Trinity.fasta.transdecoder.filtered.annot -r ../NOR_gapfilled1.fa  -w out/ -o yes -n ${SLURM_ARRAY_TASK_ID} -t ${TOTALPARTS} # -n ${SLURM_ARRAY_TASK_ID}

exit

