#!/bin/bash
#SBATCH -p blade,himem
#SBATCH -o slurm2-%a.out
#SBATCH -e slurm2-%a.err
#SBATCH -a 0-999
#SBATCH -n 1
#SBATCH -c 1

#SET THIS TO THE CORRECT NUMBER
TOTALPARTS=1000


module load BLAST
module load HMMER/2.3.2
module load php

bash doannotate_parallel.sh $TOTALPARTS ${SLURM_ARRAY_TASK_ID}


exit

