#!/bin/bash
#SBATCH -p blade
#SBATCH -o slurm-%a.out
#SBATCH -e slurm-%a.err
#SBATCH -a 0-999

#SET THIS TO THE CORRECT NUMBER
TOTALPARTS=1000


source ~/.bashrc

module load php

php parallel_raxml.php  -N ${TOTALPARTS} -f ${SLURM_ARRAY_TASK_ID} -L cleaned_fa_list.txt -D ClusteRs/





