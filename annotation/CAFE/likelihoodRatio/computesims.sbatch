#!/bin/bash
#SBATCH -p blade
#SBATCH -o slurmlogs/slurm-%a.out
#SBATCH -e slurmlogs/slurm-%a.err
#SBATCH -a 0-999

#SET THIS TO THE CORRECT NUMBER
TOTALPARTS=1000


module load php

php computesim.php  -N ${TOTALPARTS} -f ${SLURM_ARRAY_TASK_ID} 
