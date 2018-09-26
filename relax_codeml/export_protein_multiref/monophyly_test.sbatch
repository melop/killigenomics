#!/bin/bash
#SBATCH -p blade
#SBATCH -o slurmlogs/monophy3slurm-%a.out
#SBATCH -e slurmlogs/monophy3slurm-%a.err
#SBATCH -a 0-199
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4

#SET THIS TO THE CORRECT NUMBER
TOTALPARTS=200

INFO="#############################################
### Job array example - templer@age.mpg.de
### date $(date)
### hostname $(hostname)
### array ID ${SLURM_ARRAY_ID}
### task ID  ${SLURM_ARRAY_TASK_ID}
#############################################"

echo -e "$INFO" 1>&2


hostname


module load php

php monophyly_test.php -N ${TOTALPARTS} -f ${SLURM_ARRAY_TASK_ID}

exit

