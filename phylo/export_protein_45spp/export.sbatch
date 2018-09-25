#!/bin/bash
#SBATCH -p blade
#SBATCH -o slurm-%a.out
#SBATCH -e slurm-%a.err
#SBATCH -a 0-99

#SET THIS TO THE CORRECT NUMBER
TOTALPARTS=100

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
module load fastahack
mkdir -p output
php export_coding_seq.php -N ${TOTALPARTS} -f ${SLURM_ARRAY_TASK_ID}

exit

