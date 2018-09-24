#!/bin/bash
#SBATCH -p blade
#SBATCH -o logs/improvefragslurm-%a.out
#SBATCH -e logs/improvefragslurm-%a.err
#SBATCH -a 0-39

#SET THIS TO THE CORRECT NUMBER
TOTALPARTS=40


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
module load BLAST/2.2.31+

php improve_frag_models.php -N ${TOTALPARTS} -f ${SLURM_ARRAY_TASK_ID}

exit

