#!/bin/bash
#SBATCH -p blade
#SBATCH -o logs/3refineslurm-%a.out
#SBATCH -e logs/3refineslurm-%a.err
#SBATCH -a 0-199

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


module load hhvm
module load BLAST/2.2.31+

hhvm refinemodels.php -N ${TOTALPARTS} -f ${SLURM_ARRAY_TASK_ID}

exit

