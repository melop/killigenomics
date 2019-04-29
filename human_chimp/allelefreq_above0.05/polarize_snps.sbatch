#!/bin/bash
#SBATCH -p blade
#SBATCH -o slurmlog/Consurf_slurm-%a.out
#SBATCH -e slurmlog/Consurf_slurm-%a.err
#SBATCH -a 0-24
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2

#SET THIS TO THE CORRECT NUMBER

TOTALPARTS=25

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

hhvm polarize_snps.php -N ${TOTALPARTS} -f ${SLURM_ARRAY_TASK_ID}

exit


