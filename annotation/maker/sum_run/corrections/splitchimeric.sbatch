#!/bin/bash
#SBATCH -p blade
#SBATCH -o logs/split_slurm-%a.out
#SBATCH -e logs/split_slurm-%a.err
#SBATCH -a 0-199
#SBATCH -c 1
#SBATCH --mem=250G

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


#module load php/7.0.12
module load hhvm
module load BLAST/2.2.31+
module load R

sleep $[ ( $RANDOM % 10 )  + 1 ]s
hhvm splitchimeric.php  -N ${TOTALPARTS} -f ${SLURM_ARRAY_TASK_ID}

exit

