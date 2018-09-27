#!/bin/bash
#SBATCH -p blade
#SBATCH -o slurmlogs/pcocslurm-%a.out
#SBATCH -e slurmlogs/pcocslurm-%a.err
#SBATCH -a 0-199
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

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

#Load php first!
module load php
module load R/3.3.1

export PYTHONPATH=
source ~/new_cluster_profiles/.bash*

module load shifter

export R_LIBS=:/beegfs/group_dv/home/RCui/R/x86_64-redhat-linux-gnu-library/3.2/:/beegfs/group_dv/home/RCui/R/x86_64-redhat-linux-gnu-library/3.2/


php pcoc.php -N ${TOTALPARTS} -f ${SLURM_ARRAY_TASK_ID}

exit

