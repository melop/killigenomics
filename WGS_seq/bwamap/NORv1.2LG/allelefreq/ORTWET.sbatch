#!/bin/bash
#SBATCH -p blade,himem,hugemem
#SBATCH -o slurmlog/ORTWET-%a.out
#SBATCH -e slurmlog/ORTWET-%a.err
#SBATCH -a 1-100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2


TOTALPARTS=100
POP=ORTWET
INFO="#############################################
### Job array example - templer@age.mpg.de
### date $(date)
### hostname $(hostname)
### array ID ${SLURM_ARRAY_ID}
### task ID  ${SLURM_ARRAY_TASK_ID}
#############################################"

echo -e "$INFO" 1>&2


hostname
cd $POP
../../../../cpp/Package-GFE/GFE_ray -in My_in_GFE.txt -mode f -out out_f_GFE -thispart ${SLURM_ARRAY_TASK_ID} -totalpart ${TOTALPARTS} 

exit

