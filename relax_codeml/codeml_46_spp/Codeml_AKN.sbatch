#!/bin/bash
#SBATCH -p blade
#SBATCH -o slurmlog/CodemlAKN_slurm-%a.out
#SBATCH -e slurmlog/CodemlAKN_slurm-%a.err
#SBATCH -a 0-199
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

#SET THIS TO THE CORRECT NUMBER
TreeFolder="../export_protein_multiref_23_july_2017/Codeml_AKN"
OutDIR=`basename $TreeFolder`;

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
mkdir -p $OutDIR

php codeml.php -N ${TOTALPARTS} -f ${SLURM_ARRAY_TASK_ID} -T $TreeFolder -o $OutDIR

exit

