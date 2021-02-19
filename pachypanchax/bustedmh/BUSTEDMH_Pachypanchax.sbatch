#!/bin/bash
#SBATCH -p blade
#SBATCH -o slurmlog/RelaxPachy_slurm-%a.out
#SBATCH -e slurmlog/RelaxPachy_slurm-%a.err
#SBATCH -a 0-199
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2

source ~/.bashrc
module load php

source ~/cluster*/.bash*
module load gcc/6.3.0

#SET THIS TO THE CORRECT NUMBER
TreeFolder="../export_aln_plp_2019_5_28/Relax_Pachypanchax"
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

php bustedmh.php -N ${TOTALPARTS} -f ${SLURM_ARRAY_TASK_ID} -T $TreeFolder -o $OutDIR

exit

