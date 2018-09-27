#!/bin/bash
#SBATCH -p blade
#SBATCH -o slurmlog/ORTDRY_slurm-%a.out
#SBATCH -e slurmlog/ORTDRY_slurm-%a.err
#SBATCH -a 0-99
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

sPop=ORTDRY.excludehybrid
PREFIX=
PREFIX=sametheta_
PREFIX=sametheta_nopolerr_

bash run.sh $sPop ${SLURM_ARRAY_TASK_ID} $PREFIX
