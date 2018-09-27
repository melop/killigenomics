#!/bin/bash
#SBATCH -p blade
#SBATCH -o slurmlog/ORTWET_slurm-%a.out
#SBATCH -e slurmlog/ORTWET_slurm-%a.err
#SBATCH -a 0-99
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

sPop=ORTWET
PREFIX=
#PREFIX=sametheta_
#PREFIX=sametheta_nopolerr_

sleep $(( ( RANDOM % 10 )  + 1 ))
bash run.sh $sPop ${SLURM_ARRAY_TASK_ID} $PREFIX
