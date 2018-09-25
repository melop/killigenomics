#!/bin/bash
#SBATCH -p blade
#SBATCH -o slurm-%a.out
#SBATCH -e slurm-%a.err
#SBATCH -a 0-4
#SBATCH -n 1
#SBATCH -c 40

#SET THIS TO THE CORRECT NUMBER
TOTALPARTS=5


source ~/.bashrc

module load php
#PROT=allproteins.fa

php parallel_raxml.php  -N ${TOTALPARTS} -f ${SLURM_ARRAY_TASK_ID} -L cleaned_fa_list.txt -D ClusteRs/

#UPhO.py -in *.tre -m 5 -S 0.95 -ouT -iP -R../$PROT



