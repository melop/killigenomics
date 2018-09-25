#!/bin/bash
#SBATCH -p blade
#SBATCH -a 0-999
#SBATCH -o logs/fasttree_slurm-%a.out
#SBATCH -e logs/fasttree_slurm-%a.err


#SET THIS TO THE CORRECT NUMBER
TOTALPARTS=1000


source ~/.bashrc

module load php


php parallel_fasttree.php  -N ${TOTALPARTS} -f ${SLURM_ARRAY_TASK_ID} -L cleaned_fa_list.txt -D ClusteRs/




