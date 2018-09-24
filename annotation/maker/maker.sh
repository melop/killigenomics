#!/bin/bash
#SBATCH -n 3
#SBATCH -c 40
#SBATCH --mem=250G
#SBATCH -p blade

bash ~/.bashrc
perlbrew use 5.22.1t

#load mpich
source /beegfs/group_dv/software/mpich3.2.module

mkdir -p makertmp
#module load repeatmasker
module load BLAST
mpirun -n 3 maker  -fix_nucleotides  > out.log 2>err.log
