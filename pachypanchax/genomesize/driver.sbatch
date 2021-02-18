#!/bin/bash
#SBATCH -p blade
#SBATCH -n 1 
#SBATCH -c 40

module load hhvm
module load R
module load Python

hhvm getgenomesize.php


