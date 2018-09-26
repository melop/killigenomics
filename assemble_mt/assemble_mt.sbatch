#!/bin/bash
#SBATCH -n 1
#SBATCH -c 2
#SBATCH --mem=1G
#SBATCH -p blade

module load hhvm
hhvm novoplasty.php
