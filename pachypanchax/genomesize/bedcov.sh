#!/usr/bin/bash
#SBATCH --mem=500G
#SBATCH -p himem
#SBATCH -c 2

bedtools coverage -a singlecopy_coord.bed -b ./mappedbam/illuminaFish1.bam -d > perbasecov.fish1.txt
