#!/bin/bash
#SBATCH -p blade
#SBATCH -c 40

module load Python/2.7.11

export PATH=/beegfs/group_dv/software/source/progressiveCactus2/progressiveCactus/bin/:$PATH
source /beegfs/group_dv/software/source/progressiveCactus2/progressiveCactus/environment

runProgressiveCactus.sh --maxThreads 40 taxondef.txt ./ ./progressiveCactusAln_9spp.hal


