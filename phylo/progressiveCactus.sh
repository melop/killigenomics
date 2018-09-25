#!/bin/bash
#SBATCH -c 40
#SBATCH -p blade

export PATH=/beegfs/group_dv/software/source/progressiveCactus/bin/:$PATH
source /beegfs/group_dv/software/source/progressiveCactus/environment

runProgressiveCactus.sh --maxThreads 40 --over taxondef.txt ./ ./progressiveCactusAln_5spp.hal
