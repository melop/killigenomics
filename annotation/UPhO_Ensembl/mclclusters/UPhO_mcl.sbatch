#!/usr/bin/bash
#SBATCH -n 1
#SBATCH -c 40
#SBATCH -p blade,himem,hugemem

source ~/.bashrc

module load BLAST
PROT=allproteins.fa


cd ClusteRs/
paMATRAX+.sh -c

bash cptrees.sh
UPhO.py -in trees/*.tre -m 3 -S 0.10 -ouT -iP -R ../$PROT



