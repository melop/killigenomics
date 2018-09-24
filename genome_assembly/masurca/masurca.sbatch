#!/bin/bash
#SBATCH -p blade
#SBATCH -n 1
#SBATCH -c 40

MASURCAROOT=./MaSuRCA-3.2.1_08102016/bin #the root to the MasurCA executable

INFO="#############################################
### Job array example - templer@age.mpg.de
### date $(date)
### hostname $(hostname)
#############################################"

echo -e "$INFO" 1>&2


hostname


$MASURCAROOT/masurca masurca.config
bash assemble.sh

exit
