#!/bin/bash
#SBATCH -p blade
#SBATCH -o slurmlog/RelaxNothos_slurm-%a.out
#SBATCH -e slurmlog/RelaxNothos_slurm-%a.err
#SBATCH -a 0-238
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2


TOTALPARTS=40

INFO="#############################################
### Job array example - templer@age.mpg.de
### date $(date)
### hostname $(hostname)
### array ID ${SLURM_ARRAY_ID}
### task ID  ${SLURM_ARRAY_TASK_ID}
#############################################"

echo -e "$INFO" 1>&2


hostname


#module load php
TOTALPARTSminus1=$(( TOTALPARTS-1 ));
mkdir logs

for i in $(seq 0 $TOTALPARTSminus1); do 
#	hhvm map.php -N ${TOTALPARTS} -f $i > logs/map${i}.log 2>&1 &
        hhvm map_nodedup.php -N ${TOTALPARTS} -f $i > logs/map${i}.log 2>&1 & #did not include deduplication because complexity is high and GATK (buggy) failed on some samples

done

wait

exit

