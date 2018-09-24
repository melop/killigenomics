#!/bin/bash
#SBATCH -p blade
#SBATCH -n 1
#SBATCH -c 5

#SET THIS TO THE CORRECT NUMBER
TOTALPARTS=1000

INFO="#############################################
### Job array example - templer@age.mpg.de
### date $(date)
### hostname $(hostname)
### array ID ${SLURM_ARRAY_ID}
### task ID  ${SLURM_ARRAY_TASK_ID}
#############################################"

echo -e "$INFO" 1>&2


hostname


module load hhvm
module load BLAST/2.2.31+



hhvm replace_frag_with_proteinaln.php -F -N ${TOTALPARTS} -f 120 > logs/replacefragslurm-120.out 2> logs/replacefragslurm-120.err &
hhvm replace_frag_with_proteinaln.php -F -N ${TOTALPARTS} -f 181 > logs/replacefragslurm-181.out 2> logs/replacefragslurm-181.err &
hhvm replace_frag_with_proteinaln.php -F -N ${TOTALPARTS} -f 211 > logs/replacefragslurm-211.out 2> logs/replacefragslurm-211.err &
hhvm replace_frag_with_proteinaln.php -F -N ${TOTALPARTS} -f 694 > logs/replacefragslurm-694.out 2> logs/replacefragslurm-694.err &
hhvm replace_frag_with_proteinaln.php -F -N ${TOTALPARTS} -f 760 > logs/replacefragslurm-760.out 2> logs/replacefragslurm-760.err &

wait


exit

