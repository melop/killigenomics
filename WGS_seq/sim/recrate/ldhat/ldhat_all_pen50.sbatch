#!/bin/bash
#SBATCH -p blade
#SBATCH -o slurmlog/ldhat50_slurm-%a.out
#SBATCH -e slurmlog/ldhat50_slurm-%a.err
#SBATCH -a 1-2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2


arrPops=( ORTWET ORTDRY )
arrChr=( `seq 1 1` )
sChrPrefix=chr

JobCount=0;

SOURCEDIR=`realpath ldhat`
OUTDIR=ldhat_pen50
SCRIPT=_run_pen50.sh

OUTDIR=`realpath $OUTDIR`
SCRIPT=`realpath $SCRIPT`

mkdir -p $OUTDIR

for i in "${arrPops[@]}";do

   for nChr in "${arrChr[@]}";do
		JobCount=$(( JobCount + 1 ));
		if (( JobCount == SLURM_ARRAY_TASK_ID )); then
			mkdir -p ${OUTDIR}/${i}/chr${nChr}/
			ln -sf $SOURCEDIR/${i}/chr${nChr}/sites.txt ${OUTDIR}/${i}/chr${nChr}/
                        ln -sf $SOURCEDIR/${i}/chr${nChr}/loc.txt ${OUTDIR}/${i}/chr${nChr}/

			echo Job $JobCount
			echo ${OUTDIR}/${i}/chr${nChr}/
			LikFile=`realpath lktables/${i}/cat_lk.txt`;
			ln -sf $SCRIPT ${OUTDIR}/${i}/chr${nChr}/run.sh;
			cd ${OUTDIR}/${i}/chr${nChr}/;
			bash run.sh $LikFile
		fi
	done
done
wait



