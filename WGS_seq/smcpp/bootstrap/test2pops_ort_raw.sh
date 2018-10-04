#!/bin/bash
#SBATCH -p blade
#SBATCH -n 1 
#SBATCH -c 40

export PYTHONPATH=
source /beegfs/group_dv/software/source/python_virtualenv/python3.5/bin/activate

SP=ORT

sSinglePopRetRoot=est/
sOutputDIRRoot=est_2pops/
nChr=19
mu=2.63E-9

POP1=${SP}WET
POP2=${SP}DRY

source param_2pops.sh
mkdir -p ${OUTDIR}


for sRepDir in $INPUTROOT/${SP}/bootstrap/*; do
	nRep=`basename $sRepDir`
	sOutputDIR=${sOutputDIRRoot}/${SP}/${nRep}/
	sSinglePopRetDIR1=${sSinglePopRetRoot}/${POP1}/${nRep}/
	sSinglePopRetDIR2=${sSinglePopRetRoot}/${POP2}/${nRep}/

	echo =======doing replicate $nRep ===============
	echo $sRepDir
	echo $sOutputDIR
	echo $sSinglePopRetDIR1
	echo $sSinglePopRetDIR2
	echo ============================================
	#exit

	mkdir -p $sOutputDIR

        sDoneFlag=$sOutputDIR/chrALL/model.final.json;

	if [ -e $sDoneFlag ]; then
		echo chr$i skipped, already finished.
	else 


		smc++ split  $PARAM -o $sOutputDIR/chrALL  ${sSinglePopRetDIR1}/chrALL/model.final.json  ${sSinglePopRetDIR2}/chrALL/model.final.json \
		$sRepDir/chr*/smc*.gz  \
		> ${sOutputDIR}/${SP}.log 2>&1 && touch $sDoneFlag

		if [ ! -e $sDoneFlag ]; then
			echo Error occurred, exit.
			exit 10;
		fi

		smc++ plot ${sOutputDIR}/chrALL/plot.png  ${sOutputDIR}/chrALL/model.final.json

		
	fi

done

exit



