#!/bin/bash
#SBATCH -p blade
#SBATCH -n 1 
#SBATCH -c 40

export PYTHONPATH=
source /beegfs/group_dv/software/source/python_virtualenv/python3.5/bin/activate

PARAM=" --cores 10 --timepoints h  --hs 2"
sPop="ORTDRY";
sOutputDIRRoot=est/${sPop}
nChr=19
mu=2.63E-9

source params.sh


for sRepDir in $INPUTROOT/${sPop}/bootstrap/*; do
	nRep=`basename $sRepDir`
	sOutputDIR=${sOutputDIRRoot}/${nRep}/

	mkdir -p $sOutputDIR

        sDoneFlag=$sOutputDIR/chrALL/model.final.json;

	if [ -e $sDoneFlag ]; then
		echo chr$i skipped, already finished.
	else 

		smc++ estimate $PARAM  -o $sOutputDIR/chrALL $mu $sRepDir/chr*/smc.gz >  $sOutputDIR/chrALL.log.txt 2>&1
		smc++ plot $sOutputDIR/chrALL/plot.png $sOutputDIR/chrALL/model.final.json
		
	fi

done

exit
