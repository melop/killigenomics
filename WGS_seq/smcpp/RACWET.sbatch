#!/bin/bash
#SBATCH -p blade
#SBATCH -n 1 
#SBATCH -c 40

export PYTHONPATH=
source /beegfs/group_dv/software/source/python_virtualenv/python3.5/bin/activate

PARAM=" --cores 10 --timepoints h  --hs 2"
sPop="RACWET";
sOutputDIR=est/${sPop}_unfold
nChr=19
mu=2.63E-9

source params.sh

mkdir -p $sOutputDIR;


	sDoneFlag=$sOutputDIR/chrALL/model.final.json;

	if [ -e $sDoneFlag ]; then
		echo chr$i skipped, already finished.
	else 

		smc++ estimate $PARAM  -o $sOutputDIR/chrALL $mu $INPUTROOT/${sPop}/chr*/smc.gz >  $sOutputDIR/chrALL.log.txt 2>&1
		smc++ plot $sOutputDIR/chrALL/plot.png $sOutputDIR/chrALL/model.final.json
		
	fi
exit

for i in $(seq 1 $nChr); do 
	sDoneFlag=$sOutputDIR/chr$i/model.final.json;

	if [ -e $sDoneFlag ]; then
		echo chr$i skipped, already finished.
	else 

                smc++ estimate  $PARAM  -o $sOutputDIR/chr$i $mu $INPUTROOT/${sPop}/chr${i}/smc.gz >  $sOutputDIR/chr$i.log 2>&1

		smc++ plot $sOutputDIR/chr$i/plot.png $sOutputDIR/chr$i/model.final.json
		
	fi
done

smc++ plot $sOutputDIR/plotAll.png $sOutputDIR/chr*/model.final.json

