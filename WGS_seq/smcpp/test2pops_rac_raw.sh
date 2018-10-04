#!/bin/bash
#SBATCH -p himem
#SBATCH -n 1
#SBATCH -c 10
#SBATCH --mem=500G


export PYTHONPATH=
source /beegfs/group_dv/software/source/python_virtualenv/python3.5/bin/activate

SP=RAC

source param_2pops.sh
mkdir -p ${OUTDIR}

smc++ split  $PARAM -o ${OUTDIR}/${SP}  ${INDIR}/${SP}WET_unfold/chrALL/model.final.json ${INDIR}/${SP}DRY_unfold/chrALL/model.final.json \
${INPUTROOT}/${SP}/chr*/smc.*.gz \
> ${OUTDIR}/${SP}.log 2>&1



smc++ plot ${OUTDIR}/${SP}/plot.png  ${OUTDIR}/${SP}/model.final.json

