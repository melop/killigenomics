#!/usr/bin/bash
#SBATCH -n 1
#SBATCH -c 40
#SBATCH -p blade

BUSCO=/software/source/busco2beta4/BUSCO.py
LIBROOT=/software/source/busco2beta4/libraries/actinopterygii_odb9
ASSEMBLY=`realpath transcripts.fa` #transcript file to be evaluated
OUTDIR=output
THREADS=40

export PATH=/software/source/hmmer-3.1b2-linux-intel-x86_64/binaries:$PATH; #path of hmmer3
mkdir -p $OUTDIR
cd $OUTDIR

module load Python/3.4.2 
module load BLAST/2.2.31+


python $BUSCO -i $ASSEMBLY -o $OUTDIR -l $LIBROOT -m tran -c $THREADS 
