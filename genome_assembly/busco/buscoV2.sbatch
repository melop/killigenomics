#!/usr/bin/bash
#SBATCH -n 1
#SBATCH -c 40
#SBATCH -p blade

BUSCO=./busco2beta4/BUSCO.py #PATH TO BUSCO SCRIPT
LIBROOT=./busco2beta4/libraries/actinopterygii_odb9	#CHOOSE YOUR BUSCO LIBRARY
ASSEMBLY=./release/NOR/1.0/NORv1.0.fa #Which assembly to evaluate?
OUTDIR=v2NOR_v1.0 #Output folder
THREADS=40

export PATH=/source/hmmer-3.1b2-linux-intel-x86_64/binaries:$PATH; #make sure hmmer3.1b2 is in the path
mkdir -p $OUTDIR
cd $OUTDIR

module load Python/3.4.2 
module load BLAST/2.2.31+


python $BUSCO -i $ASSEMBLY -o $OUTDIR -l $LIBROOT -m genome -c $THREADS -sp NORgapfilled1 #sp is the species profile, trained for the species

