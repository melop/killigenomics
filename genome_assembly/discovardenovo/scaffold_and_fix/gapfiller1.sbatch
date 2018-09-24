#!/bin/bash
#SBATCH -p blade
#SBATCH -n 1
#SBATCH -c 40

pwd

OUTDIR=Gapfilled_1 #output folder for the gapfilled scaffolds
ASSEMBLY=./BEEST_1/scaffold/BESST_output/pass7/Scaffolds_pass7.fa #the scaffolded assembly for the gapfilling
THREADS=40
READBASE=./reads #mate pair reads
GAPFILLER=./bin/GapFiller_v1-10_linux-x86_64/GapFiller.pl #Full path to the gapfiller program
READLIBFILE=gapfiller_librarydef.txt #definition of libraries
TRIMOFFCONTIG=25 #trim off this much from the edges of the scaffolds.
OVERLAP=74 #Longest read length for mate pairs is 75
GAPCLOSEMAXDIFF=500 #If gap closed is within this number of difference from the originally estimated gap size by the scaffolder, then accept.
ITERATIONS=2

source ~/.bashrc

module load R
perlbrew switch 5.22.1t

mkdir -p $OUTDIR;

cd $OUTDIR;

ln -s `realpath $ASSEMBLY` ./scfld.assembly.fa
ln -s `realpath $READBASE` ./
ln -s `realpath ../$READLIBFILE` ./$READLIBFILE

perl $GAPFILLER  -s scfld.assembly.fa -T $THREADS -l $READLIBFILE -m $OVERLAP -t $TRIMOFFCONTIG -d $GAPCLOSEMAXDIFF -i $ITERATIONS -b OUTPUT 
