#!/bin/bash
#SBATCH -p blade
#SBATCH --mem 250G
#SBATCH -n 1 
#SBATCH -c 40
#SBATCH -o hostname_%j.out      # File to which STDOUT will be written
#SBATCH -e hostname_%j.err      # File to which STDERR will be written


ASSEMBLY=xxxx.fa #specify assembly contigs here
REDUNDANS=PATH2REDUNDANS/redundans.py #specify the redundans script here
THREADS=40 #how many threads to use?
OUTDIR=redundans_out #output folder
IDENTITY=0.8 # the minimal contig identity to merge


module load Python/2.7.9


$REDUNDANS --noscaffolding --nogapclosing --identity $IDENTITY -v -f $ASSEMBLY -o $OUTDIR -t $THREADS -q 30  --log log.txt 




