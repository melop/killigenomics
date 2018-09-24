#!/usr/bin/bash
#SBATCH -p blade,himem,hugemem
#SBATCH -c 40
#SBATCH --mem=250G

ASSEMBLY=./BEEST_RNA_5/scaffold/pass1/Scaffolds-pass1.fa # the assembly to look at
SP=NOR
REPMODEL=RepeatModeler #command of repeatmodeler
THREADS=40

ln -sf $ASSEMBLY scf.fa

$REPMODEL/BuildDatabase -name $SP -engine ncbi scf.fa

$REPMODEL/RepeatModeler -engine ncbi -pa $THREADS -database $SP
