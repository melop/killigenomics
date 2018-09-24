#!/usr/bin/bash
#SBATCH -c 40
#SBATCH -p blade,himem,hugemem

ASSEMBLY=./final_polish/pilon_out/output/nor.pilon.rename.fa # assembly to mask
REPBASELIB=`realpath ./fullreplib.fa` #repeat library, with killifish repeats
THREADS=40

ln -sf $ASSEMBLY scf.fa


RepeatMasker -gff -html -pa $THREADS -lib $REPBASELIB  scf.fa

