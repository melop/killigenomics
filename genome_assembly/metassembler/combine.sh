#!/usr/bin/bash
#SBATCH -p blade
#SBATCH -c 40

## Script that uses the metassemble.py wrapper to metassemble assembly A and B, using

### Generate configuration file that will be input for metassemble.py
ASSEM1=./masurca/BEEST_RNA_FINALDECONTAM/scaffold/pass1/Scaffolds-pass1.fa # assembly1
ASSEM2=./discovardenovo/Scaffolds_pass6.rename.fa #assembly2
READ1=./reads/mp12kb_R1.fq.gz #read 1, only accepts a single library
READ2=./reads/mp12kb_R2.fq.gz #read 2
READ1_ORIENT=r
READ2_ORIENT=f

ln -sf $ASSEM1 ./assem1.fa
ln -sf $ASSEM2 ./assem2.fa

READ_REPOS=`realpath ./preprocessed_reads/`
mkdir -p $READ_REPOS

#reverse complement if orientation is not r f
if [ "$READ1_ORIENT" == "f" ]; then
	FNAME=`basename $READ1`;
	zcat $READ1 | fastx_reverse_complement -Q33 -z -o $READ_REPOS/$FNAME &
	READ1=$READ_REPOS/$FNAME;
fi

if [ "$READ2_ORIENT" == "r" ]; then
	FNAME=`basename $READ2`;
	zcat $READ2 | fastx_reverse_complement -Q33 -z -o $READ_REPOS/$FNAME &
	READ2=$READ_REPOS/$FNAME;
fi
wait

echo -e "\
############################################\n\
###   Metassemble A and B configuration file\n\
############################################\n\
[global]\n\
\n\
bowtie2_threads=40\n\
bowtie2_read1=$READ1\n\
bowtie2_read2=$READ2\n\
bowtie2_maxins=25000\n\
bowtie2_minins=3000\n\
nucmer_l=100\n\
nucmer_c=300\n\
\n\
\n\
mateAn_A=5000\n\
mateAn_B=20000\n\
\n\
[1]\n\
\n\
fasta=$(pwd)/assem1.fa\n\
ID=Assem1\n\
\n\
[2]\n\
\n\
fasta=$(pwd)/assem2.fa\n\
ID=Assem2\n\
" > S.J.metassemble.config

### Run metassemble
/beegfs/group_dv/software/source/Metassembler/bin/metassemble --conf S.J.metassemble.config --outd ./iter01  

mv S.J.metassemble.config ./iter01
