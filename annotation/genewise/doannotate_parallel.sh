#!/usr/bin/bash
TOTALPARTS=$1
THISPART=$2
THREADSPERTASK=1

echo part $THISPART of $TOTALPARTS

#export PATH=/beegfs/group_dv/software/source/wise2.4.1/ubuntubin/:$PATH
WDIR=`pwd`
outfolder=out/
folder=hmm/*.hmm
arrFiles=($folder)
GENOME=/beegfs/group_dv/home/RCui/killifish_genomes/denovo/release/PLP/1.0/PLPv1.0.fa

#ln -s `realpath $GENOME` genome.fa
#make blast db, only need it once:

#makeblastdb -in genome.fa -dbtype nucl -parse_seqids -out genome.fa

mkdir -p $outfolder


JOBCOUNT=0


for f in "${arrFiles[@]}"
do 
  nMODULO=$(( $JOBCOUNT % $TOTALPARTS ));
  #echo modulo $nMODULO;

  let JOBCOUNT=JOBCOUNT+1;
  if (( $nMODULO != $THISPART )) ; then
	continue;
  fi 

  echo "Processing $JOBCOUNT  $f file...";

#continue;
  cd $WDIR;
  FILENAME=`basename $f` ;
  GENE="${FILENAME%.*}" ;

	if [[ -e $outfolder/par_${GENE}/out.ok ]] ; then
	echo "already exists exit.";
	continue;
	fi ;

  ./pullout3.sh genome.masked.fa genome.fa $outfolder/par_${GENE}/ prot_nogap/${GENE}.fasta hmm/${GENE}.hmm $THREADSPERTASK 2>/dev/null  ;



  
done





