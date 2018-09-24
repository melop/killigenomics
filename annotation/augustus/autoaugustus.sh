#!/usr/bin/bash
#SBATCH -n 1
#SBATCH -p long
#SBATCH -c 2

source ~/.bashrc
WD=`pwd`
ASSEMBLY=NOR_gapfilled1.fa #the genome assembly to train augustus on
TRAININGGTF=training.gtf #the gene annotations in gtf format
CDNAFASTA=Trinity.cdhitout.filtered.fa #the filtered cds
COPYFROMPROFILE=Jena #copy from this augustus profile to start the training with
AUGUSTUSSPECIESFOLDER=/source/augustus-3.2.1/config/species #root folder where all species profiles are stored.
SPECIESPROFILE=NORgapfilled1 #The new profile name
ALNCDNA=$WD/autoAug/cdna/cdna.psl #aligned cdna file name


#THE FOLLOWING CODE SHOULD ONLY BE RUN ONCE! This copies the profile over

rm -R ${AUGUSTUSSPECIESFOLDER}/${SPECIESPROFILE}
cp -R ${AUGUSTUSSPECIESFOLDER}/${COPYFROMPROFILE}  ${AUGUSTUSSPECIESFOLDER}/${SPECIESPROFILE} 
arrFILES=${AUGUSTUSSPECIESFOLDER}/${SPECIESPROFILE}/${COPYFROMPROFILE}*;
for sF in $arrFILES; do#
	sFBase=$( basename $sF );
	sNewName=$( echo $sFBase | sed "s/^${COPYFROMPROFILE}/${SPECIESPROFILE}/" );
	mv $sF ${AUGUSTUSSPECIESFOLDER}/${SPECIESPROFILE}/${sNewName}
	sed -i "s/${COPYFROMPROFILE}/${SPECIESPROFILE}/" ${AUGUSTUSSPECIESFOLDER}/${SPECIESPROFILE}/${sNewName}
done
rm ${AUGUSTUSSPECIESFOLDER}/${SPECIESPROFILE}/*.orig*

###UNTIL HERE

autoAug.pl -g ${ASSEMBLY} -t ${TRAININGGTF} --species=${SPECIESPROFILE} --cdna=${CDNAFASTA} -v  --useexisting >> autoaug.out.log 2>&1

cd ./autoAug/autoAugPred_abinitio/shells
echo "" > submitcmd.sbatch;
arrFiles=( `find ./ -maxdepth 1  ! -name "*.*" | grep aug` );
if [ -s aug1 ]; then
	for i in "${arrFiles[@]}"; do
		echo '#!/bin/bash' > $i.sh;
		cat $i >> $i.sh;
		chmod 777 $i.sh;
		echo $i.sh >> submitcmd.sbatch;
	done
fi

cat submitcmd.sbatch | ./submit -m 5G -c 1; 

cd $WD

autoAug.pl -g ${ASSEMBLY} -t ${TRAININGGTF} --species=${SPECIESPROFILE} --cdna=${CDNAFASTA} -v  --useexisting >> autoaug.out.log 2>&1 #

autoAug.pl --species=${SPECIESPROFILE} --genome=${WD}/autoAug/seq/genome_clean.fa --useexisting --hints=${WD}/autoAug/hints/hints.E.gff  -v -v -v  --index=1 >> autoaug.out.log 2>&1


cd ./autoAug/autoAugPred_hints/shells/
echo "" > submitcmd.sbatch;
arrFiles=( `find ./ -maxdepth 1  ! -name "*.*" | grep aug` );
if [ -s ./aug1 ]; then
	for i in "${arrFiles[@]}"; do
		echo '#!/bin/bash' > $i.sh;
		cat $i >> $i.sh;
		chmod 777 $i.sh;
		echo $i.sh >> ./submitcmd.sbatch;
	done
fi


cat submitcmd.sbatch | ../../../submit -m 5G -c 1; 

cd $WD

autoAug.pl --species=${SPECIESPROFILE} --genome=${WD}/${ASSEMBLY} --useexisting --hints=${WD}/autoAug/hints/hints.E.gff --estali=${WD}/autoAug/cdna/cdna.psl -v -v -v  --index=2  >> autoaug.out.log 2>&1


cd ./autoAug/autoAugPred_hints_utr/shells
echo "" > submitcmd.sbatch;
arrFiles=( `find ./ -maxdepth 1  ! -name "*.*" | grep aug` );
if [ -s ./aug1 ]; then
	for i in "${arrFiles[@]}"; do
		echo '#!/bin/bash' > $i.sh;
		cat $i >> $i.sh;
		chmod 777 $i.sh;
		echo $i.sh >> ./submitcmd.sbatch;
	done
fi


cat submitcmd.sbatch | ../../../submit -m 5G -c 1; 

cd $WD

autoAug.pl --species=${SPECIESPROFILE} --genome=${WD}/${ASSEMBLY} --useexisting --hints=${WD}/autoAug/hints/hints.E.gff  -v -v -v  --index=3  >> autoaug.out.log 2>&1






