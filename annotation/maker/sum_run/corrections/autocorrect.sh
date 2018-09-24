#!/bin/bash
#SBATCH -p blade
#SBATCH -n 1
#SBATCH -c 2

module load hhvm
module load R



sGenome=../../genome/genome_pilon.fa #Genome Assembly
sMaskedGenome=../../repeatmasker/NOR_v1.0_pilon/scf.fa.masked #Repeat masked assembly
sGFF=../genome_pilon.all.gff #The GFF output from MAKER

sCopyFrom=../run_pilon_exonerateOn/corrections #copy scripts from this folder.

sJobID="";
nMaxTry=1000;

function SubmitJob {
	sJobID=`sbatch $1`
	sJobID=`echo $sJobID | grep -oP "\d+"`
	Wait $sJobID
}

function Wait {
	nJobID=$1
	while true; do
		sCode=`squeue -j $nJobID 2>/dev/null | wc -l`
		if [[ "$sCode" != "0" ]]; then
			#not done yet
			sleep 60;
		else
			break;
		fi
	done
}


echo copying and linking files...

 cp  $sCopyFrom/*.R ./
 cp $sCopyFrom/*.php ./
 cp $sCopyFrom/*.sbatch ./
 cp $sCopyFrom/*.sh ./
 ln -sf $sCopyFrom/all.protein.evidence.fa* ./ #This is all protein evidences used in the maker annotation , cat them together into one file
 ln -sf `realpath $sGenome` ./query.fa
 ln -sf `realpath $sMaskedGenome` ./query.masked.fa
 ln -sf $sGFF ./genome.all.gff

#exit;

echo building fastahack index...
fastahack -i query.fa
fastahack -i query.masked.fa

echo spliting GFF ...
hhvm splitGff.php genome.all.gff

mkdir -p logs

if [ -s chimericfixed.gff ]; then
	echo chimaeric splitting already finished;
else
	echo splitting chimaeric models...
	SubmitJob splitchimeric.sbatch
	echo summarizing chimaeric models...
	hhvm splitchimeric_sum.php > splitchimeric.log 2>&1

fi


if [ -s maker.refined.gff ]; then
	echo refine model finished...
else

	for nTry in $(seq 1 $nMaxTry); do 

		echo refining models...
		SubmitJob refinemodels.sbatch
		echo checking leftover ones...

		echo summarizing refined models...
		hhvm refine_sum.php > refinemodel.log 2>&1
		sErr=`grep "FINISHED!" refinemodel.log | wc -l`
		if [[ $sErr == "1" ]]; then
			echo finished.	
			break;
		else 
			if [[ $nTry ==  $nMaxTry ]]; then
				echo refine models failed.
				exit 2;
			else
				echo Give another try
				continue;
			fi
		fi

	done


fi

echo keep best isoform...
hhvm keeplongestisoform.php > removeisoform.log 2>&1

echo add missed genewise models back...
hhvm add_nooverlap_candidates.php > addgenewiseback.log 2>&1;

if [ -s maker.replace_frag_aln.gff ]; then
	echo replacefrag finished
else
	echo merge fragmented genes with genewise and update annotations...
	SubmitJob replace_frag_algn.sbatch
	hhvm replace_frag_aln_sum.php > frag2aln.log 2>&1
fi

echo generating fasta files...
bash maker.refine.proteins.sh









