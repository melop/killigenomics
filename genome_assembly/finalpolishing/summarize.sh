sIn="./release/NOR/1.0/NORv1.0.fa"; #Assembly to be fixed
sMasked="./release/NOR/1.0/NORv1.0.hardmasked.fa"; #The masked version of the assembly
AGP2FASTA="./agp2fasta.pl" #the script to convert AGP+Fasta to a new Fasta
sTmpFolder="./tmp_aggresive/"; #output folder of blastmerge_cov.sh
sSumOut=summed.blastmerged.fa #output file
sSumAGP=summed.blastmerged.agp #output file
sMaskedOut=summed.blastmerged.hardmasked.fa #output masked version of the lifted genome

rm $sSumOut;
rm $sSumAGP;

arrScf=( `grep ">" $sIn | sed "s/>//"` )

for sScf in "${arrScf[@]}"
do
  sWD=$sTmpFolder/$sScf;
  sAGP=$sWD/overlap.joined.agp;
  sScf=$sWD/blast.gapmerged.fa;

  if [ ! -e $sAGP ]; then
	echo Error $sAGP not found!
	exit 2;
  fi

  if [ ! -e $sScf ]; then
	echo Error $sScf not found!
	exit 3;
  fi

  cat $sScf >> $sSumOut;
  cat $sAGP >> $sSumAGP;

done

perl $AGP2FASTA $sSumAGP $sMasked ;
