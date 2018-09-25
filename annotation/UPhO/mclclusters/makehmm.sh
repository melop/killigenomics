
HMMBUILD=/software/source/hmmer-2.3.2/src/hmmbuild 
folder=UPhO_Seqs/*_clean*.fa
outdir=hmm/
arrFiles=($folder)
THREADS=40

USEDCORE=0

mkdir -p $outdir
for f in "${arrFiles[@]}"
do 
  echo "Processing $f file...";
  FILENAME=`basename $f` ;
  filestem="${FILENAME%_clean*.fa}" ;

  if [ -s $outdir/${filestem}.hmm  ];then
	echo $outdir/${filestem}.hmm already made.
  else
  	$HMMBUILD -F $outdir/${filestem}.hmm $f & 
        let "USEDCORE=USEDCORE+1";
  fi


  if (("$USEDCORE" >= "$THREADS")); then
      let "USEDCORE=0"
      wait;
  fi
  
done

wait
