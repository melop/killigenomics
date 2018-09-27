perl /software/source/PopLD/Package-GFE/Programs/Ext_Ref_Nuc.pl ../ref.fa RefNuc.txt #convert reference fasta

PRODIR=../mapped/pro

TOTALCORES=40
USEDCORES=0
function ProcessPop {
	WD=`pwd`;
	POP=$1
	OUTDIR=$POP
	#OUT=$OUTDIR/In_MIF_$POP.txt
	arrSamples=( `ls $PRODIR/$POP*.pro` );

	mkdir -p $OUTDIR

	ln -s `realpath RefNuc.txt` $OUTDIR/


	for i in "${arrSamples[@]}" 
	do 
		sFileName=`basename $i`;
		sFileStem=${sFileName/%.pro/};
		ln -fs `realpath $i` ./$OUTDIR/

		USEDCORES=$(( USEDCORES+1 ));
		if (( USEDCORES >=  TOTALCORES )); then
			USEDCORES=0;
			wait;
		fi 

		( cd $OUTDIR;\
		  /software/source/PopLD/Package-GFE/Programs/Make_InFile RefNuc.txt 1 $sFileName $sFileStem;\
			cd $WD; ) &
	done;

	wait;

	paste $OUTDIR/RefNuc.txt $OUTDIR/In_GFE_*.txt > $OUTDIR/My_in_GFE.txt;
}

ProcessPop ORTDRY
ProcessPop ORTWET

ProcessPop RACDRY
ProcessPop RACWET


