export PYTHONPATH=
source /beegfs/group_dv/software/source/python_virtualenv/python3.5/bin/activate

SMCPP=smc++

OUTDIR=smcpp_2pops_input_reveel
sInVCFDIR=`realpath reveelimpute_2pops`
sMaskDIR=`realpath masks_2pops`
sSp=ORT
arrPops=( ${sSp}WET ${sSp}DRY )
arrDistinguishedInd=( ORTWET_528_NO-528-M_D508_N707 ORTDRY_34_NO-34-L_D508_D702 )
arrChr=( `seq 1 19` )

sChrPrefix=chr


mkdir -p $OUTDIR

OUTDIR=`realpath $OUTDIR`;



	sPop1=${arrPops[0]};
	sDistinguishInd1=${arrDistinguishedInd[0]};

	sIndList1=$sPop1:`grep -Po "$sPop1\S+" samplenames.txt | tr "\n" ","`
	sIndList1=${sIndList1%?}

	sDist1=`grep -Po "$sDistinguishInd1\S*" samplenames.txt | tr "\n" " "`
	sDist1=${sDist1%?}

	sPop2=${arrPops[1]};
	sDistinguishInd2=${arrDistinguishedInd[1]};

	sIndList2=$sPop2:`grep -Po "$sPop2\S+" samplenames.txt | tr "\n" ","`
	sIndList2=${sIndList2%?}

	sDist2=`grep -Po "$sDistinguishInd2\S*" samplenames.txt | tr "\n" " "`
	sDist2=${sDist2%?}


	echo Pop $sPop1 $sPop2

	WD=$OUTDIR/$sSp/; 
	mkdir -p $WD; 


   for nChr in "${arrChr[@]}";do

	sChr=${sChrPrefix}${nChr}
	WDChr=$WD/$sChr
	mkdir -p $WDChr;
	$SMCPP vcf2smc -m $sMaskDIR/${sSp}.mask.varqual.bed.gz -d "$sDist1" "$sDist1" "$sInVCFDIR/${sSp}/${sChr}/reveelcalled.wref.bcf" "$WDChr/smc.12.gz" "$sChr" $sIndList1 $sIndList2 > $WDChr/log.12.txt 2>&1 &

	$SMCPP vcf2smc -m $sMaskDIR/${sSp}.mask.varqual.bed.gz -d "$sDist2" "$sDist2" "$sInVCFDIR/${sSp}/${sChr}/reveelcalled.wref.bcf" "$WDChr/smc.21.gz" "$sChr" $sIndList2 $sIndList1 > $WDChr/log.21.txt 2>&1 &

	
   done

wait

