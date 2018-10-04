export PYTHONPATH=
source /beegfs/group_dv/software/source/python_virtualenv/python3.5/bin/activate

SMCPP=smc++
OUTDIR=smcpp_input_rawvcf_varqualmask
sInVCFDIR=`realpath vcfs`
sMaskDIR=`realpath masks`
arrPops=( ORTWET ORTDRY RACWET RACDRY )
arrDistinguishedInd=( ORTWET_528_NO-528-M_D508_N707 ORTDRY_34_NO-34-L_D508_D702 RACWET_96_NR-96-1_S507_D708 RACDRY_43_NR-43-A5_S508_N705 )
arrChr=( `seq 1 19` )

sChrPrefix=chr


mkdir -p $OUTDIR

OUTDIR=`realpath $OUTDIR`;

for i in "${!arrPops[@]}"; do 

	sPop=${arrPops[i]};
	sDistinguishInd=${arrDistinguishedInd[i]};

	sIndList=$sPop:`grep -Po "$sPop\S+" samplenames.txt | tr "\n" ","`
	sIndList=${sIndList%?}

	sDist=`grep -Po "$sDistinguishInd\S*" samplenames.txt | tr "\n" " "`
	sDist=${sDist%?}

	echo Pop $sPop, distinguish ind $sDist, $sIndList

	WD=$OUTDIR/$sPop/; 
	mkdir -p $WD; 


   for nChr in "${arrChr[@]}";do

	sChr=${sChrPrefix}${nChr}
	WDChr=$WD/$sChr
	mkdir -p $WDChr;
	$SMCPP vcf2smc -m $sMaskDIR/${sPop}.mask.varqual.bed.gz -d "$sDist" "$sDist" "$sInVCFDIR/${sPop}.var.vcf.gz" "$WDChr/smc.gz" "$sChr" $sIndList > $WDChr/log.txt 2>&1 & # USE variant quality filter only

	
   done
done
wait

