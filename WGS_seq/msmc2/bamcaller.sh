SAMTOOLS=/software/source/samtools-1.6/bin/bin/samtools
BCFTOOLS=/software/source/bcftools-1.6/bcftools
BAMCALLER=/software/source/msmc-tools/bamCaller.py

export PYTHONPATH=; source /software/source/python_virtualenv/python3.5/bin/activate

CONVERT=/software/source/msmc-tools/generate_multihetsep.py



arrSamples=( ORTDRY_34_NO-34-L_D508_D702 ORTWET_528_NO-528-M_D508_N707 RACDRY_43_NR-43-A5_S508_N705 RACWET_96_NR-96-1_S507_D708 )
arrCov=( 32.71566 26.10244 40.31834 26.4357 )

arrChr=( `seq 1 10` )
arrChr=( `seq 11 19` )

sOutDir=formsmc2_in/


for nSample in "${!arrSamples[@]}"; do 
	sSample=${arrSamples[$nSample]};
	nCov=${arrCov[$nSample]};

	for nChr in "${arrChr[@]}"; do 
		sChrDir=$sOutDir/$sSample/chr$nChr
		mkdir -p $sChrDir
		( if [ ! -s $sChrDir/out_mask.bed.gz ];then $SAMTOOLS mpileup -q 20 -Q 20 -C 50 -u -r chr$nChr -f ../ref.fa ../mapped_hicov/bam/$sSample.bam | $BCFTOOLS call -c -V indels | $BAMCALLER $nCov $sChrDir/out_mask.bed.gz | gzip -c > $sChrDir/out.vcf.gz 2> $sChrDir/log.txt; fi; \
		if [ ! -e $sChrDir/done.txt ]; then $CONVERT --mask $sChrDir/out_mask.bed.gz $sChrDir/out.vcf.gz > $sChrDir/formsmc2.multihetsep.txt 2> $sChrDir/formsmc2.multihetsep.log && touch $sChrDir/done.txt; fi; ) &
	done

done

wait
