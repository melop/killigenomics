SAMTOOLS=/beegfs/group_dv/software/source/samtools-1.6/bin/bin/samtools
BCFTOOLS=/beegfs/group_dv/software/source/bcftools-1.6/bcftools
BAMCALLER=/beegfs/group_dv/software/source/msmc-tools/bamCaller.py

export PYTHONPATH=; source /beegfs/group_dv/software/source/python_virtualenv/python3.5/bin/activate

CONVERT=/beegfs/group_dv/software/source/msmc-tools/generate_multihetsep.py



arrSamples=( ORTDRY ORTWET RACDRY RACWET )
arrCov=( 25 25 25 25 )

arrChr=( All )

sOutDir=formsmc2_in/


for nSample in "${!arrSamples[@]}"; do 
	sSample=${arrSamples[$nSample]};
	nCov=${arrCov[$nSample]};

	for nChr in "${arrChr[@]}"; do 
		sChrDir=$sOutDir/$sSample/chr$nChr
		mkdir -p $sChrDir
		( if [ ! -s $sChrDir/out_mask.bed.gz ];then $SAMTOOLS mpileup -q 20 -Q 20 -C 50 -u -r chr$nChr -f chr1-19.fa mapped_hicov_gatk_realigned/bam/${sSample}*.bam | $BCFTOOLS call -c -V indels | $BAMCALLER $nCov $sChrDir/out_mask.bed.gz | gzip -c > $sChrDir/out.vcf.gz 2> $sChrDir/log.txt; fi; \
		if [ ! -e $sChrDir/done.txt ]; then $CONVERT --mask $sChrDir/out_mask.bed.gz $sChrDir/out.vcf.gz > $sChrDir/formsmc2.multihetsep.txt 2> $sChrDir/formsmc2.multihetsep.log && touch $sChrDir/done.txt; fi; ) &
	done

done

wait
