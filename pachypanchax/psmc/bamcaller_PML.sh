SAMTOOLS=/beegfs_old/group_dv/software/source/samtools-1.6/bin/bin/samtools
BCFTOOLS=/beegfs_old/group_dv/software/source/bcftools-1.6/bcftools
BAMCALLER=/beegfs_old/group_dv/software/source/msmc-tools/bamCaller.py

#export PYTHONPATH=; source /beegfs_old/group_dv/software/source/python_virtualenv/python3.5/bin/activate

CONVERT=/beegfs_old/group_dv/software/source/msmc-tools/generate_multihetsep.py

sBAMFolder=../mapped_PML/realigned_bam/

arrSamples=( PML )
arrCov=( 15.1963 )

arrChr=( `seq 1 24` )

sOutDir=formsmc2_in/


for nSample in "${!arrSamples[@]}"; do 
	sSample=${arrSamples[$nSample]};
	nCov=${arrCov[$nSample]};

	for nChr in "${arrChr[@]}"; do 
		sChrDir=$sOutDir/$sSample/chrLG$nChr
		mkdir -p $sChrDir
		( if [ ! -s $sChrDir/out_mask.bed.gz ];then $SAMTOOLS mpileup -q 60 -Q 20 -C 0 -u -r chrLG$nChr -f ../ref.fa $sBAMFolder/$sSample.realigned.bam | $BCFTOOLS call -c -V indels | $BAMCALLER $nCov $sChrDir/out_mask.bed.gz | gzip -c > $sChrDir/out.vcf.gz 2> $sChrDir/log.txt; fi; \
		if [ ! -e $sChrDir/done.txt ]; then $CONVERT --mask $sChrDir/out_mask.bed.gz $sChrDir/out.vcf.gz > $sChrDir/formsmc2.multihetsep.txt 2> $sChrDir/formsmc2.multihetsep.log && touch $sChrDir/done.txt; fi; ) &
	done

done

wait
