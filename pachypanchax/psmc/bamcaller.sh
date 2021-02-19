SAMTOOLS=/beegfs_old/group_dv/software/source/samtools-1.6/bin/bin/samtools
BCFTOOLS=/beegfs_old/group_dv/software/source/bcftools-1.6/bcftools
BAMCALLER=/beegfs_old/group_dv/software/source/msmc-tools/bamCaller.py

#export PYTHONPATH=; source /beegfs_old/group_dv/software/source/python_virtualenv/python3.5/bin/activate

CONVERT=/beegfs_old/group_dv/software/source/msmc-tools/generate_multihetsep.py

sBAMFolder=../mapped/realigned_bam/

arrSamples=( PLP_A_PLP-A3_D501_D701
PLP_A_PLP-A4_D502_D702
PLP_B_PLP-B3_D504_D704
PLP_C_PLP-C1_D506_D706
PLP_C_PLP-C2_D507_D707
PLP_C_PLP-C3_D508_D708
PLP_C_PLP-C4_D502_D701
PLP_D_PLP-D0_ref
PLP_D_PLP-D1_D503_D702
PLP_E_PLP-E1_D505_D704
PLP_E_PLP-E2_D506_D705
PLP_E_PLP-E3_D507_D706
PLP_E_PLP-E4_D508_D707
PLP_E_PLP-E5_D501_D708
PLP_F_PLP-F1_D503_D701
PLP_F_PLP-F2_D504_D702
PLP_F_PLP-F3_D505_D703
PLP_G_PLP-G1_D506_D704
PLP_G_PLP-G3_D508_D706
PLP_G_PLP-G5_D502_D708
 )
arrCov=( 17.5716
14.8559
20.0662
19.7482
15.2611
12.8079
15.1832
76.1434
18.8164
22.3941
19.032
17.0744
18.546
15.6852
14.1833
18.6097
18.8801
19.515
18.1725
18.5816
 )

arrChr=( `seq 1 2` )
arrChr=( `seq 3 4` )
arrChr=( `seq 5 8` )
arrChr=( `seq 9 11` )
arrChr=( `seq 12 14` )
arrChr=( `seq 15 24` )
#arrChr=( `seq 11 19` )

sOutDir=formsmc2_in/


for nSample in "${!arrSamples[@]}"; do 
	sSample=${arrSamples[$nSample]};
	nCov=${arrCov[$nSample]};

	for nChr in "${arrChr[@]}"; do 
		sChrDir=$sOutDir/$sSample/chrLG$nChr
		mkdir -p $sChrDir
		( if [ ! -s $sChrDir/out_mask.bed.gz ];then $SAMTOOLS mpileup -q 20 -Q 20 -C 50 -u -r chrLG$nChr -f ../ref.fa $sBAMFolder/$sSample.realigned.bam | $BCFTOOLS call -c -V indels | $BAMCALLER $nCov $sChrDir/out_mask.bed.gz | gzip -c > $sChrDir/out.vcf.gz 2> $sChrDir/log.txt; fi; \
		if [ ! -e $sChrDir/done.txt ]; then $CONVERT --mask $sChrDir/out_mask.bed.gz $sChrDir/out.vcf.gz > $sChrDir/formsmc2.multihetsep.txt 2> $sChrDir/formsmc2.multihetsep.log && touch $sChrDir/done.txt; fi; ) &
	done

done

wait
