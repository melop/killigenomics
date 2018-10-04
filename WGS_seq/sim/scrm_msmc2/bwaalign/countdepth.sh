arrSamples=( ORTDRY_34_NO-34-L_D508_D702 ORTWET_528_NO-528-M_D508_N707 RACDRY_43_NR-43-A5_S508_N705 RACWET_96_NR-96-1_S507_D708 )
arrChr=( chr1 chr2 chr3 chr4 chr5 )
outdir=depthcounts
mkdir $outdir
for sSample in "${arrSamples[@]}"; do
	for sChr in "${arrChr[@]}"; do
		/beegfs/group_dv/software/source/samtools-1.6/bin/bin/samtools depth -r $sChr ../mapped_hicov/bam/${sSample}.bam | awk '{sum += $3} END {print sum / NR}' > $outdir/${sSample}_${sChr}.txt &
	done
done

wait
