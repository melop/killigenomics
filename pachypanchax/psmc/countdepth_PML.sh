arrFiles=( `ls ../mapped_PML/realigned_bam/*.realigned.bam` )
arrChr=( chrLG1 chrLG2 chrLG3 chrLG4 chrLG5 )
outdir=depthcounts
mkdir $outdir
for sBAM in "${arrFiles[@]}"; do
	sSample=`basename $sBAM`;
	sSample=${sSample/.realigned.bam/};
	for sChr in "${arrChr[@]}"; do
		/beegfs/group_dv/software/source/samtools-1.6/bin/bin/samtools depth -r $sChr $sBAM | awk '{sum += $3} END {print sum / NR}' > $outdir/${sSample}_${sChr}.txt &
	done
done

wait
