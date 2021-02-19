for i in `seq 1 19`; do
	cat ../../cactus/NFZ2.0/fixgff/cds.gff | awk -v i="$i" '{if ($1=="NFZ.chr"i) print $0}' > chr${i}_cds.gff
	Rscript rphastcon.R $i &
done

wait
