for i in ./killifishproteins/*; do
	if [ -d $i ]; then
		sp=`basename $i`
		zcat -f $i/*.fa | hhvm preprocessAA.php > $i.filtered.fa 2> $i.filtered.tab
		minreID.py $i.filtered.fa $sp \|
	fi
done

cat *.fst > allkillifishproteins.fa
