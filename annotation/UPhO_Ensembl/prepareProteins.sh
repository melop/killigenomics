for i in ../ensembl/*; do
	if [ -d $i ]; then
		sp=`basename $i`
		zcat $i/*.pep.all.fa.gz | hhvm preprocessAA.php > $i.filtered.fa 2> $i.filtered.tab
		minreID.py $i.filtered.fa $sp \|
	fi
done

cat *.fst > allproteins.fa
