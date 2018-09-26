for i in mt/*.fasta; do
	echo adjust $i ...
	hhvm avoid_breaking_genes.php $i &
done
wait
