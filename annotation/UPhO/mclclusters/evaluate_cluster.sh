echo "" > evaluate_cluster.out.txt

for i in out.seq.mci.I*; do
	hhvm evaluate_cluster.php $i 16 >> evaluate_cluster.out.txt
done
