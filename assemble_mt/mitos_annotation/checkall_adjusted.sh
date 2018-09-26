for i in mt_pos_adjusted/*.fasta; do
	echo  $i ...
	sFileName=`basename $i`;
	hhvm checkresults.php $i out_adjusted > out_adjusted/$sFileName/checklog.txt &
done
wait
