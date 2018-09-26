for i in mt/*.fasta; do
	echo  $i ...
	sFileName=`basename $i`;
	hhvm checkresults.php $i  > out/$sFileName/checklog.txt &
done
wait
