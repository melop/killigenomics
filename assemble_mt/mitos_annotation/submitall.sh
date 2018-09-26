for i in mt/*.fasta; do
	echo submit $i ...
	hhvm submit_annotation_jobs.php $i;
	sleep 2;
done
