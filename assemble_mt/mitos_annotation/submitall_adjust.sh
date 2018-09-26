for i in mt_pos_adjusted/*.fasta; do
	echo submit $i ...
	hhvm submit_annotation_jobs.php $i out_adjusted;
	sleep 2;
done
