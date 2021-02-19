for i in `seq 0 7`; do
	hhvm coverage_mask.php 8 $i > cov.$i.log &
done

wait
