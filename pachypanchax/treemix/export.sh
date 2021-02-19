mkdir -p fortreemix

for i in `seq 0 5`; do
	hhvm export_for_treemix.php $i 24 > fortreemix/$i.log &
done

wait
