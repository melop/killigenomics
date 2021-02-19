mkdir -p fortreemix_PLPonly

for i in `seq 0 23`; do
	hhvm export_for_treemix.php $i 24 > fortreemix_PLPonly/$i.log &
done

wait
