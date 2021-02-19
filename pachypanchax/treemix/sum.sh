s=""
for i in `seq 1 24`; do
	s="$s fortreemix/in_chrLG$i.txt "
done

cat $s | gzip -c > fortreemix.txt.gz
