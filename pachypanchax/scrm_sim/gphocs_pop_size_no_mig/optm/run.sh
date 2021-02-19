FILE=../fortreemix.txt.gz
treemix=/beegfs/group_dv/software/source/treemix/src/treemix


for m in {1..5}
   do
   for i in {1..20}
      do
      k=$((i*200));
$treemix -i $FILE -k $k -global -m $m -o test.${i}.${m} > log.${i}.${m}.log 2>&1 &

      done 
wait
done

m=0
for i in {1..20}
do
      k=$((i*200));
$treemix -i $FILE -k $k -global -o test.${i}.${m} > log.${i}.${m}.log 2>&1 &

done

wait
