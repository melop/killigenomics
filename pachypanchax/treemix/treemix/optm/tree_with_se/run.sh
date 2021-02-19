FILE=../../../fortreemix.txt.gz
treemix=/beegfs/group_dv/software/source/treemix/src/treemix
m=1
i=2
k=$((i*200))
$treemix -i $FILE -k $k -global -root PML -m $m -o test.${i}.${m} -se > log.${i}.${m}.log 2>&1
