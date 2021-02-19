FILE=../../fortreemix.txt.gz
treemix=/beegfs/group_dv/software/source/treemix/src/treemix
nCPU=40

nCount=0;

for m in {1..5}
   do
   for i in {1..20}
      do
      nCount=$((nCount+1));
      if [ $nCount -gt $nCPU ]; then
        wait;
        nCount=1;
      fi
      k=$((i*200));
      if [ -s test.${i}.${m}.llik ]; then
        echo test.${i}.${m}.llik done
      else
          $treemix -i $FILE -k $k -global -root PML -m $m -o test.${i}.${m} > log.${i}.${m}.log 2>&1 &
      fi
      done 
done

m=0
for i in {1..20}
do
      nCount=$((nCount+1));
      if [ $nCount -gt $nCPU ]; then
        wait;
        nCount=1;
      fi
      
      k=$((i*200));
      if [ -s test.${i}.${m}.llik ]; then
        echo test.${i}.${m}.llik done
      else
          $treemix -i $FILE -k $k -global -root PML -o test.${i}.${m} > log.${i}.${m}.log 2>&1 &
      fi
done

wait
