PARTS=13
let END=$PARTS-1;
export OMP_NUM_THREADS=3;
 
for i in $(seq 0 $END); do
 php parallel_fasttree.php  -N $PARTS -f $i -L cleaned_fa_list.txt -D ClusteRs/  > fasttree_log_par${i}_of_${PARTS}.log 2>&1 &
done

wait
