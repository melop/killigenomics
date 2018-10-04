MSMC2=/beegfs/group_dv/software/source/msmc2/build/release/msmc2_linux64bit

sOut=msmc2ret
mkdir -p $sOut

for i in msmcinput/*.multihetsep.txt; do
	sStem=`basename $i`
	sStem=${sStem/.multihetsep.txt/};
	sCmd="$MSMC2 -o $sOut/$sStem -i 40 -t 10 -r 2.61976354 $i"
	eval $sCmd &
done

wait
