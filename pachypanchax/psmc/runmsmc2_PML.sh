MSMC2=/beegfs/group_dv/software/source/msmc2/build/release/msmc2_linux64bit

arrPops=(  PML )
sIn=formsmc2_in
sOut=msmc2ret
mkdir -p $sOut


for sPop in "${arrPops[@]}"; do
	sFiles="${sIn}/${sPop}*/chr*/formsmc2.multihetsep.txt"
	sCmd="$MSMC2 -o $sOut/$sPop -i 40 -t 2 -r 2.61976354 $sFiles > /dev/null "
	eval $sCmd &
done

wait
