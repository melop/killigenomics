MSMC2=/software/source/msmc2/build/release/msmc2_linux64bit

nReps=30
nChunkSize=5000000


arrPops=( ORTDRY ORTWET RACDRY RACWET )

sOutDir=bootstrapped


for sPop in "${arrPops[@]}"; do
	for i in $(seq 1 $nReps); do 

		sBSFolder="$sOutDir/$sPop/_${i}";
		[ ! -s $sBSFolder/out.final.txt ] && $MSMC2 -o $sBSFolder/out -i 40 -t 1 -r 2.61976354 $sBSFolder/bootstrap_multihetsep.*.txt 2> /dev/null &
	done

done

wait
