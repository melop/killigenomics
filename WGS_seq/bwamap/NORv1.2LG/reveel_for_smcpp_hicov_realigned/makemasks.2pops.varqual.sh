mkdir -p masks_2pops
arrPops=( ORT RAC )


for i in "${!arrPops[@]}"; do 
	sPop=${arrPops[$i]};
	#call variants using the FTH outgroup as the reference.
	hhvm makemask_2pops_varqual.php $sPop  > masks_2pops/$sPop.varqual.log 2>&1 &

done
wait
