mkdir -p masks
arrPops=( ORTWET ORTDRY RACWET RACDRY )


for i in "${!arrPops[@]}"; do 
	sPop=${arrPops[$i]};
	#call variants using the FTH outgroup as the reference.
	hhvm makemask_varqual.php $sPop  > masks/$sPop.varqual.log 2>&1 &

done
wait
