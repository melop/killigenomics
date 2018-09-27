sPOP=RACDRY.excludehybrid
sPrefix=
sPrefix=sametheta_
sPrefix=sametheta_nopolerr_

for i in ${sPOP}/${sPrefix}out.rep*.txt;do 
	echo $i

	if [ "$i" == "${sPOP}/${sPrefix}out.rep0.txt" ]; then
		cp $i ${sPOP}/sum_${sPrefix}out.txt; 
#		echo "" >> ${sPOP}/sum_${sPrefix}out.txt; 
	else	
		tail -n1 $i >> ${sPOP}/sum_${sPrefix}out.txt; 
#		echo "" >> ${sPOP}/sum_${sPrefix}out.txt; 
	fi
done
