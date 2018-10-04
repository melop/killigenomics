
arrPops=( ORTWET ORTDRY )
arrChr=( `seq 1 1` )

sChrPrefix=chr


for i in "${arrPops[@]}";do

   for nChr in "${arrChr[@]}";do

		hhvm convertformat.php $i chr${nChr} > ${i}.chr${nChr}.convert.log &
	done
done
wait



