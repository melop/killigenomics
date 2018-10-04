mkdir lnbams

arrPops=( ORTWET ORTDRY RACWET RACDRY )

for i in "${arrPops[@]}";do
	mkdir lnbams/${i}
	ln -s `realpath ../mapped_hicov/realigned_bam/${i}*.bam` lnbams/${i}/

done

rm lnbams/*/ORTDRY_414_NO-414-S7_D502_D705*
rm lnbams/*/RACWET_92_NR-92-R11_S503_D707*
rm lnbams/*/RACWET_96_NR-96-Y4_S503_D709*
rm lnbams/*/RACDRY_215_NR-215-UP_S507_N701*
