totalchr=240
chrpop1=120
locuslen=10000000
mu=2.63E-9
Rscript table2step_fun2ms_command.R demography_RAC_RACWET_table_for_smc2ms.txt demography_RAC_RACDRY_table_for_smc2ms.txt $totalchr 1 $chrpop1 $mu $locuslen 258202.574154806 > forms_RAC.txt
Rscript table2step_fun2ms_command.R demography_ORT_ORTWET_table_for_smc2ms.txt demography_ORT_ORTDRY_table_for_smc2ms.txt $totalchr 1 $chrpop1 $mu $locuslen 205352.6086637803 > forms_ORT.txt

