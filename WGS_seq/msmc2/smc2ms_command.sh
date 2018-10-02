totalchr=240
chrpop1=120
locuslen=1
mu=2.63E-9
Rscript table2step_fun2ms_command.R RAC 258.2e3 $totalchr 1 $chrpop1 $mu $locuslen  > forms_RAC.txt
Rscript table2step_fun2ms_command.R ORT 205.35e3 $totalchr 1 $chrpop1 $mu $locuslen  > forms_ORT.txt

