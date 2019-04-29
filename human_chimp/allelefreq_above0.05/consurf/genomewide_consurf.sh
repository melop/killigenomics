arrPops=( ACB ASW BEB CDX CEU CHB CHS CLM ESN FIN GBR GIH GWD ) 
arrPops=( IBS ITU JPT KHV LWK MSL MXL PEL PJL PUR STU TSI YRI )

for sPop in "${arrPops[@]}"; do
	Rscript genomewide_consurf.2.R $sPop > $sPop.log &
done

wait
