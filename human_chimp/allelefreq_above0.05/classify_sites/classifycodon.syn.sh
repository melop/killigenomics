arrPops=( ACB ASW BEB CDX CEU CHB CHS CLM ESN FIN GBR GIH GWD )

#arrPops=( IBS ITU JPT KHV LWK MSL MXL PEL PJL PUR STU TSI YRI ) 
for pop in "${arrPops[@]}"; do

hhvm classifycodon.syn.php -p $pop > $pop.syn.classify.log 2>&1 &

done
wait
