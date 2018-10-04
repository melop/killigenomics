#/beegfs/group_dv/software/source/LDhat/rhomap -seq sites.txt -loc loc.txt -lk lk_n114_t0.001 -bpen 20  -samp 10000 -its 3000000 -seed 33333 -burn 300000 -prefix rhomap > log.txt 2>&1 &
LIKELIHOODTABLE=$1

ln -sf `realpath $LIKELIHOODTABLE` ./lktable.txt
nTestFlag=0
[ -e log.txt* ] && nTestFlag=`zcat -f log.txt* | grep "Final block length"  | wc -l`; #check if done already
if (( nTestFlag == 1 )); then
	echo interval finished already;
else
	/beegfs/group_dv/software/source/LDhat/interval -seq sites.txt -loc loc.txt -lk lktable.txt -its 4000000 -burn 30000 -samp 2000 -bpen 50  > log.txt 2>&1 
fi

if [ ! -e rates.txt ]; then
	if [ ! -e rates.txt.gz ]; then 
		realpath `pwd`
		echo Error occured interval program;
		exit
	fi
fi

if [ -e rates.txt.gz ]; then  
	mkfifo rates.txt;
	zcat  rates.txt.gz > rates.txt &
fi

/beegfs/group_dv/software/source/LDhat/stat -input rates.txt -burn 1000 -loc loc.txt

[ -e type_table.txt ] && rm type_table.txt
[ -e new_lk.txt ] && rm new_lk.txt

if [ ! -e bounds.txt.gz ]; then  
	gzip bounds.txt
fi
if [ ! -e rates.txt.gz ]; then  
	gzip rates.txt
fi
if [ ! -e log.txt.gz ]; then  
gzip log.txt
fi



