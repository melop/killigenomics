
num_reads=83586399 #900M*25/(142+151)  25X coverage
read1len=142
read2len=151
fragsize=327
seqerror=0.0033
errindelprop=0.001235

out=hicovreads_err$seqerror

arrPops=( ORTDRY ORTWET RACDRY RACWET )

sFADir=ind_fa

mkdir -p $sFADir
mkdir -p $out

for sPop in "${arrPops[@]}";do
	sSp=${sPop:0:3}
	sSimFile=sim${sSp}.out.txt.gz;	
	[ ! -e $sFADir/$sPop.1.fa.done ] && zcat -f $sSimFile | grep -m1 -A1 "$sPop\\.hap1$"  > $sFADir/$sPop.1.fa && \
        zcat -f $sSimFile | grep -m1 -A1 "$sPop\\.hap2$"  >> $sFADir/$sPop.1.fa && fold -w100 $sFADir/$sPop.1.fa > $sFADir/$sPop.1.fa.fold && mv $sFADir/$sPop.1.fa.fold  $sFADir/$sPop.1.fa && touch $sFADir/$sPop.1.fa.done
	seed=$RANDOM
	sCMD="./wgsim -N$num_reads -1$read1len -2$read2len -d$fragsize -S$seed -e$seqerror -r0 -R$errindelprop \"$sFADir/$sPop.1.fa\" ./$out/${sPop}1_read1.fastq ./$out/${sPop}1_read2.fastq";
	echo $sCMD;
	( eval $sCMD ;\
	  gzip ./$out/${sPop}1_read1.fastq ;gzip ./$out/${sPop}1_read2.fastq ) &
done

wait

