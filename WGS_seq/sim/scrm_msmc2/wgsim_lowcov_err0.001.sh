out=lowcovreads_err0.001
mkdir -p $out
read1len=142
read2len=151
fragsize=327
seqerror=0.001
errindelprop=0.001235
locuslen=72567360
totalreadlen=$(( read1len + read2len ));

for i in pop1_and_pop2_reads.fq_parsed_*/indiv*.fa; do
	sleep 5;
	sFolder=`dirname $i`;
	sFAName=`basename $i`
	sFANameStem=${sFAName/.fa/}	
	sPop=${sFolder/pop1_and_pop2_reads.fq_parsed_/}	
	sCov=`grep -A1 $sPop coverages.txt | tail -n 1`
	nNumReads=`Rscript getreads_locuslen.R $sCov $locuslen $totalreadlen`
	seed=$RANDOM
	echo wgsim -N$nNumReads -1$read1len -2$read2len -d$fragsize -S$seed -e$seqerror -r0 -R$errindelprop $i ./$out/${sPop}_${sFANameStem}_read1.fastq ./$out/${sPop}_${sFANameStem}_read2.fastq
	( ./wgsim -N$nNumReads -1$read1len -2$read2len -d$fragsize -S$seed -e$seqerror -r0 -R$errindelprop $i ./$out/${sPop}_${sFANameStem}_read1.fastq ./$out/${sPop}_${sFANameStem}_read2.fastq ;\
	  gzip ./$out/${sPop}_${sFANameStem}_read1.fastq ;gzip ./$out/${sPop}_${sFANameStem}_read2.fastq) &

done

wait

