out=hicovreads_err0.001
mkdir -p $out
num_reads=6191754 #72567360*25/(142+151)  25X coverage
read1len=142
read2len=151
fragsize=327
seqerror=0.001
errindelprop=0.001235

for i in pop1_and_pop2_reads.fq_parsed_*/indiv1.fa; do
	sFolder=`dirname $i`;
	sPop=${sFolder/pop1_and_pop2_reads.fq_parsed_/}	
	seed=$RANDOM
	( ./wgsim -N$num_reads -1$read1len -2$read2len -d$fragsize -S$seed -e$seqerror -r0 -R$errindelprop $i ./$out/${sPop}1_read1.fastq ./$out/${sPop}1_read2.fastq ;\
	  gzip ./$out/${sPop}1_read2.fastq ;gzip ./$out/${sPop}1_read2.fastq ) &
done

wait

