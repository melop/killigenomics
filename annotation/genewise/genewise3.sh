FILES=./*.fas
THREADS=$2

USEDCORE=0
SEQCOUNT=0

getArray() {
    local i=0
    arrStrands=() # Clear array
    while IFS= read -r line # Read a line
    do
        arrStrands+=( "$line" ) # Append line to the array
    done < "$1"
}

getArray "../strands.txt"



for f in $FILES
do
  echo "Processing $f file for hmm $1 ..."
  # take action on each file. $f store current file name
  STRANDCMD="-tfor"
  if [ "${arrStrands[$SEQCOUNT]}" == "+" ]; then
	STRANDCMD="-tfor"; 
	echo "On plus strand"; 
  else
 	STRANDCMD="-trev";
	echo "On reverse strand"; 
  fi 

  (genewise -kbyte 500000 -sum -cdna -gff $STRANDCMD -quiet -hmmer ../../../$1 $f > $f.wise.txt ; \
   csplit --digits=2  --quiet --prefix=$f.wise.split.  $f.wise.txt "/\/\//+1" "{*}"; \
   cat $f.wise.split.00 | tr "//" "\n\n"  >> sum.ret.txt ; \
   cat $f.wise.split.01 | tr "//" "\n\n" >> all.ret.fasta ; \
   cat $f.wise.split.02 | tr "//" "\n\n" >> all.ret.gtf ; \
  ) & 

  let "USEDCORE=USEDCORE+1"
  let "SEQCOUNT=SEQCOUNT+1"
  
  if (("$USEDCORE" >= "$THREADS")); then
      let "USEDCORE=0"
      wait;
  fi
done

wait



