#!/usr/bin/bash
#SBATCH -c 40
#SBATCH -p blade


COUNTDUP=/beegfs/group_dv/software/source/countduppe/dist/Release/GNU-Linux/countduppe

arrFiles=(../trim/trimmed/*_*_R1.fq)

NUMTASKS=40
USEDCORE=0

for f in "${arrFiles[@]}"
do 
  
  ( samplename=`basename $f` ;
  samplename=${samplename%'_R1.fq'};
  read2=${f%'_R1.fq'}_R2.fq;


  echo doing $samplename $readcount reads...;

  $COUNTDUP $f $read2 ${samplename}_counts.txt 100 9 0 > ${samplename}_sampledtotalread.txt ; 
  
  preseq lc_extrap -s 3e+06 -o ${samplename}_extrapolate.txt -V ${samplename}_counts.txt ;

   ) &

  let "USEDCORE=USEDCORE+1";

  if (("$USEDCORE" >= "$NUMTASKS")); then
      let "USEDCORE=0"
      wait;
  fi

done




wait

#--log ${samplename}_log.txt


