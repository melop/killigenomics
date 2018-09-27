sPop=$1 #RACDRY.excludehybrid
PREFIX=$3
REP=$2

ANAVAR=/beegfs/group_dv/software/source/anavar_src/anavar

mkdir $sPop

cd $sPop
sControlFile=${PREFIX}.bsrep${REP}.control.txt

cat ../${PREFIX}configheader.txt > $sControlFile
cat "../bs/$sPop/rep${REP}.txt" >> $sControlFile
cat ../${PREFIX}configtail.txt >> $sControlFile


#RANDOMSEED

sOutFile=${PREFIX}out.bsrep${REP}.txt
sLogFile=${PREFIX}out.bsrep${REP}.log
filesize=0
[ -e $sOutFile ] && filesize=$( stat -c%s "$sOutFile" )


(( filesize < 1000  )) && $ANAVAR $sControlFile $sOutFile $sLogFile $(( 333333 + REP ))
