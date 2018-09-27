sPop=$1 #RACDRY.excludehybrid
PREFIX=$3
REP=$2

ANAVAR=/beegfs/group_dv/software/source/anavar_src/anavar

mkdir $sPop

cd $sPop
if [ ! -e ${PREFIX}control.txt ]; then
   cat ../${PREFIX}configheader.txt > ${PREFIX}control.txt
   cat "../../$sPop.SFS.for.anavar.txt" >> ${PREFIX}control.txt
   cat ../${PREFIX}configtail.txt >> ${PREFIX}control.txt
   chmod 555 ${PREFIX}control.txt
fi

#RANDOMSEED

filesize=0
[ -e ${PREFIX}out.rep${REP}.txt ] && filesize=$( stat -c%s "${PREFIX}out.rep${REP}.txt" )


(( filesize < 1000  )) && $ANAVAR ${PREFIX}control.txt ${PREFIX}out.rep${REP}.txt ${PREFIX}log.rep${REP}.txt $(( 333333 + REP ))
