#!/bin/bash
#export PATH=$PATH:`pwd`/../fastahack/:`pwd`/../hmmer-2.3.2/src/:`pwd`/../wise2.2.0/src/bin/
#export WISECONFIGDIR=`pwd`/../wise2.2.0/wisecfg/
GENOMEMASK=$1
GENOME=$2
WORKDIR=$3
QUERY_PROBE=$4
HMM_FILE=$5
THREADS=$6
KEEPMAXLOGRATIO=250

#rm -r $WORKDIR
mkdir -p $WORKDIR

WORKDIR=`realpath $WORKDIR`


if [[ -e $WORKDIR/blast.ok ]]; then
  echo blast ok;
else
	( tblastn -num_threads $THREADS -evalue 1e-10 -db $GENOMEMASK -outfmt 7 -query $QUERY_PROBE > $WORKDIR/blast_ret.txt ) &&  ( touch $WORKDIR/blast.ok ) #use the masked genome for blast
fi ;

php getbesthit.php -i $WORKDIR/blast_ret.txt -o $WORKDIR/coordinates.txt -c $WORKDIR/strands.txt -x 200000 -s 15000 -f $QUERY_PROBE

if [[ -s $WORKDIR/coordinates.txt ]] ; then
echo "coordinate found";
else
echo "blast not found.";
exit;
fi ;

if [[ -e $WORKDIR/genewise.ok ]]; then
  echo genewise ok;
else
 rm -r  $WORKDIR/splitted ;
 fastahack $GENOME -c < $WORKDIR/coordinates.txt > $WORKDIR/fasthackseq.txt


 paste <(sed -e 's/^/\>/' $WORKDIR/coordinates.txt) $WORKDIR/fasthackseq.txt | tr "\t" "\n"    > $WORKDIR/fasthackseq.fas

 mkdir $WORKDIR/splitted; split -l 2 $WORKDIR/fasthackseq.fas $WORKDIR/splitted/seq;
 for f in $WORKDIR/splitted/seq*;  do mv "$f" "$f.fas"; done

 cp genewise3.sh $WORKDIR/splitted/genewise3.sh
 chmod 777 $WORKDIR/splitted/genewise3.sh
 cd $WORKDIR/splitted/
 ( bash genewise3.sh $HMM_FILE $THREADS ) && ( touch $WORKDIR/genewise.ok );

fi;

 cd $WORKDIR/../../

genename=`grep -Po "([A-z_0-9]+)$" $QUERY_PROBE | head -n 1`
( php pickbestmodel2.php -g $genename -f $WORKDIR/splitted/all.ret.gtf -s $WORKDIR/splitted/sum.ret.txt -o $WORKDIR/out.gtf -k 100 -S $KEEPMAXLOGRATIO ) && touch $WORKDIR/out.ok 

