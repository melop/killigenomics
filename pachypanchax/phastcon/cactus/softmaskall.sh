arrDir=( `ls -d */` )

for i in "${arrDir[@]}"; do
	sSp=`basename $i`;
	hhvm softmask.php $sSp > $sSp/softmask.log &
done

wait
