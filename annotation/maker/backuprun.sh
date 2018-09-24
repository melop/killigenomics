
backupname=genome.maker.output.round1.snap

if [ -e $backupname ]; then
	echo The backup name $backupname is already existing. error;
	exit -1;
fi

mv genome.maker.output $backupname
cp -R $backupname genome.maker.output
