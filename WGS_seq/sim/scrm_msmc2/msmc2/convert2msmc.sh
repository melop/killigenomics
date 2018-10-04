export PYTHONPATH=; source /beegfs/group_dv/software/source/python_virtualenv/python3.5/bin/activate

CONVERT=/beegfs/group_dv/software/source/msmc-tools/generate_multihetsep.py
sOUTDIR=msmcinput
mkdir -p $sOUTDIR
for i in vcfs/*.gz;do
	sStem=`basename $i`;
	sStem=${sStem/.vcf.gz/};
	$CONVERT $i > $sOUTDIR/$sStem.multihetsep.txt 2> $sOUTDIR/$sStem.multihetsep.log &
done

wait
