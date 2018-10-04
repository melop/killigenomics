
REVEEL=/beegfs/group_dv/software/source/reveel/release/reveel_caller
BCFTOOLS=/beegfs/group_dv/software/source/samtools-1.6/bin/bin/bcftools
BEAGLE=/beegfs/group_dv/software/source/reveel/beagle.jar
OUTDIR=reveelimpute_2pops
arrPops=( ORT RAC )
arrChr=( `seq 1 10` )
arrChr=( `seq 11 19` )
sChrPrefix=chr
VARQUAL=60

VCFDIR=`realpath vcfs_2pops`

mkdir -p $OUTDIR

OUTDIR=`realpath $OUTDIR`;

for i in "${arrPops[@]}";do

	WD=$OUTDIR/$i/; 
	mkdir -p $WD; 


   for nChr in "${arrChr[@]}";do

	sChr=${sChrPrefix}${nChr}
	WDChr=$WD/$sChr
	mkdir -p $WDChr;

	(	
		cd $WDChr; \
		[ ! -f qry.vcf.gz ] && $BCFTOOLS filter -i 'QUAL>=60' -Oz -r $sChr:1-999999999 $VCFDIR/${i}.var.vcf.gz > ./qry.vcf.gz; \
		[ ! -f qry.vcf.gz.csi ] && $BCFTOOLS index ./qry.vcf.gz; \
		[ ! -f qry.vcf.gz.blk ] && $REVEEL index -b 1000000 qry.vcf.gz; \
		[ ! -f reveelcalled.candidate_num ] && $REVEEL shortlist -n -M 2 -r $sChr:1-999999999 qry.vcf.gz qry.vcf.gz reveelcalled ;\
		[ ! -f reveelcalled.wref.bcf ] && java -Djava.io.tmpdir=$TMPDIR -jar $BEAGLE like=reveelcalled.$sChr.iter10.pr markers=reveelcalled.$sChr.iter10.marker out=prob;\
		[ ! -f reveelcalled.wref.bcf ] && $REVEEL merge qry.vcf.gz reveelcalled.iter10 reveelcalled prob.reveelcalled.$sChr.iter10.pr.dose.gz ;\
		[ ! -f reveelcalled.wref.bcf.csi ] && $BCFTOOLS index reveelcalled.wref.bcf ) >>  "$WDChr/log.txt" 2>&1 &

	done
done
wait



