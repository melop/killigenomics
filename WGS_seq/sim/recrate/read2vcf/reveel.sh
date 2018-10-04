
REVEEL=/beegfs/group_dv/software/source/reveel/release/reveel_caller
BCFTOOLS=/beegfs/group_dv/software/source/samtools-1.6/bin/bin/bcftools
BEAGLE=/beegfs/group_dv/software/source/reveel/beagle.jar
OUTDIR=reveelimpute
arrPops=( ORTWET ORTDRY )
arrChr=( `seq 1 1` )

sChrPrefix=chr

VCFDIR=`realpath vcfs`

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
		$BCFTOOLS view -Oz -r $sChr:1-999999999 $VCFDIR/${i}.filtered.vcf.gz > ./qry.vcf.gz; \
		$BCFTOOLS index ./qry.vcf.gz; \
		$REVEEL index -b 1000000 qry.vcf.gz; \
		$REVEEL shortlist -n -M 2 -r $sChr:1-999999999 qry.vcf.gz qry.vcf.gz reveelcalled ;\
		java -jar $BEAGLE like=reveelcalled.$sChr.iter10.pr markers=reveelcalled.$sChr.iter10.marker out=prob;\
		$REVEEL merge qry.vcf.gz reveelcalled.iter10 reveelcalled prob.reveelcalled.$sChr.iter10.pr.dose.gz ) >  "$WDChr/log.txt" 2>&1 &

	done
done
wait



