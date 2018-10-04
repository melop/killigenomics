<?php
//produce a coverage mask for smc++
//$BCFTOOLS  mpileup  -q 60 /beegfs/group_dv/home/RCui/killifish_genomes/denovo/discovardenovo/NOR/LG/v1.2/FTH.chr.fa lnbams/${i}/*.bam | less

$BCFTOOLS = "/beegfs/group_dv/software/source/samtools-1.6/bin/bin/bcftools";

$sOutgroup = "/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/SimNORLG1/ref.fa";
$sPop = $argv[1];//"ORTWET";
$nMinVarQual = 60;

$sOutDIR = "masks_2pops";

exec("mkdir -p $sOutDIR");

$sOutMask = "$sOutDIR/$sPop.mask.varqual.bed";
$hOutBed = fopen($sOutMask , 'w');

#ORTDRY 2.229171        59      263.042178
#ORTWET 2.230603        60      267.67236
#RACDRY 2.572075        59      303.50485
#RACWET 2.251057        57      256.620498



$sCmd = "$BCFTOOLS  view ./vcfs_2pops/$sPop.var.vcf.gz";
echo($sCmd."\n");
$hIn = popen($sCmd , 'r');
$sCurrentChr = "";
$arrCovForChr = array(); //key is position, value is coverage.
$arrChrLens = array();

do {
	$sLn = fgets($hIn);
	$arrF = ($sLn!==false)? explode("\t", trim($sLn)):array();
	if ($sLn !== false && $sLn[0] == '#') {
		//##contig=<ID=chr1,length=72567360>
		if (substr($sLn, 0, 9 ) == '##contig=') {
			preg_match("/ID=(\S+),length=(\d+)/", trim($sLn), $arrM);
			$arrChrLens[$arrM[1]] = intval($arrM[2]);
		}
		continue;
	}

	if ($sLn === false || $sCurrentChr !=$arrF[0] ) {

		if ($sLn === false) break;

		$sCurrentChr = $arrF[0];
		$arrCovForChr = array(); //key is position, value is coverage.
	}


	if ( intval($arrF[5]) < $nMinVarQual || $arrF[3] == 'N') {
		fwrite($hOutBed , "$arrF[0]\t".($arrF[1] - 1)."\t$arrF[1]\n");
                continue;
	}


	if (intval($arrF[1]) % 10000 == 0) {
		echo($arrF[1]."\n");
	}

	
} while(true);

exec("rm $sOutMask.gz; bgzip $sOutMask; tabix -p bed $sOutMask.gz ");



?>
