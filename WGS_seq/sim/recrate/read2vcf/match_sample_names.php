<?php

$sRealBamList = "/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/NORv1.2LG/mapped/ORTWET.bamlist.txt";
$sSimPrefix = "/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/sim/recrate/pop1_and_pop2_reads.fq_parsed_ORTWET/ORTWET";

$sRealBamList = "/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/NORv1.2LG/mapped/ORTDRY.bamlist.txt";
$sSimPrefix = "/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/sim/recrate/pop1_and_pop2_reads.fq_parsed_ORTDRY/ORTDRY";


$sSimSuffix1 = "_read1.fastq.gz";
$sSimSuffix2 = "_read2.fastq.gz";
$nIndStartID = 1;
$sOutDIR = "simreads";
//$sOutDIR = "simreads_30x";

exec("mkdir -p $sOutDIR");

$arrIndIDs = fnLoadIDFromBamList($sRealBamList);

$nIDCount = $nIndStartID;

foreach($arrIndIDs as $sIDStem => $bDump) {
	exec("ln -sf $sSimPrefix$nIDCount$sSimSuffix1 $sOutDIR/$sIDStem"."_R1.fq.gz");
	exec("ln -sf $sSimPrefix$nIDCount$sSimSuffix2 $sOutDIR/$sIDStem"."_R2.fq.gz");
	$nIDCount++;
}

function fnLoadIDFromBamList($sRealBamList) {
	$arrRet = array();

	$h = fopen($sRealBamList  , 'r');
	while(false !== ($sLn = fgets($h))) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;

		$sBaseName = basename($sLn , ".bam");
		$arrRet[$sBaseName] = true;
	}

	return $arrRet;
}
?>
