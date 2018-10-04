<?php

$sFQDIR = "./simreads/";
$sOutDIR = "./mapped";

$nThisPart = 0;
$nTotalPart = 1;

$nPerJobThread = 2;

$sRefGenome = "ref.fa";
$nMinMapQ = 30; // skip low quality alignments


$sBWA = "/beegfs/group_dv/software/bin/bwa.kit/bwa";
$SAMTOOLS = "samtools1.2";
$SAM2PRO = "/beegfs/group_dv/software/source/PopLD/Sam2pro_0.6/sam2pro";


while(count($argv) > 0) {
	$arg = array_shift($argv);
	switch($arg) {
		case '-q':
			$sFQDIR = trim(array_shift($argv));
			break;
		case '-o':
			$sOutDIR  = trim(array_shift($argv));
			break;
		case '-N':
			$nTotalPart  = trim(array_shift($argv));
			break;
		case '-f':
			$nThisPart  = trim(array_shift($argv));
			break;

	}
}

$sFQDIR = realpath($sFQDIR);
$arrR1 = glob("$sFQDIR/*_R1.fq.gz");

for($i=0;$i<count($arrR1);$i++) {
	if ( $i % $nTotalPart != $nThisPart) continue;

	$sBAMDir = "$sOutDIR/bam";
	$sMPileUpDir = "$sOutDIR/mpileup";
	$sProfileDir = "$sOutDIR/pro";

	exec("mkdir -p $sBAMDir ");
	exec("mkdir -p $sMPileUpDir ");
	exec("mkdir -p $sProfileDir ");

	$sR1 = $arrR1[$i];

	$sR1FileName = basename($sR1);
	$sR2FileName = str_replace("_R1.fq.gz" , "_R2.fq.gz" , $sR1FileName);
	$sR2 = $sFQDIR."/$sR2FileName";

	$sFileStem = str_replace("_R1.fq.gz"  , '', $sR1FileName);

	if (!file_exists($sR2)) {
		die("ERROR: $sR2 not found");
	}

	$bMapDoneFlag = "$sBAMDir/$sFileStem.map.done";
	$sSpSortedBamPrefix = "$sBAMDir/$sFileStem";

	$sCmd = "if [ -e $bMapDoneFlag ];then echo $bMapDoneFlag was finished. skip.; else $sBWA mem -t $nPerJobThread \"$sRefGenome\" \"$sR1\" \"$sR2\" | $SAMTOOLS view -Su - | $SAMTOOLS sort - $sSpSortedBamPrefix && touch $bMapDoneFlag; fi;";
	echo($sCmd."\n");
	exec($sCmd);

	if (!file_exists($bMapDoneFlag)) {
		die("$bMapDoneFlag job failed.\n");
	}

	$bMPileUpDoneFlag = "$sMPileUpDir/$sFileStem.mpileup.done";
	$sMPileUpFile = "$sMPileUpDir/$sFileStem.mpileup";

	$sCmd = "if [ -e $bMPileUpDoneFlag ];then echo $bMPileUpDoneFlag was finished. skip.; else $SAMTOOLS mpileup -f $sRefGenome -q $nMinMapQ \"$sBAMDir/$sFileStem.bam\" >  $sMPileUpFile && touch $bMPileUpDoneFlag; fi;";
	echo($sCmd."\n");
	exec($sCmd);


	$bProfileFlag = "$sProfileDir/$sFileStem.pro.done";
	$sProfileFile= "$sProfileDir/$sFileStem.pro";

	$sCmd = "if [ -e $bProfileFlag ];then echo $bProfileFlag was finished. skip.; else cat $sMPileUpFile | $SAM2PRO -c 6 -m 1 >  $sProfileFile && touch $bProfileFlag; fi;"; //fixed bug 29-1-2018, -m coverage by default was 4, should set this to 1
	echo($sCmd."\n");
	exec($sCmd);

}

?>
