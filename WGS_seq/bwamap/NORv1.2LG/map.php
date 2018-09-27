<?php

$sFQDIR = "../../joinruns/joined/";
$sOutDIR = "./mapped";

$nThisPart = 0;
$nTotalPart = 1;

$nPerJobThread = 20;

$sRefGenome = "ref.fa";
$nMinMapQ = 30; // skip low quality alignments


$sBWA = "/beegfs/group_dv/software/bin/bwa.kit/bwa";
$SAMTOOLS = "samtools1.2";
$SAM2PRO = "/beegfs/group_dv/software/source/PopLD/Sam2pro_0.6/sam2pro";
$PICARDDIR="/beegfs/group_dv/software/source/picard-tools-1.119/";
$GATKPATHNEW="/beegfs/group_dv/software/source/gatk3.4.46/";

$SCRATCHDIR=getcwd()."/tmp/";
exec("mkdir -p $SCRATCHDIR");


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

$arrR1 = glob("$sFQDIR/*_R1.fq.gz");

for($i=0;$i<count($arrR1);$i++) {
	if ( $i % $nTotalPart != $nThisPart) continue;

	$sBAMDir = "$sOutDIR/bam";
	$sRealignedBamDir = "$sOutDIR/realigned_bam";
	$sMPileUpDir = "$sOutDIR/mpileup";
	$sProfileDir = "$sOutDIR/pro";

	exec("mkdir -p $sBAMDir ");
	exec("mkdir -p $sRealignedBamDir");
	exec("mkdir -p $sMPileUpDir ");
	exec("mkdir -p $sProfileDir ");

	$sR1 = realpath($arrR1[$i]);
	$sFQDIR = realpath($sFQDIR);

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

	$bRealignedBamDoneFlag = "$sRealignedBamDir/$sFileStem.realigned.done";
	$sSpSortedBam = $sSpSortedBamPrefix.".bam";
	$sSpRGBam = "$sRealignedBamDir/$sFileStem.sorted.RG.bam";
	$sSpMarkDupBam = "$sRealignedBamDir/$sFileStem.sorted.markdup.bam";
	$sSpDupMetrics = "$sRealignedBamDir/$sFileStem.dupmetrics.log";
	$sSpRealnIntervals = "$sRealignedBamDir/$sFileStem.realign.intervals";
	$sSpRealignedBam = "$sRealignedBamDir/$sFileStem.realigned.bam";

	//indel realignment
	$sCmd = "if [ -e $bRealignedBamDoneFlag ]; then echo $bRealignedBamDoneFlag was finished. skip; else ";
	$sCmd .= " module load Java; ";
	$sCmd .= "java -Dsnappy.disable=true -Xmx12g -jar $PICARDDIR/AddOrReplaceReadGroups.jar I=$sSpSortedBam O=$sSpRGBam  SORT_ORDER=coordinate RGID=$sFileStem RGLB=$sFileStem RGPL=illumina RGSM=$sFileStem RGPU=$sFileStem CREATE_INDEX=True VALIDATION_STRINGENCY=SILENT TMP_DIR=$SCRATCHDIR ; ";
	$sCmd .= "java -Dsnappy.disable=true -Xmx12g -jar $PICARDDIR/MarkDuplicates.jar MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=900  INPUT=$sSpRGBam OUTPUT=$sSpMarkDupBam ASSUME_SORTED=true METRICS_FILE=$sSpDupMetrics VALIDATION_STRINGENCY=SILENT TMP_DIR=$SCRATCHDIR; "; //markdup
	$sCmd .= "$SAMTOOLS index $sSpMarkDupBam; ";
	$sCmd .= "java -Xmx12g -jar $GATKPATHNEW/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $sRefGenome -I $sSpMarkDupBam -o $sSpRealnIntervals ; "; //markdup
	$sCmd .= "java -Xmx16g -jar $GATKPATHNEW/GenomeAnalysisTK.jar -T IndelRealigner -R $sRefGenome -I $sSpMarkDupBam -targetIntervals $sSpRealnIntervals -o $sSpRealignedBam --maxReadsForRealignment 12000 "; //markdup
	$sCmd .= " && touch $bRealignedBamDoneFlag; fi;";


	echo($sCmd."\n");
	exec($sCmd);

	if (file_exists($sSpRGBam) ) { unlink($sSpRGBam);}
	if (file_exists($sSpMarkDupBam)) {unlink($sSpMarkDupBam);}

	if (!file_exists($bRealignedBamDoneFlag)) {
		die("$bRealignedBamDoneFlag job failed.\n");
	}

/*
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
*/

}

?>
