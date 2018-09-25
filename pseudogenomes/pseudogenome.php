<?php
$sMainList = "list.txt";
//$sMainList = "list2.txt";
$sOutDIR = "out";
$BWA = "bwa";
$SAMTOOLS = "samtools1.2";
$BCFTOOLS = "bcftools1.2";
$nPerJobThread = 40;
$PICARDDIR="/software/source/picard-tools-1.119/";
$GATKPATHNEW="/software/source/gatk3.4.46/";
$SCRATCHDIR=getcwd()."/tmp/";

exec("mkdir -p $sOutDIR");
exec("mkdir -p $SCRATCHDIR");

$sBamDIR = "$sOutDIR/bams";
$sBcfDIR = "$sOutDIR/bcfs";
$sPseudoGenomeDIR = "$sOutDIR/pseudogenomes";

exec("mkdir -p $sBamDIR");

list($arrRefs , $arrSpecies) = fnLoadList($sMainList); 

//print_r($arrRefs);
//print_r($arrSpecies);

$arrJobIDs = array();
//spawn mapping jobs:
foreach($arrSpecies as $sRef => $arrSppInRef) {
	$sRefGenome = $arrRefs[$sRef];

	foreach( $arrSppInRef as $sSp => $arrReads) {
		$sSpBAMDIR = "$sBamDIR/$sRef/$sSp";
		$sSpSortedBamPrefix = "$sSpBAMDIR/sorted";
		$sBWABatch = "$sSpBAMDIR/bwa.sbatch";
		$bMapDoneFlag = "$sSpBAMDIR/map.done";
		exec("mkdir -p $sSpBAMDIR");
		$sCmd = "";
		$sR1File = "";
		$sR2File = "";
		if (count($arrReads['r1']) == 1 ) {
			$sR1File = $arrReads['r1'][0];
			$sR2File = $arrReads['r2'][0];
		} else {
			$sR1File = "$sSpBAMDIR/concat_R1.fq";
			$sR2File = "$sSpBAMDIR/concat_R2.fq";
			$sCmd .= "zcat -f " . implode(" " , $arrReads['r1']) . " > $sR1File & \n" ;
			$sCmd .= "zcat -f " . implode(" " , $arrReads['r2']) . " > $sR2File & \nwait;\n" ;
		}

		$sCmd .= "if [ -e $bMapDoneFlag ];then echo $bMapDoneFlag was finished. skip.; else $BWA mem -t $nPerJobThread \"$sRefGenome\" \"$sR1File\" \"$sR2File\" | $SAMTOOLS view -Su - | $SAMTOOLS sort - $sSpSortedBamPrefix && touch $bMapDoneFlag; fi;";
		if (!file_exists($bMapDoneFlag)) {
			fnSubmitSlurm($sCmd , $nPerJobThread, '50G', $sBWABatch, false, $arrJobIDs); //append job id to job ids
		} else {
			echo("$bMapDoneFlag is done, skip.\n");
		}

	}
}

echo("Mapping jobs submitted, waiting for return...\n");
fnWaitForJobs($arrJobIDs);
echo("All Mapping jobs returned...\n");

//now the mapping should be done. get bcf files!
foreach($arrSpecies as $sRef => $arrSppInRef) {
	$sRefGenome = $arrRefs[$sRef];

	foreach( $arrSppInRef as $sSp => $arrReads) {
		$sSpBAMDIR = "$sBamDIR/$sRef/$sSp";
		$sSpSortedBam = "$sSpBAMDIR/sorted.bam";
		$bMapDoneFlag = "$sSpBAMDIR/map.done";

		if (!file_exists($bMapDoneFlag )) {
			echo("Mapping has failed for $bMapDoneFlag.\n ");
			continue;
		}
		
		$sSpBCFDIR = "$sBcfDIR/$sRef/$sSp";
		$sSpRGBam = "$sSpBCFDIR/sorted.RG.bam";
		$sSpMarkDupBam = "$sSpBCFDIR/sorted.markdup.bam";
		$sSpDupMetrics = "$sSpBCFDIR/dupmetrics.log";
		$sSpRealnIntervals = "$sSpBCFDIR/realign.intervals";
		$sSpRealignedBam = "$sSpBCFDIR/realigned.bam";
		$sSpRawBcf = "$sSpBCFDIR/raw.bcf"; 
		$sSpCalledBcf = "$sSpBCFDIR/called.bcf"; 
		$sBCFBatch = "$sSpBCFDIR/bcf.sbatch";
		$sBCFDoneFlag = "$sSpBCFDIR/bcf.done";

		exec("mkdir -p $sSpBCFDIR");
		$sCmd = "module load Java;\n";

		$sCmd .= "java -Dsnappy.disable=true -Xmx12g -jar $PICARDDIR/AddOrReplaceReadGroups.jar I=$sSpSortedBam O=$sSpRGBam  SORT_ORDER=coordinate RGID=$sSp RGLB=$sSp RGPL=illumina RGSM=$sSp RGPU=$sSp CREATE_INDEX=True VALIDATION_STRINGENCY=SILENT TMP_DIR=$SCRATCHDIR ;\n";
		$sCmd .= "java -Dsnappy.disable=true -Xmx12g -jar $PICARDDIR/MarkDuplicates.jar MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=900  INPUT=$sSpRGBam OUTPUT=$sSpMarkDupBam ASSUME_SORTED=true METRICS_FILE=$sSpDupMetrics VALIDATION_STRINGENCY=SILENT TMP_DIR=$SCRATCHDIR;\n"; //markdup
		$sCmd .= "$SAMTOOLS index $sSpMarkDupBam;\n";
		$sCmd .= "java -Xmx12g -jar $GATKPATHNEW/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $sRefGenome -I $sSpMarkDupBam -o $sSpRealnIntervals ;\n"; //markdup
		$sCmd .= "java -Xmx16g -jar $GATKPATHNEW/GenomeAnalysisTK.jar -T IndelRealigner -R $sRefGenome -I $sSpMarkDupBam -targetIntervals $sSpRealnIntervals -o $sSpRealignedBam --maxReadsForRealignment 120000\n"; //markdup
		$sCmd .= "$SAMTOOLS mpileup --count-orphans -t DP,DV,DPR,INFO/DPR,DP4,SP --adjust-MQ 0 --min-MQ 20 --min-BQ 25 -go $sSpRawBcf -f $sRefGenome $sSpRealignedBam ;\n";
		$sCmd .= "$BCFTOOLS call -m -M -A  -O b -o $sSpCalledBcf $sSpRawBcf;\n$BCFTOOLS index $sSpCalledBcf && touch $sBCFDoneFlag\n";
		if (!file_exists($sBCFDoneFlag)) {
			fnSubmitSlurm($sCmd , $nPerJobThread, '50G', $sBCFBatch, false, $arrJobIDs); //append job id to job ids
		} else {
			echo("$sBCFDoneFlag exists, skip\n");
		}

	}
}
echo("Variant calling jobs submitted, waiting for return...\n");
fnWaitForJobs($arrJobIDs);
echo("All Variant calling jobs returned...\n");


//finally make the pseudogenomes
//now the mapping should be done. get bcf files!
foreach($arrSpecies as $sRef => $arrSppInRef) {
	$sRefGenome = $arrRefs[$sRef];
	$sPGDir = "$sPseudoGenomeDIR/$sRef";
	exec("mkdir -p $sPGDir");

	$sAllTaxa = "$sPGDir/alltaxa.txt";
	$hAllTaxa = fopen($sAllTaxa , 'w');

	fwrite($hAllTaxa , "$sRef\t$sRefGenome\n");

	foreach( $arrSppInRef as $sSp => $arrReads) {
		$sSpBCFDIR = "$sBcfDIR/$sRef/$sSp";
		$sBCFDoneFlag = "$sSpBCFDIR/bcf.done";
		$sSpCalledBcf = realpath("$sSpBCFDIR/called.bcf"); 

		if (!file_exists($sBCFDoneFlag )) {
			echo("Variant call has failed for $sBCFDoneFlag.\n ");
			continue 2;
		}
		fwrite($hAllTaxa , "$sSp\t$sSpCalledBcf\n");

	}

	//now link stuff
	exec("cp `pwd`/getAlignmentAllowN.php $sPGDir/");
	exec("cp `pwd`/phy2fasta.php $sPGDir/");
	exec("cp `pwd`/lib.php $sPGDir/");
	exec("cp `pwd`/config.php $sPGDir/");
	exec("cp `pwd`/codons.php $sPGDir/");
	exec("cp `pwd`/getpseudogenome.sh $sPGDir/");
	exec("cp charRequirements.txt $sPGDir/");

	$sPGBatch = "$sPGDir/sub.sbatch"; 
	fnSubmitSlurm("bash getpseudogenome.sh" , 40, '250G', $sPGBatch, $sPGDir, $arrJobIDs); //append job id to job ids
}

echo("producing pseudogenomes.\n");
fnWaitForJobs($arrJobIDs);
echo("All done.\n");

function fnLoadList($sMainList) {

	$hIn = fopen($sMainList , 'r');
	fgets($hIn); //discard header line
	$arrRetRef = array();
	$arrRetSp = array();
	while( false !== ( $sLn = fgets($hIn) )) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		$arrF = explode("\t" , $sLn);
		if (count($arrF) != 7) continue;
		if ($arrF[1] == $arrF[2]) { // ref def.
			$arrRetRef[$arrF[1]] = $arrF[5];
			if (!file_exists($arrF[5].".sa") ) {
				die("Reference genome not yet indexed by bwa, please do that! : " . $arrF[5] ."\n");
			} 
		} else {
			if (!array_key_exists($arrF[1] , $arrRetSp )) {
				$arrRetSp[$arrF[1]] = array();
			}

			if (array_key_exists($arrF[2], $arrRetSp[$arrF[1]]) ) {
				die("Error: species " . $arrF[2] . " already defined!\n");
			}

			$arrRetSp[$arrF[1]][$arrF[2]] = array('r1' => explode(',',$arrF[5]) , 'r2' => explode(',' , $arrF[6]) );
		}
	}

	return array($arrRetRef , $arrRetSp);
}

function fnSubmitSlurm($sCmd , $nPerJobThread, $sMem, $sBatch, $sWD = false, &$arrJobIDs) {
	if ($sWD === false) {
		$sWD = getcwd();
	}
	$sHeader = "#!/bin/bash\n#SBATCH -n 1\n#SBATCH -c $nPerJobThread\n#SBATCH --mem=$sMem\n#SBATCH -p blade\ncd " . $sWD ."\n";
	$hBatch = fopen($sBatch , 'w');
	fwrite($hBatch , $sHeader . $sCmd ."\n");
	//return;
	exec("sbatch $sBatch" , $arrOut);
	if (count($arrOut) < 1) {
		die("Failed to submit slurm job.\n$sCmd\n");
	}

	preg_match("/\d+/" , $arrOut[0], $arrM);
	if (count($arrM) != 1) {
		die("Cannot understand sbatch output\n");
	}

	$sJobID = $arrM[0];
	echo("submitted: $sJobID\n");
	$arrJobIDs[$sJobID] = true;
}

function fnWaitForJobs(&$arrJobIDs) {
	while(count($arrJobIDs) > 0 ) {
		foreach($arrJobIDs as $sJobID => $bDum) {
			exec("squeue -j $sJobID 2>/dev/null | wc -l " , $arrOut);
			if (trim($arrOut[0]) == '0') { // if job is done
				unset($arrJobIDs[$sJobID]);
			}
		}

		sleep(30);
	}
}
?>
