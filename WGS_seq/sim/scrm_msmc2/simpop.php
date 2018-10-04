<?php

/*$sTest = '';
for($i=0;$i<1000000;$i++) $sTest .= fnGetAltBase('T');
echo(substr_count($sTest, 'A')."\n");
echo(substr_count($sTest, 'T')."\n");
echo(substr_count($sTest, 'C')."\n");
echo(substr_count($sTest, 'G')."\n");
die();
*/
$sScrm = "/beegfs/group_dv/software/source/scrm/scrm_ray_moredigits"; //I changed the version to print more digits in the positions to avoid collisions. or specify -p digits in the original version.
$sPop1Name = trim($argv[1]);
$sPop2Name = trim($argv[2]);
$nPop1Chr = intval($argv[3]);
$nPop2Chr = intval($argv[4]);
$nPersiteTheta = floatval($argv[5]);
$nPersiteRho = floatval($argv[6]);
$sAncChrFile = trim($argv[7]);
$sDemography = trim($argv[8]);
$sOutFile = trim($argv[9]);
$sGTRFile = trim($argv[10]);

$nTotalChr = $nPop1Chr + $nPop2Chr;
$sAncSeq = fnReadAnc($sAncChrFile);
$arrNormMutProbList = fnParseGTR($sGTRFile);
$nLocusLen = strlen($sAncSeq);

//die();

$hOut = popen("gzip -c > $sOutFile.gz" , 'w');
$hOutRaw = popen("gzip -c > ".$sOutFile.".raw.txt.gz" , 'w');

echo("Calling scrm (pop1 $sPop1Name $nPop1Chr pop2 $sPop2Name $nPop2Chr locus length $nLocusLen)...\n");
$sCmd = "$sScrm $nTotalChr 1 -t ".($nPersiteTheta * $nLocusLen). " -r ".($nPersiteRho * $nLocusLen)." $nLocusLen -I 2 $nPop1Chr $nPop2Chr $sDemography";
echo($sCmd."\n");
$hIn = popen($sCmd , 'r' );

$bResultReached = false;
$arrRelativePos = array();
$arrAbsPos = array();
$arrDerivedBases = array();
$bStartReadInd = false;
$nReadChr = 0;
while( false !== ($sLn = fgets($hIn) ) ) {
	$sLn = trim($sLn);

	if ($sLn == '') continue;

	if ( $sLn == '//') {
		$bResultReached = true;
		echo("scrm responded, waiting for output...\n");
		continue;
	}

	if (!$bResultReached) {
		echo($sLn."\n");
		continue;
	}

	//echo($sLn."\n");
	$arrF = explode(":" , $sLn);

	if (count($arrF) >=2 ) { // information field
		$sFieldName = trim($arrF[0]);
		$sFieldValue = trim($arrF[1]);
		if ($sFieldName == 'segsites') {
			echo("Simulated segregating sites: $sFieldValue\n");
			continue;
		}
		if ($sFieldName == 'positions') {
			$arrRelativePos = explode(' ', $sFieldValue);
			$arrAbsPos = fnRel2Abs($arrRelativePos);
			$arrDerivedBases = fnGetDerivedBases($sAncSeq,  $arrAbsPos);
			echo("Parsed positions of :". count($arrAbsPos )." sites\n");
			fwrite($hOutRaw, ">Positions\n".implode(',', $arrAbsPos)."\n" );
			fwrite($hOutRaw, ">DerivedBases\n".implode(',', $arrDerivedBases)."\n" );
			$bStartReadInd = true;
			continue;
		}
	} 

	if ($bStartReadInd) {
		$nReadChr++;
		$sPopName = ($nReadChr <= $nPop1Chr)? $sPop1Name : $sPop2Name; 
		$nChrInPop = ($nReadChr <= $nPop1Chr)? $nReadChr : $nReadChr - $nPop1Chr;
		echo("converting chr $nReadChr\n");
		$sNewChr = fnConvertChr2DNA($sLn,$sAncSeq, $arrAbsPos, $arrDerivedBases);
		//echo($sNewChr);
		fwrite($hOut , ">$sPopName.hap$nChrInPop\n$sNewChr\n");
		fwrite($hOutRaw, ">$sPopName.hap$nChrInPop\n$sLn\n");
	}

}

function fnReadAnc($sFile) {
	$hIn = fopen($sFile, 'r');
	$nRecordCount = 0;
	$sSeq = '';
	echo("Reading $sFile...\n");
	while( false !== ($sLn = fgets($hIn)) ) {
		$sLn = trim($sLn);
		if ($sLn[0] == '>') {
			$nRecordCount++;
			if ($nRecordCount > 1) {
				die("Error: the anc.fa file should contain a single fasta record for the ancestral sequence.\n");
			}
			continue;
		}

		if (preg_replace("/[ATCGatcg]/", "", $sLn) != '') {
			die("Error: DNA sequence can only contain ATCG, no ambiguous character is allowed.\n$sLn\n");
		}
		$sSeq .= $sLn;
	}

	return strtoupper($sSeq);
}

function fnRel2Abs($arrRelativePos) {
	global $nLocusLen;

	$arrRet = array();

	foreach($arrRelativePos as $nPos) {
		$arrRet[] = round(floatval($nPos) *  $nLocusLen);
	}

//check for conflicts due to numerical issues.

	$nCollisionPos = -1;
	$nCollisionCount = 0;
	foreach($arrRet as $nIdx => $nPos) {
		if ($nIdx == 0) {
			if ($nPos == 0) {
				$arrRet[$nIdx] = 1; //position is 1-based
			}

			continue;
		}

		if ( $arrRet[$nIdx-1] >= $nPos ) {
			if ($nCollisionPos == -1) {
				$nCollisionPos = $arrRet[$nIdx-1];
				$nCollisionCount = 0;
			}

			$nCollisionCount++;

			if (($arrRet[$nIdx]+$nCollisionCount) >= $nLocusLen) {
				echo("$nPos conflicts with the previous position, and it has reached the end of chromosome, discard\n");
			} else {
				echo("warning: position conflict with previous position $nCollisionPos, add $nCollisionCount to $nPos\n");
				$arrRet[$nIdx] += $nCollisionCount;
			}
		} else {
			$nCollisionPos = -1;
		}
	}

	return $arrRet;
}

function fnConvertChr2DNA($sGenotype,$sAncSeq, $arrAbsPos, $arrDerivedBases) {
	$sNewSeq = $sAncSeq;

	/*if (strlen($sGenotype) != count($arrAbsPos) ) {
		die("Error: the produced genotype has a different length compared to the 'positions'\n");
	}*/

	for($i=0;$i<count($arrAbsPos);$i++) {
		$nAllele = $sGenotype[$i];
		if ($nAllele == 0) {
			continue; //ancestral, no need to change.
		}

		$sNewSeq[intval($arrAbsPos[$i])-1] = $arrDerivedBases[$i];
	}

	return $sNewSeq;
}

function fnGetDerivedBases($sAncSeq,  $arrAbsPos) {
	$arrRet = array();

	foreach($arrAbsPos as $nIdx => $nPos) {
		$sBase = $sAncSeq[intval($nPos)-1];
		$arrRet[] = fnGetAltBase($sBase);
	}

	return $arrRet;
}

function fnGetAltBase($sAncBase) { //change it randomly to a different base
	global $arrNormMutProbList;

	if (!array_key_exists($sAncBase , $arrNormMutProbList) ) {
		$arrMutProb = array('A'=>0.25, 'T'=>0.25, 'C'=>0.25 , 'G'=>0.25);
		$arrNormMutProb = array();
		$nSumProb = 0;
		foreach($arrMutProb as $sBase => $nMutProb) {
			if ($sBase != $sAncBase) {
				$nSumProb += $nMutProb;
				$arrNormMutProb[$sBase] = $nMutProb;
			}
		}

		//normalize
		$nAccProb = 0;
		foreach($arrNormMutProb as $sBase => $nMutProb) {
			$nAccProb+= $nMutProb;
			$arrNormMutProb[$sBase] = $nAccProb / $nSumProb;
		}

		$arrNormMutProbList[$sAncBase] = $arrNormMutProb;
	}

	$arrNormMutProb = $arrNormMutProbList[$sAncBase];

	$nRd = mt_rand(0, mt_getrandmax()) / mt_getrandmax();
	//now decide which base to return:

	$sLastBase = '';
	foreach($arrNormMutProb as $sBase => $nAccProb) {
		$sLastBase = $sBase;
		if ($nRd <= $nAccProb) {
			return $sBase;
		}
	}

	return $sLastBase;
}

function fnParseGTR($sGTRFile) {
	$hIn = fopen($sGTRFile, 'r');
	$arrBases = array('A', 'T', 'G', 'C');
	$arrGTR = array();
	$arrFreq = array();

	foreach($arrBases as $sBase1) {
		$arrFreq[$sBase1] = 0;
		$arrGTR[$sBase1] = array();
		foreach($arrBases as $sBase2) {
			if ($sBase1 == $sBase2) continue;
			$arrGTR[$sBase1][$sBase2] = 0;
		}
	}
	
	while(false !== ($sLn = fgets($hIn)) ) {
		$sLn = trim($sLn);
		$arrF = explode("\t" , $sLn);
		if (count($arrF) == 3) {
			$arrGTR[$arrF[0]][$arrF[1]] = floatval($arrF[2]);
			$arrGTR[$arrF[1]][$arrF[0]] = floatval($arrF[2]);
		}

		if (count($arrF) == 2) {
			$arrFreq[$arrF[0]] = floatval($arrF[1]);
		}
	}

	//print_r($arrGTR);
	//print_r($arrFreq);

	$arrGTRWithFreq = array();
	foreach($arrGTR as $sB1 => $arrB2) {
		$arrGTRWithFreq[$sB1] = array();
		foreach($arrB2 as $sB2 => $nRate) {
			$arrGTRWithFreq[$sB1][$sB2] = $arrFreq[$sB2] * $nRate;
		}
	}

	//print_r($arrGTRWithFreq);

	$arrAccNorm = array();

	foreach($arrGTRWithFreq as $sB1 => $arrB2) {
		$nSumRate = array_sum($arrB2);
		$nAcc = 0;
		$arrAccNorm[$sB1] = array();
		foreach($arrB2 as $sB2 => $nRate) {
			$nAcc += $nRate;
			$arrAccNorm[$sB1][$sB2] = $nAcc / $nSumRate;
		}
	}

	//print_r($arrAccNorm);

	return($arrAccNorm);
}

?>
