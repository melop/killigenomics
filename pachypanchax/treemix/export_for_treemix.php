<?php
//prepare a coverage mask file, read in bcf files, and find the peak coverage. exclude regions with abnormal sequence coverage.
//print_r(fnCountAlleleReads('A', 'T', 'GT:DP:DV:SP:DP4:DPR' , array(' 0/0:222:0:0:164,58,1,2:222') ) );
//die();
$sBCFTOOLS = "/beegfs/group_dv/software/source/bcftools-1.6/bcftools";
$sSampleDef = "samples.txt";
$sOutDir = "fortreemix/";

$sSampleDef = "samples_PLPonly.txt";
$sOutDir = "fortreemix_PLPonly/";

$nRecRate = 100*1000* 24*2/860e6  ; //cM/kb. 24 chromsomes, 2 recombinations per meisosis per chromosome


$nBlockSize = 1e6;

exec("mkdir -p $sOutDir");

$sChrFilter = 'chrLG';
$nDoPart = $argv[1];
$nTotalPart = $argv[2];
$nMinQUAL = 30; //min variant quality.

$arrSamples = fnLoadSamples($sSampleDef);
$arrRefNames = array_keys($arrSamples['ref']);
$arrChrInfo = fnLoadChrInfo($arrSamples['ref'][$arrRefNames[0]][0]);

//print_r($arrChrInfo );

$nChrCount = -1;

$nRecRate = $nRecRate / 100 / 1000; //Morgan/bp

foreach($arrChrInfo as $sChrName => $nChrLen) {
	$nChrCount++;
	if ( ($nChrCount % $nTotalPart) != $nDoPart) {
		continue;
	}
	echo("Doing $sChrName\t$nChrLen\n");

	$nPrevPos = 0;
	$hOut = fopen("$sOutDir/in_$sChrName.txt", 'w');

	for($nBlockStart=1;$nBlockStart<$nChrLen;$nBlockStart+=$nBlockSize) {
		$nBlockEnd = $nBlockStart + $nBlockSize - 1;
		if ($nBlockEnd > $nChrLen) {
			$nBlockEnd = $nChrLen;
		}
		//read in refs
		echo("\t$sChrName:$nBlockStart-$nBlockEnd\n");
		$arrRefData = array();
		foreach($arrSamples['ref'] as $sRefPop => $arrRefBCFs) {
			$arrRefData[$sRefPop] = fnLoadRefData($arrRefBCFs , "$sChrName:$nBlockStart-$nBlockEnd" , "cov_ref_$sRefPop.cov.cutoffs.txt");
		}

		if ($nChrCount == 0 && $nBlockStart==1) {
			fwrite($hOut , implode("\t", array_keys($arrSamples['ref']) )."\n" );
		}


		$arrCommonPos = array();
		foreach($arrRefData as $sRefPop => $arrRefPopData) {
			if (count($arrCommonPos) == 0) {
				$arrCommonPos = array_keys($arrRefPopData);
				continue;
			}
			$arrCommonPos = array_intersect(array_keys($arrRefPopData), $arrCommonPos);
		}


		echo("Found ".count($arrCommonPos)." common positions on chr $sChrName\n" );


		foreach($arrCommonPos as $nPos) {
			$arrOut = array();
			//$arrOut[] = $sChrName;
			//$arrOut[] = $nPos;
			$arrAlleles = array();

			foreach($arrRefData as $sRefPop => $arrRefPopData) {
				list($sThisAllele1, $sThisAllele2) = array_keys($arrRefPopData[$nPos]);
				if ($sThisAllele1!='.' && $sThisAllele1!='N' && $arrRefPopData[$nPos][$sThisAllele1] != 0) {
					$arrAlleles[$sThisAllele1] = true;
				} 

				if ($sThisAllele2!='.' && $sThisAllele2!='N' && $arrRefPopData[$nPos][$sThisAllele2] != 0) {
					$arrAlleles[$sThisAllele2] = true;
				} 

			}

			if (count($arrAlleles) != 2 ) {
				continue; //monomorphic, skip
			}


			$arrAlleles = array_keys($arrAlleles);
			//echo("Pos: $nPos\n");
			foreach($arrRefData as $sRefPop => $arrRefPopData) { 
				//echo($sRefPop);
				//print_r($arrRefPopData[$nPos]);
				$nCount1 = array_key_exists($arrAlleles[0], $arrRefPopData[$nPos])? $arrRefPopData[$nPos][$arrAlleles[0]] : 0;
				$nCount2 = array_key_exists($arrAlleles[1], $arrRefPopData[$nPos])? $arrRefPopData[$nPos][$arrAlleles[1]] : 0;
				$arrOut[] = "$nCount1,$nCount2";
			}

			//die();
			//$arrOut[] = ($nPos - $nPrevPos) * $nRecRate;
			//$nPrevPos = $nPos; 


			fwrite($hOut, implode("\t", $arrOut )."\n");
		}

	}

}

function fnLoadRefData($arrSamples, $sChrName, $sCovCutoffFile) {
	global $sBCFTOOLS, $nMinQUAL;
	list($nDPMinCutoff, $nDPMaxCutoff) = fnLoadCovCutoffs($sCovCutoffFile);
	$nDPMinCutoff = floor($nDPMinCutoff);
	$nDPMaxCutoff = ceil($nDPMaxCutoff);

	$sCmd = '';
	if (count($arrSamples) > 1) {
		$sCmd .= "$sBCFTOOLS merge ";
	} else {
		$sCmd .= "$sBCFTOOLS view ";
	}

	$sCmd .= " -r $sChrName ";
	$sCmd .= implode(' ', $arrSamples);
	//$sCmd .= " | head -n " . ($nSampleLines + 10000);
	echo($sCmd."\n");

	$descriptorspec = array(
	   0 => array("pipe", "r"),  // stdin is a pipe that the child will read from
	   1 => array("pipe", "w"),  // stdout is a pipe that the child will write to
	   2 => array("pipe", "w") // stderr is a file to write to
	);

	$hProc = proc_open($sCmd, $descriptorspec, $pipes);

	$arrRet = array();
	$nLnCount=0;
	while( false !== ($sLn = fgets($pipes[1]) ) ) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		if ($sLn[0] == '#') continue;
		$arrF = explode("\t", $sLn);
		if (count($arrF) < 8) continue;

		if ($arrF[5] < $nMinQUAL) continue; //quality filter

		preg_match('/DP=([0-9]+)/', $arrF[7], $arrM);
		if (count($arrM) != 2) continue;
		$nDP = intval($arrM[1]);
		if ($nDP < $nDPMinCutoff || $nDP > $nDPMaxCutoff) {
			continue;
		}

		//now parse genotypes
		if (strlen($arrF[3])>1 || strlen($arrF[4])>1 ) {
			continue; //ignore indels.
		}

		$arrRet[$arrF[1]] = fnCountAlleles($arrF[3], $arrF[4], $arrF[8], array_slice($arrF, 9) );
		//echo($nLnCount."\n");

		$nLnCount++;
		if ($nLnCount % 0.1e6 == 0 ) {
			echo("Proccessed ".($nLnCount / 1e6). "M positions\n");
		}

	}

	return $arrRet;
}

function fnLoadSampleData($arrSamples , $sChrName , $sCovCutoffFile) {
	global $sBCFTOOLS, $nMinQUAL;
	list($nDPMinCutoff, $nDPMaxCutoff) = fnLoadCovCutoffs($sCovCutoffFile);
	$nDPMinCutoff = floor($nDPMinCutoff);
	$nDPMaxCutoff = ceil($nDPMaxCutoff);

	$sCmd = '';
	if (count($arrSamples) > 1) {
		$sCmd .= "$sBCFTOOLS merge ";
	} else {
		$sCmd .= "$sBCFTOOLS view ";
	}

	$sCmd .= " -r $sChrName ";
	$sCmd .= implode(' ', $arrSamples);
	//$sCmd .= " | head -n " . ($nSampleLines + 10000);
	echo($sCmd."\n");

	$descriptorspec = array(
	   0 => array("pipe", "r"),  // stdin is a pipe that the child will read from
	   1 => array("pipe", "w"),  // stdout is a pipe that the child will write to
	   2 => array("pipe", "w") // stderr is a file to write to
	);

	$hProc = proc_open($sCmd, $descriptorspec, $pipes);

	$arrRet = array();
	$nLnCount=0;
	while( false !== ($sLn = fgets($pipes[1]) ) ) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		if ($sLn[0] == '#') continue;
		$arrF = explode("\t", $sLn);
		if (count($arrF) < 8) continue;

		if ($arrF[5] < $nMinQUAL) continue; //quality filter

		preg_match('/DP=([0-9]+)/', $arrF[7], $arrM);
		if (count($arrM) != 2) continue;
		$nDP = intval($arrM[1]);
		if ($nDP < $nDPMinCutoff || $nDP > $nDPMaxCutoff) {
			continue;
		}

		//now parse genotypes
		if (strlen($arrF[3])>1 || strlen($arrF[4])>1 ) {
			continue; //ignore indels.
		}

		$arrRet[$arrF[1]] = fnCountAlleleReads($arrF[3], $arrF[4], $arrF[8], array_slice($arrF, 9) );
		//echo($nLnCount."\n");
		$nLnCount++;
		if ($nLnCount % 1e6 == 0 ) {
			echo("Proccessed ".($nLnCount / 1e6). "M positions\n");
		}

	}

	return $arrRet;
}

function fnLoadSamples($sSampleDef) {
	$h = fopen($sSampleDef, 'r');

	$arrRet = array('ref' => array(), 'samples' => array());
	while(false !== ($sLn =fgets($h))) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		if ($sLn[0] == '#') continue;

		list($sRef, $sPop, $sSampleID, $sBCF) = explode("\t", $sLn );

		$sType = ($sRef == 'Ref')? 'ref':'samples'; 
		if (!array_key_exists($sPop, $arrRet[$sType]) ) {
			$arrRet[$sType][$sPop] = array();
		}

		$arrRet[$sType][$sPop][] = $sBCF;		
	}

	return $arrRet;
}

function fnCountAlleles($sRefAllele, $AltAllele, $sDataDef, $arrSampleData ) {
	$arrDataDefs = array_flip(explode(':', $sDataDef));
	if (!array_key_exists('GT',$arrDataDefs ) ) {
		return false;
	}

	$nRefCount = 0;
	$nAltCount = 0;
	foreach($arrSampleData as $sDat) {
		$arrF = explode(':', $sDat);
		$sGT = $arrF[$arrDataDefs['GT']];
		$arrAlleles = preg_split('/[\/|]/', $sGT);
		if (count($arrAlleles) != 2) {
			return false;
		}

		if ($arrAlleles[0] == '.' || $arrAlleles[1] == '.') {
			continue;
		}

		$nAlt = $arrAlleles[0] + $arrAlleles[1];
		$nRef = 2 - $nAlt;

		$nRefCount += $nRef;
		$nAltCount += $nAlt;
	}

	return array($sRefAllele => $nRefCount , $AltAllele => $nAltCount);
}

function fnCountAlleleReads($sRefAllele, $AltAllele, $sDataDef, $arrSampleData ) {
	$arrDataDefs = array_flip(explode(':', $sDataDef));
	if (!array_key_exists('DP4',$arrDataDefs ) ) {
		return false;
	}

	$nRefCount = 0;
	$nAltCount = 0;
	foreach($arrSampleData as $sDat) {
		$arrF = explode(':', $sDat);
		$sGT = $arrF[$arrDataDefs['DP4']];
		$arrAlleles = explode(',', $sGT);
		if (count($arrAlleles) != 4) {
			return false;
		}

		if ($arrAlleles[0] == '.' || $arrAlleles[1] == '.' || $arrAlleles[2] == '.' || $arrAlleles[3] == '.' ) {
			continue;
		}

		$nRef = $arrAlleles[0] + $arrAlleles[1];
		$nAlt = $arrAlleles[2] + $arrAlleles[3];

		$nRefCount += $nRef;
		$nAltCount += $nAlt;
	}

	return array($sRefAllele => $nRefCount , $AltAllele => $nAltCount);
}


function fnLoadChrInfo($sBCF) {
	//echo("$sBCF\n");
	global $sBCFTOOLS, $sChrFilter;
	$sCmd = "$sBCFTOOLS view -h $sBCF";
	$h = popen($sCmd, 'r');
	$arrRet = array();
	while(false !== ($sLn = fgets($h) ) ) {
		$sLn = trim($sLn);
		if (substr($sLn, 0, 9) == '##contig=') {
			preg_match('/ID=([^,]+).*length=([0-9]+)/', $sLn, $arrM);
			if (count($arrM) != 3 ) continue;
			$sChrName = $arrM[1];
			$nChrLen = $arrM[2];
			if (strpos($sChrName , $sChrFilter) === false ) {
				continue;
			}

			$arrRet[$sChrName] = $nChrLen;
		}
	}

	pclose($h);
	return $arrRet;
}

function fnLoadCovCutoffs($sCovCutoffFile) {
	return array_slice(explode("\t", trim(file_get_contents($sCovCutoffFile))), 2,2);
}

?>
