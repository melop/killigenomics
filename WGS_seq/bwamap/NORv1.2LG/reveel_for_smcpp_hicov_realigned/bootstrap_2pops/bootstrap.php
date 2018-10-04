<?php

$sRawDir = $argv[1];
$nBlockSize = floatval($argv[2]); //5e6; //cut into 5Mb blocks, like the PSMC paper
$nReps = intval($argv[3]);
$sBaseOutDir  = trim( $argv[4]);
$sOutDir = $sBaseOutDir."/".basename($sRawDir);

if ($nBlockSize < 1e6) {
	die("Block size is too small $nBlockSize\n");
}

exec("mkdir -p $sOutDir");

$bSplitDoneFlag = "$sOutDir/split.done";
$sSplitDIR = "$sOutDir/split";

if (!file_exists($bSplitDoneFlag) ) {
	//now split each chromosome into pieces

	$arrSMC = glob("$sRawDir/*/smc*.gz");

	foreach($arrSMC as $sSMCFile) {
		$sChr = basename(dirname($sSMCFile));
		$sChrOutDir = $sSplitDIR."/$sChr/";
		exec("mkdir -p $sChrOutDir"); 

		$h = popen("zcat -f $sSMCFile", 'r');
		$sHeader = "";
		$nPart = 0;
		$sContent = "";
		$nPartLen = 0;
		$nTotalLen = 0;
		
		$sSMCBase = basename($sSMCFile, '.gz');

		do {
			$sLn = fgets($h);
			if ($sLn === false || $nPartLen >= $nBlockSize) {
				if ($sContent != '') {
					fnWriteSplitPart($sChrOutDir, $nPart++, $sHeader, $sContent, $sSMCBase);
				}

				$nPartLen = 0;
				$sContent = "";
				if ($sLn === false) {
					exec("echo $nTotalLen > $sChrOutDir/length.txt");
					break;
				} else {
					continue;
				}
			}

			$sLn = trim($sLn);

			if ($sLn == '') continue;

			if ($sLn[0] == '#') {
				$sHeader = $sLn;
				continue;
			}

			$arrF = explode(' ' , $sLn);
			$nPartLen += intval($arrF[0]);
			$nTotalLen += intval($arrF[0]);
			$sContent .= "$sLn\n";

		} while (true);
	}

	exec("touch $bSplitDoneFlag");

}

//start bootstraping.
$arrSMCStems = fnGetSMCStems();

for($nRep=0;$nRep<$nReps;$nRep++) {
	$sBSOut = "$sOutDir/bootstrap/$nRep";
	exec("mkdir -p $sBSOut");

	$arrLens = glob("$sSplitDIR/*/length.txt");

	//die();
	foreach($arrLens as $sChrLenFile) {
		$sChr = basename(dirname($sChrLenFile));
		$nChrLen = intval(file_get_contents($sChrLenFile));
		$sChrOutDir = $sBSOut."/$sChr";
		exec("mkdir -p $sChrOutDir");

		//foreach($arrSMCStems as $sSMCStem) {
		$arrParts = fnLoadParts($arrLens, $arrSMCStems);
		//print_r($arrParts);

		$arrHOUT = array();
		$arrWrittenHeader = array();		

		foreach($arrSMCStems as $nIdx => $sSMCStem) {
			$arrHOUT[$nIdx] = popen("gzip -c > $sChrOutDir/$sSMCStem.gz", 'w');
			$arrWrittenHeader[$nIdx] = false;
			exec("echo components > $sChrOutDir/$sSMCStem.log.txt");
		}

		$nWrittenLen = 0;

		$nTargetLen = $nChrLen * count($arrSMCStems) ;
		do {

			$arrNewPart = fnDrawPart($arrParts);
			foreach($arrSMCStems as $nIdx => $sSMCStem) {
				list($sNewPart, $hNewPart) = $arrNewPart[$nIdx];
				$hOut = $arrHOUT[$nIdx];
				exec("echo $sNewPart >> $sChrOutDir/$sSMCStem.log.txt");
				$bDone = false;

				do {
					$sLn = fgets($hNewPart); 

					if ($sLn === false || $nWrittenLen >= $nTargetLen ) {
						$bDone = ($nWrittenLen >= $nTargetLen);
						break;
					}

					if (trim($sLn) == '') continue;

					if ($sLn[0] == '#' ) {
						if (!$arrWrittenHeader[$nIdx] ) {
							fwrite($hOut , $sLn);
							$arrWrittenHeader[$nIdx]  = true;
						}

						continue;
					}

					$arrF = explode(' ' , $sLn);
					$nWrittenLen += intval($arrF[0]);
					fwrite($hOut , $sLn);
				
				} while(true);

				pclose($hNewPart);
			} //foreach
		} while(!$bDone);

		foreach($arrSMCStems as $nIdx => $sSMCStem) {
			pclose($arrHOUT[$nIdx]);
			
		}


		//}
	}
}

function fnWriteSplitPart($sChrOutDir, $nPart, $sHeader, $sContent, $sSMCBase) {
	$sOutFile = "$sChrOutDir/$sSMCBase.part.$nPart.gz";
	$hOut = popen("gzip -c > $sOutFile", 'w');
	fwrite($hOut, $sHeader."\n");
	fwrite($hOut, $sContent);
	pclose($hOut);
}

function fnLoadParts($arrLens, $arrSMCStems) {
	echo("Loading parts...\n");
	$arrRet = array();

	foreach($arrSMCStems as $nIdx => $sSMCStem) {
		$arrRet[$nIdx] = array();
	}

	foreach($arrLens as $sChrLenFile) {
		$sChrDir = (dirname($sChrLenFile));
		foreach($arrSMCStems as $nIdx => $sSMCStem) {
			$arrRet[$nIdx] = array_merge($arrRet[$nIdx], glob($sChrDir."/$sSMCStem.part*.gz"));
		}
	}

	return $arrRet;
}

function fnDrawPart($arrPartsARR) {

	$nMax = count($arrPartsARR[0])-1;
	$nDrawn = random_int(0, $nMax);
	$arrRet = array();

	foreach($arrPartsARR as $nStemIdx => $arrParts) {
		$arrRet[$nStemIdx] = array($arrParts[$nDrawn], popen( "zcat -f ". $arrParts[$nDrawn], 'r'));
	}

	return $arrRet;
}

function fnGetSMCStems() {
	global $sRawDir;
	$arrSMCFiles = glob("$sRawDir/*/smc*.gz");
	$arrStems = array();
	foreach($arrSMCFiles as $sSMCF) {
		$sStem = basename($sSMCF, '.gz');
		$arrStems[$sStem] = true;
	}

	return array_keys($arrStems);
} 

?>
