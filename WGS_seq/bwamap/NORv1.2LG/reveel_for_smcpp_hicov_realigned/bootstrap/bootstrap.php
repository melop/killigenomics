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

	$arrSMC = glob("$sRawDir/*/smc.gz");

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

		do {
			$sLn = fgets($h);
			if ($sLn === false || $nPartLen >= $nBlockSize) {
				if ($sContent != '') {
					fnWriteSplitPart($sChrOutDir, $nPart++, $sHeader, $sContent);
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

for($nRep=0;$nRep<$nReps;$nRep++) {
	$sBSOut = "$sOutDir/bootstrap/$nRep";
	exec("mkdir -p $sBSOut");

	$arrLens = glob("$sSplitDIR/*/length.txt");
	$arrParts = fnLoadParts($arrLens);
	print_r($arrParts);

	//die();
	foreach($arrLens as $sChrLenFile) {
		$sChr = basename(dirname($sChrLenFile));
		$nChrLen = intval(file_get_contents($sChrLenFile));
		$sChrOutDir = $sBSOut."/$sChr";
		exec("mkdir -p $sChrOutDir");
		$hOut = popen("gzip -c > $sChrOutDir/smc.gz", 'w');
		exec("echo components > $sChrOutDir/log.txt");
		$nWrittenLen = 0;
		$bWrittenHeader = false;
		do {
			list($sNewPart, $hNewPart) = fnDrawPart($arrParts);
			exec("echo $sNewPart >> $sChrOutDir/log.txt");
			$bDone = false;

			do {
				$sLn = fgets($hNewPart); 

				if ($sLn === false || $nWrittenLen >= $nChrLen ) {
					$bDone = ($nWrittenLen >= $nChrLen);
					break;
				}

				if (trim($sLn) == '') continue;

				if ($sLn[0] == '#' ) {
					if (!$bWrittenHeader) {
						fwrite($hOut , $sLn);
						$bWrittenHeader = true;
					}

					continue;
				}

				$arrF = explode(' ' , $sLn);
				$nWrittenLen += intval($arrF[0]);
				fwrite($hOut , $sLn);
				
			} while(true);

			pclose($hNewPart);

		} while(!$bDone);
	}
}

function fnWriteSplitPart($sChrOutDir, $nPart, $sHeader, $sContent) {
	$sOutFile = "$sChrOutDir/smc.part.$nPart.gz";
	$hOut = popen("gzip -c > $sOutFile", 'w');
	fwrite($hOut, $sHeader."\n");
	fwrite($hOut, $sContent);
	pclose($hOut);
}

function fnLoadParts($arrLens) {
	echo("Loading parts...\n");
	$arrRet = array();
	foreach($arrLens as $sChrLenFile) {
		$sChrDir = (dirname($sChrLenFile));
		$arrRet = array_merge($arrRet, glob($sChrDir."/smc.part*.gz"));
	}

	return $arrRet;
}

function fnDrawPart($arrParts) {
	$nMax = count($arrParts)-1;
	$nDrawn = random_int(0, $nMax);
	return array($arrParts[$nDrawn], popen( "zcat -f ". $arrParts[$nDrawn], 'r'));
}

?>
