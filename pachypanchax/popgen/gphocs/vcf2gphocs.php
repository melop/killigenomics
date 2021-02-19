<?php
$sVCF = "../PLP_merged.vcf.gz";

$arrExclude = array("PML");
$nBlockSize = 1000;//5000;
$nBreakBlockDist = 20; //if this many bases are excluded, break the block.
$nMinDist = 50000; //100000;
$nMinCov = 196;
$nMaxCov = 572;
$nMinQual = 285;
$nTotalBlocks = 0;

$sOut = "PLP.gphocs.$nBlockSize.$nMinDist.in";

$hIn = popen("zcat -f $sVCF", 'r');
$hOut = fopen($sOut, 'w');

$arrExclude = array_flip($arrExclude);
$bHeaderRead = false;

$nBlockTotalCov = 0;
$arrBlockSeq = array();
$arrBlockSeqTemplate = array();
$nBlockStart = 0;
$nBlockEnd = 0;
$nBlockLen = 0;
$nPrevGoodPos = 0;
$sPrevScf = '';
$nSkipUntil = 0;
$arrIUPACMAP = array('CT' => 'Y', 'TC' => 'Y', 'AG' => 'R', 'GA' => 'R', 'AT' => 'W', 'TA' => 'W' , 'CG' => 'S', 'GC' => 'S', 'GT' => 'K', 'TG' => 'K', 'AC' => 'M', 'CA' => 'M');

while(false !== ($sLn = fgets($hIn)) ) {
	$sLn = trim($sLn);
	$arrF = explode("\t", $sLn);
	if (substr($sLn, 0, 6) == '#CHROM' ) {
		$arrColNames = array_flip($arrF);
		$arrSampleNames = array_slice($arrColNames, 9);
		foreach($arrSampleNames as $sName => $nIdx) {
			if (array_key_exists($sName, $arrExclude) ) {
				unset($arrSampleNames[$sName]);
			}
		}
		
		print_r($arrSampleNames);
		die();
		$bHeaderRead = true;
		$arrBlockSeqTemplate = array_fill_keys(array_keys($arrSampleNames), '');
		continue;
	}
	
	if (!$bHeaderRead) {
		continue;
	}
	
	$sScf = $arrF[0];
	$nPos = $arrF[1];
	$sRefBase = $arrF[3];
	$sAltBase = $arrF[4];
	$nQual = $arrF[5];
	$sInfo = $arrF[7];
	
	if ($sScf == $sPrevScf && $nSkipUntil > $nPos) {
		continue; //skipping
	}
	if (strlen($sRefBase) != 1 || strlen($sAltBase)!=1) {//doesn't break block
		continue;
	}
	
	$oInfo = fnParseInfo($sInfo);
	
	//print_r($oInfo);
	if (!array_key_exists('DP', $oInfo) || $oInfo['DP'] > $nMaxCov || $oInfo['DP'] < $nMinCov || $nQual < $nMinQual) {//doesn't break block
	
		continue;
	}
	
	if ($sScf != $sPrevScf || ($nPos - $nPrevGoodPos )>$nBreakBlockDist ) {//break and discard previous block
		$arrBlockSeq = $arrBlockSeqTemplate;
		$nBlockStart = $nPos;
		$sPrevScf = $sScf;
		$nBlockLen=0;
		$nSkipUntil = 0;
	}
	
	
	foreach($arrSampleNames as $sTax => $nIdx) {
		$arrBlockSeq[$sTax] .= fnGenotype2Base( $sRefBase, $sAltBase, $arrF[$nIdx]);
	} 
	
	$nBlockLen++;
	$nPrevGoodPos = $nPos;
	
	if ($nBlockLen >= $nBlockSize) {
		$nBlockEnd = $nPos;
		fnWriteBlock("$sScf.$nBlockStart.$nBlockEnd", $arrBlockSeq);
		$nSkipUntil = $nBlockEnd + $nMinDist;
		$arrBlockSeq = array();
		$nBlockStart = 0;
		$nBlockEnd = 0;
		$nBlockLen = 0;
		$nPrevGoodPos = 0;

	}
}

fclose($hOut);
file_put_contents("$sOut.tmp", $nTotalBlocks."\n\n");
exec("cat $sOut >> $sOut.tmp; mv $sOut.tmp $sOut");

function fnWriteBlock($sBlockName, $arrBlockSeq) {
	global $hOut, $nTotalBlocks;
	$nSeq = count($arrBlockSeq);
	$arrTax = array_keys($arrBlockSeq);
	$nTax = count($arrTax);
	$nLen = strlen($arrBlockSeq[$arrTax[0]]);
	echo("$sBlockName\n");
	//print_r($arrBlockSeq);
	//die();
	fwrite($hOut, "$sBlockName $nTax $nLen\n");
	foreach($arrTax as $sTax) {
		fwrite($hOut, "$sTax\t".$arrBlockSeq[$sTax]."\n");
	}
	fwrite($hOut, "\n");
	$nTotalBlocks++;
}

function fnGenotype2Base( $sRefBase, $sAltBase, $sGeno) {
	global $arrIUPACMAP;
	
	if ($sGeno[0] == '.' || $sGeno[2] == '.') {
		return 'N';
	}
	
	$sHetBase = $sRefBase;
	if ($sAltBase  == '.') {
		$sAltBase = $sRefBase;
	} else {
		$sCombBase = $sRefBase.$sAltBase;
		if (!array_key_exists($sCombBase, $arrIUPACMAP)) {
			die("Error $sCombBase not found in IUPAC TABLE\n");
		}
		$sHetBase = $arrIUPACMAP[$sRefBase.$sAltBase];
	}
	
	$nGeno = intval($sGeno[0]) + intval($sGeno[2]);
	if ($nGeno == 0) {
		return $sRefBase;
	} else if ($nGeno == 2) {
		return $sAltBase;
	} else if ($nGeno == 1) {
		return $sHetBase;
	}
	
	die("Genotype error: $nGeno\n");
}


function fnParseInfo($sInfo) {
	$arrF = explode(';', $sInfo) ;
	$arrRet = array();
	foreach($arrF as $sF) {
		list($sN, $sV) = explode('=', $sF);
		$arrRet[$sN] = $sV;
	}
	return $arrRet;
}




?>
