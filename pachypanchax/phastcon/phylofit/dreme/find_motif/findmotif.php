<?php
$sIn = "../plp.nfzref.acc.nonCDS.norepeat.fa";
$sPattern = $sPattern1 = "TAATTA";
$sPattern = $sPattern1 = "CAGCTG";

$sOutDetail = "$sPattern.detail.gff";
$sOutBrief = "$sPattern.brief.gff";

$h1 = fopen($sOutDetail, 'w'); 
$h2 = fopen($sOutBrief, 'w'); 

$bReverseComplement = true;
$arrFasta = fnParseFasta($sIn);

$arrIUPAC = array("A" => "A", "G" => "G" , "T" => "T", "C" => "C",
"Y" => "[CT]", 
"R" => "[AG]", 
"W" => "[AT]", 
"S" => "[GC]", 
"K" => "[TG]", 
"M" => "[CA]", 
"D" => "[AGT]", 
"V" => "[ACG]", 
"H" => "[ACT]", 
"B" => "[CGT]", 
"N" => "[ATCG]" );

if (fnRevComp($sPattern) == $sPattern && $bReverseComplement) {
	echo("Warning: $sPattern is palindromic, the same match will be reported twice if you turn on the reverse complement option.\n");
}

$nPrimerLen = strlen($sPattern);
$sPattern = str_replace(array_keys($arrIUPAC), array_values($arrIUPAC) , $sPattern);


foreach($arrFasta as $sName => $sSeq) {
	list($sType, $sChr, $nLeftCoord, $nRightCoord) = explode('|', $sName);
	$arrRet = array();
	 fnGetAllPos($arrRet, $sPattern, $sSeq, $nLeftCoord, $nRightCoord, false);
	if ($bReverseComplement) {
	 fnGetAllPos($arrRet, $sPattern, fnRevComp($sSeq), $nLeftCoord, $nRightCoord, true);
	}


	foreach($arrRet as $nMatchIdx => $oRet) {
		$nStart = $nEnd = 0;
		$sStrand = '+';
		if ($oRet[2]) { // is rev strand
			$sStrand = '-';
			$nEnd = $oRet[1];
			$nStart = $nEnd - ($nPrimerLen-1);
		} else {
			$sStrand = '+';
			$nStart = $oRet[1];
			$nEnd = $nStart + ($nPrimerLen - 1);
		}

		$s = "$sChr\tDREME\t$sType.$sPattern1\t$nStart\t$nEnd\t$sStrand\t.\tParent=$sName;ID=$sName.$nMatchIdx;Match=$oRet[0]\n";
		fwrite($h1, $s);
	}

	$nTotalMatches = count($arrRet);
	if ($nTotalMatches>0) {
		$s = "$sChr\tDREME\t$sType.$sPattern1\t$nStart\t$nEnd\t$sStrand\t.\tParent=$sName;ID=$sName.0;MatchCount=$nTotalMatches\n";
		fwrite($h2, $s);
	}
}

function fnParseFasta($sFas) {
	$arrRet = array();
	$sSeq = '';
	$sName = '';
	$h = fopen($sFas, 'r');
	while(true) {
		$sLn = fgets($h);
		if ($sLn === false || (trim($sLn)!='' && $sLn[0] == '>') ) {
			if ($sName != '') {
				if (array_key_exists($sName, $arrRet) ) {
					die("$sName is repeated in fasta\n");
				}
				$arrRet[$sName] = $sSeq;
				$sSeq = '';
				$sName = '';
			}

			if ($sLn === false) break;
			$sName = substr(trim($sLn),1);
			continue;
		}

		$sSeq .= strtoupper(trim($sLn));
	}

	return $arrRet;
}

function fnRevComp($s) {
        return strrev( strtr($s, 'ACBDGHKMNSRUTWVYacbdghkmnsrutwvy', 'TGVHCDMKNSYAAWBRTGVHCDMKNSYAAWBR'));

}

function fnGetAllPos(&$arrOffSets, $sPattern, $sSeq, $nLeftCoord, $nRightCoord, $bIsRevStrand) {
	preg_match_all("/$sPattern/", $sSeq, $arrRet, PREG_OFFSET_CAPTURE);

	foreach($arrRet[0] as $arrCapture) {
		$nOffSet = $arrCapture[1];
		if (!$bIsRevStrand) {
			$nOffSet = $nLeftCoord + $nOffSet;
		} else {
			$nOffSet = $nRightCoord - $nOffSet;
		}

		$arrOffSets[] = array($arrCapture[0], $nOffSet, $bIsRevStrand);
	}

}

?>
