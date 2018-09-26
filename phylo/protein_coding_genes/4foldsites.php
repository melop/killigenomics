<?php
require_once(dirname(__FILE__) . "/lib.php");

$sIn = "aln/all.aln.fas";
$sOut = "aln/all.4fold.aln.phy";

$arrIn = fnFasta2Array(file_get_contents($sIn) );

//check len
$nLen = -1;
$arr4FoldSites = array();
foreach($arrIn as $sSp => &$sSeq) {
	$nThisLen = strlen($sSeq);
	$arr4FoldSites[$sSp] = "";
	if ($nLen == -1) {
		$nLen = $nThisLen;
		if ($nLen % 3 != 0) {
			die("Error, sequence length not a multipication of 3!\n");
		}
		continue;
	}

	if ($nLen != $nThisLen) {
		die("Seq len not equal\n");
	}
}



for($nPos=0;$nPos <= ($nLen-3);$nPos+=3 ) {
	//first check to make sure that all taxa have the same AA here. if not, skip
	$sAA = "";
	$sCodon = "";
	$arr4FoldToAppend = array();
	foreach($arrIn as $sSp => &$sSeq) {
		list($sThisAA, $sThisCodon) = fnDNA2Peptide(substr($sSeq , $nPos  , 3), false, false);

		$arr4FoldToAppend[$sSp] = $sThisCodon[2];

		if ($sCodon == '') {
			$sCodon = $sThisCodon;
			$sAA = $sThisAA;
			if (fnIs4FoldDegen($sCodon) !== true) {
				continue 2;
			}

			continue;
		}

		if ($sThisAA != $sAA || substr($sCodon, 0,2) != substr($sThisCodon, 0,2)  ) {
			continue 2; 
		}
	}

	echo("Position: " . ($nPos+1+2). " Codon: $sCodon; AA: $sAA\n" );

	foreach($arr4FoldToAppend as $sSp => $sBase) {
		$arr4FoldSites[$sSp] .= $sBase;
	}
	
}

$hOut = fopen($sOut , 'w');
$bHeaderWritten = false;
foreach($arr4FoldSites as $sSp => $sSeq) {
	if (!$bHeaderWritten) {
		fwrite($hOut , " ".count($arr4FoldSites ) . " ".strlen($sSeq)."\n");
		$bHeaderWritten = true;
	}
	fwrite($hOut , "$sSp\t\t$sSeq\n");
}

function fnFasta2Array($str) {

	$arrLns = explode("\n", $str);

	$sSeq = "";
	$sSeqName = "";
	$arrRet = array(); //key is gene name
	do {
		$sLn = current($arrLns );
		next($arrLns);
		if ( (!$sLn) || $sLn[0] == '>') {
			if ($sSeqName != '') {

				$arrRet[$sSeqName] = $sSeq;
			}
			if ( !$sLn ) {
				break;
			}

			$sSeq = "";
			$sSeqName = substr($sLn , 1);
			continue;
		}

		$sSeq .= trim($sLn);
	} while (true);

	return $arrRet;
	
}


?>
