<?php
//concatenate cleaned protein alignments into one big alignment for RAxML

$sDIR = "../FullTaxaOrthologs";
$arrFiles = glob($sDIR."/Group_*_clean.fa");
$sOut = "concat.aln.phy";

$hOut = fopen($sOut, "w");

$arrFullAln = array();
$nTotalAlnLen = 0;

foreach($arrFiles as $sFile) {
	$arrAln = fnReadFasta($sFile);
	$arrVal = array_values($arrAln); 
	$nAlnLen = strlen($arrVal[0]);
	//print_r($arrAln);
	//break;

	if (count($arrFullAln) == 0) {
		$arrTaxaName = fnToTaxaName(array_keys($arrAln) );
		$arrFullAln = array_combine($arrTaxaName , array_fill( 0 , count($arrTaxaName) , "" ) );
	}

	//print($arrFullAln);
	$nTotalAlnLen += $nAlnLen;
	$arrAlnConverted = array_combine( array_keys($arrFullAln) , array_fill( 0 , count($arrFullAln) , str_repeat('-', $nAlnLen) ) ); //key is taxon.
	$arrAddedTaxa = array();

	foreach($arrAln as $sSeqName => $sSeq) {
		$arrSeqName = explode('|' , $sSeqName);
		$sTaxon = $arrSeqName[0];
		if (!array_key_exists($sTaxon, $arrAlnConverted)) {
			die("Taxon $sTaxon not found!\n");
		}
		if (array_key_exists($sTaxon,$arrAddedTaxa) ) {
			die("Taxon $sTaxon is repeated!\n");
		}
		$arrAddedTaxa[$sTaxon] = true;
		if (strlen($sSeq) != $nAlnLen ) {
			die("Sequence length differ!\n");
		}
		$arrAlnConverted[$sTaxon] = $sSeq;
	}

	foreach($arrFullAln as $sTaxonName => $sSeq) {
		$arrFullAln[$sTaxonName] .= $arrAlnConverted[$sTaxonName];
	}
}
echo("Total alignment length: $nTotalAlnLen\n");
$nTotalAlnLen = 0;

	foreach($arrFullAln as $sTaxonName => $sSeq) {
		if ($nTotalAlnLen == 0) {
			$nTotalAlnLen = strlen($sSeq);
			fwrite($hOut, " ".count($arrFullAln)." ".$nTotalAlnLen."\n");
		}


		fwrite($hOut, "$sTaxonName\t\t" . $sSeq . "\n" );
	}

echo("Total alignment length: $nTotalAlnLen\n");

function fnReadFasta($sFasta) {
	$hF = fopen($sFasta , "r");
	$arrRet = array();

	$sName = "";
	$sSeq = "";
	do {
		$sLn = fgets($hF); 
		if (false===$sLn || ($sLn[0] == '>' && $sName!='') ) {
			$arrRet[$sName] = $sSeq;
			if (false===$sLn) break;
		}
		$sLn = trim($sLn);
		if ($sLn[0] == '>') {
			$sName = substr($sLn, 1);
			$sSeq = "";
			continue;
		}
		$sSeq .= $sLn;
		
	} while(true);

	return $arrRet;
}

function fnToTaxaName($arrSeqName) {
	$arrRet = array();
	foreach($arrSeqName as $sSeqName) {
		$arrF = explode('|' , $sSeqName);
		$arrRet[] = $arrF[0];
	}

	return $arrRet;
}

?>
