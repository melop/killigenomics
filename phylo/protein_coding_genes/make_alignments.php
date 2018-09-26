<?php
require_once(dirname(__FILE__) . "/lib.php");

//This script reads in the fasta files from the CodeML fasta files
//it concatenates all protein-coding genes into one big alignment.

$sMitosOut = "./lnkfasta";


$sOut = "aln";

$arrMitosOut = glob("$sMitosOut/*.fas");

$arrGeneAlns = array(); //key is gene name, value is array( sp => seq )

exec("mkdir -p $sOut");
$arrSpecies = array();

foreach($arrMitosOut as $sMitoOut) {
	$sFileName = basename( $sMitoOut );
	list($sGene) = explode(".fas" , $sFileName);

	$arrSeqs = fnParseFasta($sMitoOut);
	//print_r($arrSeqs);

	//die();

	foreach($arrSeqs as $sSpName => $sSeq) {

		if (!array_key_exists($sGene, $arrGeneAlns)) {
			$arrGeneAlns[$sGene] = array();
		}

		$arrGeneAlns[$sGene][$sSpName] = $sSeq;
		if (!array_key_exists( $sSpName ,$arrSpecies) ) {
			$arrSpecies[$sSpName] = true;
		}

	}

}

$arrSpecies = array_keys($arrSpecies);
$arrAAConcat = array();
//perform TranslatorX alignment:
foreach($arrGeneAlns as $sGene => $arrAln) {

	$arrIncludedSpecies = array_keys($arrAln);
	$nLen = strlen($arrAln[$arrIncludedSpecies[0]]);

	if ($nLen % 3 != 0 ) {
		echo("Discard $sGene, length not multiplication of 3!\n");
		continue;
	}

	foreach($arrAln as $sSp => $sSeq) {
		if (!array_key_exists($sSp, $arrAAConcat)) {
			$arrAAConcat[$sSp] = "";
		}	
			$arrAAConcat[$sSp] .= $sSeq;
	}

	foreach($arrSpecies as $sSp) {

		if (!array_key_exists($sSp, $arrAAConcat)) {
			$arrAAConcat[$sSp] = "";
		}	
		if (!array_key_exists($sSp , $arrAln)) {
			$arrAAConcat[$sSp] .= str_repeat('-' , $nLen);
		}
	}

}

$hSumAA = fopen($sOut."/all.aln.fas" , 'w');

foreach($arrAAConcat as $sSp => $sSeq) {
	fwrite($hSumAA , ">$sSp\n". wordwrap( $sSeq , 100,  "\n" , true )."\n");
}

function fnParseFasta($sMitoOut) {
	echo("Reading $sMitoOut...\n");
	$hMtGenes = fopen($sMitoOut , 'r');
	$sSeq = "";
	$sSeqName = "";
	$arrRet = array(); //key is gene name
	do {
		$sLn = fgets($hMtGenes );
		if ($sLn === false || $sLn[0] == '>') {
			if ($sSeqName != '') {
				$arrRet[$sSeqName] = $sSeq;
			}
			if ($sLn === false) {
				break;
			}

			$sSeq = "";
			$sSeqName = trim(substr($sLn , 1));
			continue;
		}

		$sSeq .= trim($sLn);
	} while (true);

	return $arrRet;
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
			$sSeqName = trim(substr($sLn , 1));
			continue;
		}

		$sSeq .= trim($sLn);
	} while (true);

	return $arrRet;
	
}


?>
