<?php
require_once(dirname(__FILE__) . "/lib.php");

//This script reads in the fasta file from the MITOS output
//and makes one alignment file per protein-coding gene (trimmed with TranslatorX)
//it also concatenates all protein-coding genes into one big alignment.
//if there are multiple copies of a gene in a species (suffixes _a, _b) etc present
//the species is going to be excluded.
//if you don't want that, inspect the files mannually, and decide which one to use
//and delete the suffix for that gene.
//often this is due to a frameshift in the mt assembly. you can also fix that and rerun MITOS.

$sMitosOut = "../out_adjusted";


$sOut = "aln";

$arrMitosOut = glob("$sMitosOut/*.fasta");

$arrGeneAlns = array(); //key is gene name, value is array( sp => seq )

exec("mkdir -p $sOut");
$arrSpecies = array();

foreach($arrMitosOut as $sMitoOut) {
	$sFileName = basename( $sMitoOut );
	list($sSpName) = explode(".fasta" , $sFileName);
	$arrSpecies[] = $sSpName;

	$arrSeqs = fnParseFasta($sMitoOut);
	//print_r($arrSeqs);

	//die();

	foreach($arrSeqs as $sGene => $sSeq) {
		if (substr($sGene, 0, 3) == 'trn' ) {
			continue;
		}

		if (!array_key_exists($sGene, $arrGeneAlns)) {
			$arrGeneAlns[$sGene] = array();
		}

		$arrGeneAlns[$sGene][$sSpName] = $sSeq;
	}

}

$arrAAConcat = array();
//perform TranslatorX alignment:
foreach($arrGeneAlns as $sGene => $arrSeqs) {
	$arrFileHandles = array();

	$arrFileHandles[0] = fopen($sOut."/$sGene.raw.fas" , 'w');
	$arrFileHandles[1]  = fopen($sOut."/$sGene.aln.fas" , 'w');

	if (substr($sGene , 0, 3) != 'rrn') {

		$arrFileHandles[2] = fopen($sOut."/$sGene.AA.raw.fas" , 'w');
		$arrFileHandles[3] = fopen($sOut."/$sGene.AA.aln.fas" , 'w');
	} else {
		//foreach($arrSeqs as $sSp => $sSeq) {
		//	fwrite($arrFileHandles[0] , ">$sSp\n$sSeq\n");
		//}
		//continue;//do not align the rRNA
	}

	$arrTranslatorX = fnTranslatorX(fnExcludeLastStop($arrSeqs), 0.75, 0.85, 2, 5, 0.1, 2);
	//print_r($arrTranslatorX);
	//die();
	if ($arrTranslatorX === false) {
		echo("TranslatorX failed:  $sGene\n");
		continue;
	}

	foreach($arrTranslatorX as $nIdx => $sAln) {
		if (substr($sGene , 0, 3) == 'rrn' && $nIdx > 1) {
			continue;
		} 

		if (substr($sGene , 0, 3) != 'rrn' &&  strpos($sGene , '_') == false &&  strpos($sGene , '-') ==false && $nIdx == 1) {
			$arrAln = fnFasta2Array($sAln);
			$arrIncludedSpecies = array_keys($arrAln);
			$nLen = strlen($arrAln[$arrIncludedSpecies[0]]);
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

		fwrite($arrFileHandles[$nIdx] , $sAln);
		
	}

}

$hSumAA = fopen($sOut."/all.aln.fas" , 'w');

foreach($arrAAConcat as $sSp => $sSeq) {
	fwrite($hSumAA , ">$sSp\n". wordwrap( $sSeq , 100,  "\n" , true )."\n");
}

function fnParseFasta($sMitoOut) {
	echo("Reading $sMitoOut...\n");
	$hMtGenes = fopen($sMitoOut."/ret.fas" , 'r');
	$sSeq = "";
	$sSeqName = "";
	$arrRet = array(); //key is gene name
	do {
		$sLn = fgets($hMtGenes );
		if ($sLn === false || $sLn[0] == '>') {
			if ($sSeqName != '') {
				$arrSeqName = explode(';' , $sSeqName);
				if (count($arrSeqName)!= 4) {
					die("Sequence name $sSeqName is not in MITOS output format ($sMitoOut).\n");
				}
				$sGeneName = trim($arrSeqName[3]);
				$arrRet[$sGeneName] = $sSeq;
			}
			if ($sLn === false) {
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
