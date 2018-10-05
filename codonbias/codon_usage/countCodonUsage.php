<?php
//This script counts codon usage per gene for the specified species
$nShiftFrame = 0; //how many frames to shift, this is used for generating null datasets
$bRevComp = false; //rev comp the cds, this is used for null dataset generation

$sPremadeFasta = "/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/export_protein_multiref_23_july_2017/ret_improved_cleanDNA.fasta";

$hPremadeFasta = fopen($sPremadeFasta , "r");

$arrCountSpecies = array(
"AAU", 
"ACG", 
"ACL", 
"ACM", 
"ACY", 
"AGM", 
"AKM", 
"AKN1", 
"AKN2", 
"APC", 
"AMS", 
"APT", 
"ARG", 
"CHW", 
"CMR", 
"CSB", 
"CTO", 
"EDG", 
"EGH", 
"EGN", 
"ELM", 
"EMF", 
"ESP", 
"ETO", 
"SBT", 
"SCV", 
"SGG", 
"SSC", 
"FAM", 
"FDL", 
"FFL", 
"FSC", 
"FSJ", 
"FGR", 
"FTH", 
"NFS", 
"NKF", 
"NOC", 
"NOR", 
"NRC", 
"NVG", 
"NVS", 
"NFZ", 
"PSY", 
"APL", 
"PLP");

$sLn = "";
$nShiftFrame = 0;
$bRevComp = false;

$sOut = "codon_use_counts";

if ($nShiftFrame != 0) {
	$sOut .= "_frameshift_".$nShiftFrame;
} 

if ($bRevComp) {
	$sOut .= "_revcomp";
}

$sOut .= ".txt";

$hOut = fopen($sOut , "w");

$sSeqName = "";
$sSeq = "";

fwrite($hOut , "GeneName\tCodon");
foreach($arrCountSpecies as $sSpecies) {
	fwrite($hOut , "\t$sSpecies");
}
fwrite($hOut , PHP_EOL);


do {
	$sLn = fgets($hPremadeFasta);

	if ($sLn==false) {
		fnProcessRecord($sSeqName, $sSeq);
		break;
	}

	$sLn = trim($sLn);
	if (substr($sLn , 0,1) == ">") {
		fnProcessRecord($sSeqName, $sSeq);
		$sSeqName = substr($sLn , 1);
		$sSeq = "";
		
	} else {
		$sSeq .= $sLn;
	}

} while (true);
fnProcessAln($sCurrentGene,$arrFullCodingSeq); //don't forget the last one

$arrFullCodingSeq = array();
$sCurrentGene = "";
function fnProcessRecord($sSeqName, $sSeq) {
	global $arrFullCodingSeq , $sCurrentGene;

	if ($sSeqName == '') return;

	$arrF1 = explode(';',$sSeqName);
	$arrParams = array();
	foreach($arrF1 as $s) {
		$arrKV = explode(':' , $s);
		if (count($arrKV) == 2) {
			$arrParams[$arrKV[0]] = $arrKV[1];
		}
	}
	/*echo($sSeqName);
	print_r($arrParams);
	die();*/
	$sGeneName = $arrParams['GroupId'];
	$sSpecies = $arrParams['Sp'];

	if ($sCurrentGene == $sGeneName ) {
		if (!array_key_exists($sSpecies , $arrFullCodingSeq)) {
			$arrFullCodingSeq[$sSpecies] = "";
		}

		$arrFullCodingSeq[$sSpecies] .= $sSeq;
	} else {
		if (count($arrFullCodingSeq) > 0) {
			fnProcessAln($sCurrentGene,$arrFullCodingSeq);
		}
		$arrFullCodingSeq = array();
		$sCurrentGene = $sGeneName;
		$arrFullCodingSeq[$sSpecies] = $sSeq;

	}



}



function fnProcessAln($sCurrentGene,$arrFullCodingSeq) {
	global $hOut, $arrCountSpecies;
	echo("Doing $sCurrentGene ...\n");
	//print_r($arrFullCodingSeq);
	$arrCountRet = array();
	$arrCodonStates = array();
	foreach($arrCountSpecies as $sSpecies ) {
		$sSeq = "";
		if (!array_key_exists($sSpecies , $arrFullCodingSeq)) {
			//this species is missing, no data
		} else {
			$sSeq = $arrFullCodingSeq[$sSpecies];
		}

		$arrCountRet[$sSpecies]  = fnCountCodons($sSeq);
		$arrCodonStates = array_keys($arrCountRet[$sSpecies]);
	
	}

	foreach($arrCodonStates as $sCodon) {
		fwrite($hOut , $sCurrentGene."\t".$sCodon);

		foreach($arrCountRet as $sSpecies => $arrCounts) {
			fwrite($hOut , "\t".$arrCounts[$sCodon]);
		}
		fwrite($hOut , PHP_EOL);
	}
}

function fnCountCodons($sSeq) {
	global $nShiftFrame , $bRevComp; 

	if ($bRevComp) {
		$sSeq = fnReverseComplement($sSeq);
	}

	if ($nShiftFrame > 0) {
		$sSeq = str_repeat("N" , $nShiftFrame ) . substr($sSeq, 0, strlen($sSeq) - $nShiftFrame ); //shift the sequence
	}

	$arrRet = fnInitCountArray();

	if (strlen($sSeq) % 3 !=0 ) {
		echo($sSeq);
		die("Length of seq cannot be divided by 3\n");
	}
	for($i=0;$i<=strlen($sSeq)-3;$i+=3) {
		$sCodon = strtoupper(substr($sSeq, $i, 3));
		if (strpos($sCodon , "N") !==false || strpos($sCodon , "-") !==false) {
			continue;
		}

		if (!array_key_exists($sCodon, $arrRet) ) {
			die($sCodon." is a codon unexpected.\n");
		}

		$arrRet[$sCodon] += 1;
	}

	return $arrRet;
}

function fnInitCountArray() {
	$arrBases = array('A','T','G','C');
	$arrRet = array();

	foreach($arrBases as $sBase1) {
	foreach($arrBases as $sBase2) {
	foreach($arrBases as $sBase3) {
		$arrRet[$sBase1.$sBase2.$sBase3] = 0;
	}
	}
	}

	return $arrRet;
}

function fnReverseComplement($sIn) {


	return strrev( strtr($sIn, 'ACBDGHKMNSRUTWVYacbdghkmnsrutwvy', 'TGVHCDMKNSYAAWBRTGVHCDMKNSYAAWBR'));
	
	
}


?>
