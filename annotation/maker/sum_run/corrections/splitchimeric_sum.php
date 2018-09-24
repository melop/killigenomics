<?php
require_once(dirname(__FILE__) . "/lib.php");
$sMakerGFF = "split__main.gff";
$sWorkDir = "./splitchimeric_improved/";
$sGenome = "query.fa";
$sAppendGFF = "_append.gff";
$sOutGFF = "chimericfixed.gff";

while(count($argv) > 0) {
    $arg = array_shift($argv);
    switch($arg) {
        case '-R':
            $sGenome = trim(array_shift($argv));
            break;
	case '-M':
	    $sMakerGFF  = trim(array_shift($argv));
	    break;
	case '-o':
	    $sWorkDir  = trim(array_shift($argv));
	    break;

    }
}

if (!file_exists($sWorkDir)) {
	die("Cannot find $sWorkDir\n");
}

if (!file_exists($sMakerGFF)) {
	die("Cannot find $sWorkDir\n");
}

$hAppendGFF = fopen($sAppendGFF , "w");
$hOutGFF = fopen($sOutGFF , "w");

$oMakerGFF = new MappedCDSGFF();


$oMakerGFF->LoadGFF($sMakerGFF , $sGenome);


$arrGenes = $oMakerGFF->GetGenes('maker'); // get all maker genes

$arrDeleteGenes = array();
foreach($arrGenes as $sGeneID => &$oGene) {


	echo("Checking $sGeneID ...\n");

	$sGeneDIR = "$sWorkDir/".$oGene['scf']."/".md5($sGeneID);
	$bRevised = false;


	foreach($oGene['mRNAs'] as $smRNAID => &$omRNA) {

		$sRNADIR = $sGeneDIR."/".md5($smRNAID);
		$sFixedGFF = "$sRNADIR/fixed.gff";
		if ( file_exists("$sRNADIR/fixed.gff") ) { //if fixes have been made to the mRNA
			$oFixedGFF = new MappedCDSGFF();
			$oFixedGFF->LoadGFF($sFixedGFF , $sGenome);
			$arrSplitGenes = $oFixedGFF->GetGenes('maker');
			if (count($arrSplitGenes) == 0) continue;
			$arrDeleteGenes[$sGeneID ] = true; //flag for deletion later.
			fwrite($hAppendGFF, file_get_contents($sFixedGFF));
			echo("To be deleted $sGeneID ...\n");
			break;
		}

	}

}

echo(count($arrDeleteGenes) . " genes are chimeric, to be deleted...\n");

$hMakerGFF = fopen($sMakerGFF , "r");
$arrDeletemRNA = array();

while( false !== ($sLn = fgets($hMakerGFF) )) {

	$sLn = trim($sLn);
	if ($sLn == "") continue;
	if ($sLn[0] == '#') {
		fwrite($hOutGFF , $sLn."\n");
	}

	$arrF = explode("\t" , $sLn);

	if ($arrF[2] == "gene") {
		preg_match("/ID=([^;]+)/", $arrF[8] , $arrM);
		if (count($arrM) == 2) {
			$sID = trim($arrM[1]);
			if (array_key_exists($sID , $arrDeleteGenes)) {
				echo("remove gene $sID\n");
				continue;
			}
		}
	} else if ($arrF[2] == "mRNA") {
		preg_match("/Parent=([^;]+)/", $arrF[8] , $arrM);
		if (count($arrM) == 2) {
			$sParentID = trim($arrM[1]);
			if (array_key_exists($sParentID , $arrDeleteGenes)) {
				
				echo("remove mRNA of $sParentID\n");
				preg_match("/ID=([^;]+)/", $arrF[8] , $arrM);
				$arrDeletemRNA[$arrM[1]] = true;
				continue;
			}
		}

	} else {
		preg_match("/Parent=([^;]+)/", $arrF[8] , $arrM);
		if (count($arrM) == 2) {
			$sID = trim($arrM[1]);
			if (array_key_exists($sID , $arrDeletemRNA)) {
				echo("remove part $sID\n");
				continue;
			}
		}
	}


	fwrite($hOutGFF , $sLn."\n");

}

fwrite($hOutGFF, file_get_contents($sAppendGFF) );

?>
