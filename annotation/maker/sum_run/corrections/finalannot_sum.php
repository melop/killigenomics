<?php
require_once(dirname(__FILE__) . "/lib.php");
require_once(dirname(__FILE__) . "/scores.lib.php");
$sMakerGFF = "maker.genesymbol.supple.gff"; //this gff must have ensemblorthologs annotation. the completeness will be assessed by alignment the gene model to the ensembl ortholog
$sOutGFF = "maker.finalannot.gff";
$sGenome = "query.fa";
$sWorkDir = "./finalannot/";

while(count($argv) > 0) {
    $arg = array_shift($argv);
    switch($arg) {
        case '-R':
            $sGenome = trim(array_shift($argv));
            break;
	case '-M':
	    $sMakerGFF  = trim(array_shift($argv));
	    break;
	case '-O':
	    $sOutGFF  = trim(array_shift($argv));
	    break;
	case '-o':
	    $sWorkDir  = trim(array_shift($argv));
	    break;

    }
}


$oMakerGFF = new MappedCDSGFF();
$oMakerGFF->LoadGFF($sMakerGFF , $sGenome);
$hOutGFF = fopen($sOutGFF , 'w');



$arrGenes = $oMakerGFF->GetGenes('maker'); // get all maker genes
$arrUpdatedAnnot = array();


foreach($arrGenes as $sGeneID => &$oGene) {

	echo("Looking at $sGeneID ...\n");

	$sGeneDIR = "$sWorkDir/".$oGene['scf']."/".md5($sGeneID);

	foreach($oGene['mRNAs'] as $smRNAID => &$omRNA) {

		$smRNAIDMD5= md5($smRNAID);

		//check if annotation already exists
		$oAnnot = $omRNA['annot'];
		if (array_key_exists( 'score',$oAnnot)) {
			continue;
		} else if ( (!array_key_exists( 'ensemblorthologs',$oAnnot)) || $oAnnot['ensemblorthologs'] == '') {
			continue;
		}


		$sRNADIR = $sGeneDIR."/".$smRNAIDMD5;
		$sScoreNotes = "$sRNADIR/annotation.notes";
		$hScoreNotes = fopen($sScoreNotes , 'r');
		$sLn = fgets($hScoreNotes);
		if ($sLn === false) {
			echo("Warning: $sScoreNotes cannot be read.\n");
			continue;
		}

		$arrF = explode("\t", trim($sLn));
		if (count($arrF) != 3) {
			echo("Warning: $sScoreNotes format err.\n");
			continue;
		}

		if ($arrF[1] == 'annot') {
			$arrUpdatedAnnot[$arrF[0]] = fnParseAnnotation($arrF[2]);
		}
		
		
	}
}

echo(count($arrUpdatedAnnot) . " mRNA annotations have been updated.\n" );

$hMakerGFF = fopen($sMakerGFF , 'r');

while( false !== ($sLn = fgets($hMakerGFF) )) {

	$sLn = trim($sLn);
	if ($sLn == "") continue;
	if ($sLn[0] == '#') {
		fwrite($hOutGFF , $sLn."\n");
		continue;
	}

	$arrF = explode("\t" , $sLn);
	$arrOut = $arrF;

	if ($arrF[2] == "mRNA") {
		preg_match("/ID=([^;]+)/", $arrF[8] , $arrM);
		$sID = trim($arrM[1]);
		if (array_key_exists($sID, $arrUpdatedAnnot) ) {
			$arrOut[8] = fnArr2Annot($arrUpdatedAnnot[$sID] + fnParseAnnotation($arrOut[8]));
		}

	}

	fwrite($hOutGFF , implode("\t", $arrOut)."\n");

}


?>
