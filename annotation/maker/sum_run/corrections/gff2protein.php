<?php
require_once(dirname(__FILE__) . "/lib.php");
$sMakerGFF = "chimericfixed.gff";
$sGenome = "query.fa";
$sOutFasta = "chimericfixed.protein.fa";
$sOutTransFasta = "chimericfixed.transcript.fa";

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
	    $sOutFasta  = trim(array_shift($argv));
	    break;
	case '-t':
	    $sOutTransFasta  = trim(array_shift($argv));
	    break;

    }
}


$hOutFasta = fopen($sOutFasta , "w");
$hOutTransFasta = fopen($sOutTransFasta , "w");

$oMakerGFF = new MappedCDSGFF();


$oMakerGFF->LoadGFF($sMakerGFF , $sGenome);


$arrGenes = $oMakerGFF->GetGenes('maker'); // get all maker genes

$arrDeleteGenes = array();
foreach($arrGenes as $sGeneID => &$oGene) {

	foreach($oGene['mRNAs'] as $smRNAID => &$omRNA) {
		$oExtracted = $oMakerGFF->ExtractmRNASequence('maker' , $smRNAID);
		fwrite($hOutFasta , ">$smRNAID\n".wordwrap($oExtracted['AA'] , 75, "\n", true)."\n");
		fwrite($hOutTransFasta , ">$smRNAID\n".wordwrap($oExtracted['mRNA'] , 75, "\n", true)."\n");
	}

}

?>
