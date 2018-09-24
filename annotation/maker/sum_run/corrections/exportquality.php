<?php
require_once(dirname(__FILE__) . "/lib.php");
$sMakerGFF = "maker.refined.removeisoform.gff";
$sGenome = "query.fa";
$sOutTab = "quality.refined.removeisoform.txt";


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
	    $sOutTab  = trim(array_shift($argv));
	    break;


    }
}


$hOutTab = fopen($sOutTab , "w");

$oMakerGFF = new MappedCDSGFF();


$oMakerGFF->LoadGFF($sMakerGFF , $sGenome);


$arrGenes = $oMakerGFF->GetGenes('maker'); // get all maker genes

$arrFields = array( 	
	'score', 
	'5primemiss',
	'5primemissperc',
	'3primemiss',
	'3primemissperc',
	'5primeuncertain',
	'3primeuncertain',
	'completeness',
	'hicov',
	'smallinsertionlen',
	'smallinsertioncount',
	'biginsertionlen',
	'biginsertioncount',
	'hugeinsertionlen',
	'hugeinsertioncount',

	'smalldeletionlen',
	'smalldeletioncount',
	'bigdeletionlen',
	'bigdeletioncount',
	'hugedeletionlen',
	'hugedeletioncount'
);

fwrite($hOutTab ,"GeneID\tRNAID");
foreach($arrFields as &$sField) {
	fwrite($hOutTab , "\t$sField");
} 
	fwrite($hOutTab , "\n");

foreach($arrGenes as $sGeneID => &$oGene) {

	foreach($oGene['mRNAs'] as $smRNAID => &$omRNA) {
		$oA = $omRNA['annot'];
		$arrOut = array( $sGeneID ,
				 $smRNAID );

		foreach($arrFields as &$sField) {
			$arrOut[] = fnShow($oA, $sField);
		} 
			
				

		fwrite($hOutTab , implode("\t" , $arrOut)."\n");
	}

}

function fnShow(&$oA, $sName) {
	return array_key_exists($sName, $oA)? $oA[$sName] : "NA";
}
?>
