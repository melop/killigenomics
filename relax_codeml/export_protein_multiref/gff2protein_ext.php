<?php
require_once(dirname(__FILE__) . "/lib.php");
$sMakerGFF = "chimericfixed.gff";
$sGenome = "query.fa";
$sOutFasta = "chimericfixed.protein.fa";
$sOutTransFasta = "chimericfixed.transcript.fa";

//This script also delete 5' and 3'  uncertainties and will use N to mask all huge gaps. small and big gaps are left in 

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
		$sTrimAA = $oExtracted['AA'];
		$sTrimDNA = $oExtracted['mRNA'];

		//mask huge gaps:
		if (array_key_exists('hugeinsertionlist', $omRNA['annot'])) {
			$arrHugeInsList = explode(",", $omRNA['annot']['hugeinsertionlist']);
			foreach($arrHugeInsList as $sHugeGap) {
				$sHugeGap = trim($sHugeGap);
				if ($sHugeGap == "") continue;
				list($nGapStart, $nGapEnd) = explode('-', $sHugeGap);
				$nGapLen = abs($nGapEnd-$nGapStart) + 1;
				$nDNAGapLen = $nGapLen * 3;
				$sTrimAA = substr_replace($sTrimAA , str_repeat('X', $nGapLen) , ($nGapStart-1) ,  $nGapLen);
				$sTrimDNA = substr_replace($sTrimDNA , str_repeat('N', $nDNAGapLen) , ($nGapStart-1)*3 ,  $nDNAGapLen);
			}
		}

		$n5PrimeUncertain = 0;
		$n3PrimeUncertain = 0;
		if (array_key_exists('5primeuncertain', $omRNA['annot'])) {
			$n5PrimeUncertain = intval($omRNA['annot']['5primeuncertain']);
		}

		if (array_key_exists('3primeuncertain', $omRNA['annot'])) {
			$n3PrimeUncertain = intval($omRNA['annot']['3primeuncertain']);
		}

		//trim off 3'
		$sTrimAA = substr($sTrimAA , 0 , strlen($sTrimAA) - $n3PrimeUncertain );
		$sTrimDNA = substr($sTrimDNA , 0 , strlen($sTrimDNA) - ($n3PrimeUncertain*3) );
		//trim off 5'
		$sTrimAA = substr($sTrimAA , $n5PrimeUncertain );
		$sTrimDNA = substr($sTrimDNA ,  $n5PrimeUncertain * 3 );

		fwrite($hOutFasta , ">$smRNAID\n".wordwrap($sTrimAA  , 75, "\n", true)."\n");
		fwrite($hOutTransFasta , ">$smRNAID\n".wordwrap($sTrimDNA , 75, "\n", true)."\n");
	}

}

?>
