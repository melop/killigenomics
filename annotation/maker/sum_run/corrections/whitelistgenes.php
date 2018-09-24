<?php
require_once(dirname(__FILE__) . "/lib.php");

//produce a list of mRNA ids that are likely complete genes
// for partial genes, look to see if there is another gene assigned the same gene symbol, if so, add their "completeness" and round the the nearest integer
// for example, "genesymbol (1 of 2) completeness=0.80" and "genesymbol (2 of 2) completeness=0.40", genesymbol sum_completeness = 1.2, estimate of gene copy = round(1.2) = 1.
// keep only genesymbol (1 of 2) for CAFE analysis.
// if a gene symbol only as 1 gene model, keep it regardless of completeness
//$sMakerGFF = "maker.finalannot.gff";
//$sOut = "cafe_whitelist.txt";

$sMakerGFF = "maker.finalannot.improved.gff";
$sOut = "cafe_whitelist.improved.txt";

$sGenome = "query.fa";
$nCompletnessThreshUnknownGenes = 0.5;

$oMakerGFF = new MappedCDSGFF();
$oMakerGFF->LoadGFF($sMakerGFF , $sGenome);

$arrGenes = $oMakerGFF->GetGenes('maker'); // get all maker genes

$hOut = fopen($sOut, 'w');

$arrIdxGeneSymbol = array();
$nUnknownCount = 0;
foreach($arrGenes as $sGeneID => &$oGene) {


	foreach($oGene['mRNAs'] as $smRNAID => &$omRNA) {

		//check if annotation already exists
		$oAnnot = $omRNA['annot'];
		if (!array_key_exists("completeness" , $oAnnot)) { //completeness cannot be computed, usually are low quality models, discard
			continue; 
		}

		if (!array_key_exists("gene" , $oAnnot)) {
			continue;
		}

		$sCleanGeneSymbol = fnParseGeneSymbol($oAnnot['gene']);
		if ($sCleanGeneSymbol == 'unkown' || $sCleanGeneSymbol == 'unknown') { //i misspelled "unknown" as "unkown" oops...
			if ($oAnnot['completeness'] < $nCompletnessThreshUnknownGenes) {
				continue;
			} else {
				$sCleanGeneSymbol = "unknown".($nUnknownCount++);
			}
		}

		if (!array_key_exists($sCleanGeneSymbol  , $arrIdxGeneSymbol) ) {
			$arrIdxGeneSymbol[$sCleanGeneSymbol] = array(); 
		}

		$arrIdxGeneSymbol[$sCleanGeneSymbol][] = array('spgeneid' => $oAnnot['spgeneid'] , 'completeness' => $oAnnot['completeness'] );
	}

}

echo("Gene symbol count: " . count($arrIdxGeneSymbol) . PHP_EOL );

foreach($arrIdxGeneSymbol as $sGeneSymbol => &$arrGenes) {
	if ( count($arrGenes) == 1 ) {
		
		fwrite($hOut , $arrGenes[0]['spgeneid']."\t$sGeneSymbol\t".$arrGenes[0]['completeness']. "\n");
		continue;
	} 

	$nSumComplete = 0;
	$arrID2Comp = array();
	foreach($arrGenes as $oGene) {
		$nSumComplete += $oGene['completeness'];
		$arrID2Comp[$oGene['spgeneid']] = $oGene['completeness'];
	}

	arsort($arrID2Comp);
	$nKeepGenes = round($nSumComplete );
	$arrIDs = array_keys($arrID2Comp);
	$arrComps = array_values($arrID2Comp);
	for($nKeep=0;$nKeep<$nKeepGenes; $nKeep++) {
		fwrite($hOut , $arrIDs[$nKeep]."\t$sGeneSymbol\t".$arrComps[$nKeep]. "\n");
	}
}

function fnParseGeneSymbol($s) {
	$arrM = explode(' ', $s);

	return $arrM[0];
}
?>
