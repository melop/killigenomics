<?php
$sIn = $argv[1];
$nTaxonCount = 16;
$hIn = fopen($sIn ,"r");

$nAllTaxaIncluded = 0;
$nExactOrthologs = 0;
$nSingletons = 0;
$nClusters = 0;

while(false !== ($sLn = fgets($hIn) ) ) {
	$sLn = trim($sLn);
	if ($sLn == '') continue;
	$arrF = explode("\t" , $sLn);

	$arrTaxa = array();
	$nNumSeq = 0;


	foreach($arrF as $sID) {
		$arrF2 = explode('|' , $sID);
		$arrTaxa[$arrF2[0]] = true;
		$nNumSeq++;
	}

	if (count($arrTaxa) ==  $nTaxonCount) {
		$nAllTaxaIncluded++;
	}

	if (count($arrTaxa) ==  $nTaxonCount && $nNumSeq==$nTaxonCount ) {
		$nExactOrthologs++;
	}

	if (count($arrTaxa) == 1 && $nNumSeq==1) {
		$nSingletons++;
	} else {
		$nClusters++;
	}
}

echo("$sIn Clusters (non-singletons): $nClusters AllTaxaRepresented: $nAllTaxaIncluded ExactOrthoGroup: $nExactOrthologs  Singletons: $nSingletons\n");
?>
