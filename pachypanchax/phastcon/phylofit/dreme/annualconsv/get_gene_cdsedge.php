<?php
//outputs the edges of the first and the last CDS of any gene.

$sGFF = "cds.sorted.gff";
$sOut = "gene.cds_edges.gff";
$h = fopen($sGFF, 'r');
$hO = fopen($sOut , 'w');

$arrRet = array();
while( false !== ($sLn = fgets($h) ) ) {
	$sLn = trim($sLn);
	$arrF = explode("\t", $sLn );
	if (count($arrF) < 9 ) continue;

	preg_match('/Parent=([^;]+)/', $arrF[8], $arrM);
	if (count($arrM) != 2 ) continue;
	$sRNAID = $arrM[1];
	if (!array_key_exists($sRNAID , $arrRet) ) {
		$arrRet[$sRNAID] = array('annot' => $arrF[8], 'coords' => array(), 'chr'=>$arrF[0], 'strand' => $arrF[6], 'annotator' => $arrF[1]);
	}

	$arrRet[$sRNAID]['coords'][] = $arrF[3];
	$arrRet[$sRNAID]['coords'][] = $arrF[4];
}

foreach($arrRet as $sRNAID => $arrInfo) {
	$arrOut = array();
	$arrOut[] = $arrInfo['chr'];
	$arrOut[] = $arrInfo['annotator'];
	$arrOut[] = 'CDS_range';
	$arrOut[] = min($arrInfo['coords']);
	$arrOut[] = max($arrInfo['coords']);
	$arrOut[] = '.';
	$arrOut[] = $arrInfo['strand'];
 	$arrOut[] = '.';
	$arrOut[] = $arrInfo['annot'];
	fwrite($hO, implode("\t", $arrOut)."\n");
}

?>
