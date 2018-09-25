<?php
$nNonNothosTaxa = 14;
$arrSingleCopyIn = array("zebrafish" => true, "AAU" => true, "PLP" => true, "xmac" => true, "poecilia" => true );
$arrMissingIn = array("NOR" => true, "NFZgenbank" => true);
$sUPhoClusters = "UPhO_orthogroups.genesymbol.txt";
$sNothosMissing = "nothos_missing.txt";

$arrIndexByGeneSymbols = array();
$nUnkownCount = 0;

$hIn = fopen($sUPhoClusters ,'r');
$hOut = fopen($sNothosMissing , 'w'); 
while( false !== ($sLn = fgets($hIn) ) ) {
	$sLn = trim($sLn);

	list($sGroupID, $sGeneSymbol, $sDescription, $sGenes) = explode("\t", $sLn); 
	$arrGenes = explode(',' , $sGenes);
	if ($sGeneSymbol == 'unknown') {
		$sGeneSymbol .= ($nUnkownCount++);
	}
 
	if (!array_key_exists($sGeneSymbol , $arrIndexByGeneSymbols)) {
		$arrIndexByGeneSymbols[$sGeneSymbol] = array('ln' => array() , 'genes' => array());
	}

	$arrIndexByGeneSymbols[$sGeneSymbol]['ln'][] = $sLn;
	$arrIndexByGeneSymbols[$sGeneSymbol]['genes'] = array_merge($arrIndexByGeneSymbols[$sGeneSymbol]['genes'] , $arrGenes);

}

echo('read in '. count($arrIndexByGeneSymbols) . ' gene symbols. in which '.$nUnkownCount." are unknown.\n");

foreach($arrIndexByGeneSymbols as $sGeneSymbol => $arrInfo) {

	$arrGenes = $arrInfo['genes'];
	$arrTaxa = array();
	$nGeneCount = 0;
	foreach( $arrGenes as $sGene) {
		if (substr($sGene,0,7) == '#trees_') continue;
		list($sSp, $sNum, $sGeneID) = explode('|' , $sGene);
		if (array_key_exists($sSp , $arrMissingIn)) {
			continue 2; //this family contains nothos, skip.
		}

		$nGeneCount++;
		if (!array_key_exists($sSp, $arrTaxa) ) {
			$arrTaxa[$sSp] = 0;
		}
		$arrTaxa[$sSp] += 1;
	}

	foreach($arrSingleCopyIn as $sSp => $bDum) {
		if (!array_key_exists($sSp , $arrTaxa) ) {
			continue 2;
		} 

		if ($arrTaxa[$sSp]!=1) {
			continue 2;
		}
	}

	foreach($arrInfo['ln'] as $sLn) {
		fwrite($hOut , $sLn . "\n");
	}
}
	


?>
