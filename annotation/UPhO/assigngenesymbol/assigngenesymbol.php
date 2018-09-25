<?php
$sOrthoGroups = "UPhO_orthogroups.csv";
$sOut = "UPhO_orthogroups.genesymbol.txt";

$arrGeneSymbolMap = glob("*.id2genesymbol.map.txt");
$arrUnknownKeyWord = array("/unknown/" );

$arrGeneSymbolMaps = array();
$arrGeneDescriptionMaps = array();

$hOut = fopen($sOut , "w");

foreach($arrGeneSymbolMap as $sMap) {
	$arrM = explode('.', $sMap);
	$sSp = $arrM[0];
	list( $arrGeneSymbolMaps[$sSp] , $arrGeneDescriptionMaps[$sSp]) = fnParseMap( $sMap ); 
}

//print_r($arrGeneSymbolMaps);

$hOrthoGroups = fopen($sOrthoGroups , 'r');
while( false !== ($sLn = fgets($hOrthoGroups) )) {
	$sLn = trim($sLn);
	if ($sLn == '') continue;
	$arrF = explode(','  , $sLn);
	$sGroupID = substr($arrF[0], 7);
	$arrSymbols = array();
	$arrDescriptions = array();

	foreach($arrF as $sProtID) {
		if (substr($sProtID , 0,7) == '#trees_') continue;
		list($sSp, $sNum, $sID ) = explode('|', $sProtID);
		$arrID = explode('_' , $sID);
		$sID = $arrID[0];

		if (!array_key_exists($sSp , $arrGeneSymbolMaps) ) { //if doesn't exist
			continue;
		}
		if (!array_key_exists($sID , $arrGeneSymbolMaps[$sSp]) ) { //if doesn't exist
			continue;
		}
		$sSymbol = strtolower($arrGeneSymbolMaps[$sSp][$sID]);
		
		$sSymbol = trim(preg_replace("/\(\d+ of .+/", "", $sSymbol));
		$sDes = $arrGeneDescriptionMaps[$sSp][$sID];

		foreach($arrUnknownKeyWord as $sUnknownKeyWord) {
			if (preg_match($sUnknownKeyWord , $sSymbol) ===1 ) {
				continue 2;
			}
			if (!array_key_exists($sSymbol , $arrSymbols)) {
				$arrSymbols[$sSymbol] = 0;
			}
			$arrSymbols[$sSymbol] += 1;

		}

		if (!array_key_exists($sSymbol,  $arrDescriptions)) {
			$arrDescriptions[$sSymbol] = $sDes;
		}

		$arrDescriptions[$sSymbol] = (strlen($arrDescriptions[$sSymbol]) < strlen($sDes) && strlen($arrDescriptions[$sSymbol]) > 10 )? $arrDescriptions[$sSymbol] : $sDes;
	}



	arsort( $arrSymbols , SORT_NUMERIC);
	$arrSymbolsFlip = array_keys($arrSymbols);

	if (count($arrSymbolsFlip) > 1) {
		echo("warning: more than one possible gene symbols found for group $sGroupID : " . implode(', ' , $arrSymbolsFlip ) . "\n");
		
	}
	
	$sFinalSymbol = 'unknown';
	if (count($arrSymbolsFlip) > 0) {
		$sFinalSymbol = $arrSymbolsFlip[0];
	}

	fwrite($hOut , "$sGroupID\t".$sFinalSymbol."\t".$arrDescriptions[$sFinalSymbol]."\t".$sLn."\n");
}


function fnParseMap( $sMap ) {
	$hMap = fopen($sMap, 'r');
	$arrRet = array(array(), array());
	while( false !== ($sLn = fgets($hMap) )) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;

		$arrF = explode("\t" , $sLn);
		$arrRet[0][$arrF[0]] = $arrF[1];
		$arrRet[1][$arrF[0]] = $arrF[2];
	}

	return $arrRet;
}

?>
