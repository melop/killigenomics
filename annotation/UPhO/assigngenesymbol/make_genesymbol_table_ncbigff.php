<?php
$sList = "ncbi_gffs.txt";

$hList = fopen($sList, 'r');

while( false !== ($sLn = fgets($hList) ) ) {
	$sLn = trim($sLn);
	if ($sLn == '') continue;
	list($sTaxon, $sGff) = explode("\t" , $sLn);
	$hOut = fopen($sTaxon.".id2genesymbol.map.txt" , 'w');
	$arrID2Symbol = fnProcess($hOut, $sGff); 

}


function fnProcess($hOut, $sGff) {
	$hGff = popen("zcat -f $sGff", 'r');

	while( false !== ($sLn = fgets($hGff ) ) ) {
		if ($sLn == '#') continue;
		$sLn = trim($sLn);
		if ($sLn == '') continue;

		$arrF = explode("\t" , $sLn);

		if (count($arrF) != 9) continue;
		
		if ($arrF[2] == 'mRNA') {
			$arrAnnot = fnParseFields($arrF[8]);
			if (!array_key_exists("ID" , $arrAnnot) ) continue;
			$sGeneSymbol =  (array_key_exists("gene" , $arrAnnot) )? $arrAnnot['gene'] : 'unknown';
			$sGeneFullName =  (array_key_exists("product" , $arrAnnot) )? $arrAnnot['product'] : 'unknown';
			fwrite($hOut , $arrAnnot['ID'] . "\t" . $sGeneSymbol . "\t" . $sGeneFullName  . "\n"  );
			continue;
		}


	}



	
}

function fnParseFields($s) {
	$arrF1 = explode(";" , $s);
	$arrRet = array();
	foreach($arrF1 as $sF) {
		$arrF2 = explode("=" , $sF);
		$arrRet[trim($arrF2[0]) ] = trim($arrF2[1]);
	}

	return $arrRet;
}

?>
