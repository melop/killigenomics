<?php
$sList = "ensembl_gffs.txt";

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
	$arrGenes = array(); //key is gene id , values are gene symbol and gene full name
	$arrRNA2GeneMap  = array();
	$arrProtein2RNAMap = array();
	while( false !== ($sLn = fgets($hGff ) ) ) {
		if ($sLn == '#') continue;
		$sLn = trim($sLn);
		if ($sLn == '') continue;

		$arrF = explode("\t" , $sLn);

		if (count($arrF) != 9) continue;
		
		if ($arrF[2] == 'mRNA') {
			$arrAnnot = fnParseFields($arrF[8]);
			if (!array_key_exists("ID" , $arrAnnot) ) continue;
			if (!array_key_exists("Parent" , $arrAnnot) ) continue;
			$arrRNA2GeneMap[$arrAnnot['ID']] = $arrAnnot['Parent'];
			continue;
		}

		if ($arrF[2] == 'CDS') {
			$arrAnnot = fnParseFields($arrF[8]);
			if (!array_key_exists("protein_id" , $arrAnnot) ) continue;
			if (!array_key_exists("Parent" , $arrAnnot) ) continue;
			$arrProtein2RNAMap[$arrAnnot['protein_id']] = $arrAnnot['Parent'];
			continue;
		}


		if ($arrF[2] != 'gene' ) {
			continue;

		} 
		$arrAnnot = fnParseFields($arrF[8]);
		if (!array_key_exists("ID" , $arrAnnot) ) continue;
		$arrGenes[$arrAnnot['ID']] = array('unknown','unknown');
		if (array_key_exists("Name" , $arrAnnot)) {
			$arrGenes[$arrAnnot['ID']][0] = $arrAnnot['Name'];
		} 

		if (array_key_exists("description" , $arrAnnot)) {
			$arrGenes[$arrAnnot['ID']][1] = $arrAnnot['description'];
		} 

	}

	foreach($arrProtein2RNAMap as $sProteinID => $sTranscriptID) {
		//find the gene
		if (!array_key_exists($sTranscriptID , $arrRNA2GeneMap) ) {
			echo("warning: Transcript $sTranscriptID not found in rna2gene map.\n");
			continue;
		}

		$sGeneID = $arrRNA2GeneMap[$sTranscriptID];
		if (!array_key_exists($sGeneID , $arrGenes) ) {
			echo("warning: Gene ID $sGeneID undefined.\n");
			continue;
		}

		fwrite($hOut , $sProteinID . "\t" . $arrGenes[$sGeneID][0] . "\t" . $arrGenes[$sGeneID][1] . "\n"  );
		
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
