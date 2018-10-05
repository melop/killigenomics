<?php
//translates the featurecount mRNA IDs into orthogroup ids
//if a gene is not included in the orthogroup, do not output it.
$sRefDef = "refs.txt";
$sOrthologList = "orthologs.improved.txt";

$sDir = "../";
$sOutPrefix = "counts_";

/*============ Load the ref genomes ======================================*/
$hRefDef = fopen($sRefDef , "r");
$arrRefGTFs = array();

$nLn = -1;
echo("Parsing reference genome definitions...\n");
while( ($sLn=fgets($hRefDef))!==false ) {
        $sLn = trim($sLn);
        if ($sLn == "") {
                continue;
        }
        $nLn++;
        if ($nLn==0) {
                continue; // skip header
        }
        
        $arrFields = explode("\t", $sLn);
        $arrRefGTFs[$arrFields[0]] = $arrFields[2];
}

$arrRNAID2SpGeneIDMap = array();

foreach($arrRefGTFs as $sSp => $sGFF) {
	$hGFF = popen('grep -P "\tmRNA\t" '.$sGFF , 'r');
	$arrRNAID2SpGeneIDMap[$sSp] = array();
	while( ($sLn=fgets($hGFF))!==false ) {
		$arrF = explode("\t" , trim($sLn));
		$sAnn = $arrF[8];
		preg_match("/ID=([^;]+)/", $sAnn, $arrM);
		if (count($arrM ) != 2) continue;
		$sRNAID = $arrM[1];
		preg_match("/spgeneid=([^;]+)/", $sAnn, $arrM);
		if (count($arrM ) != 2) continue;
		$sSpGeneID = $arrM[1];
		$arrRNAID2SpGeneIDMap[$sSp][$sRNAID] = $sSpGeneID;
	}

	echo("Loaded ".count($arrRNAID2SpGeneIDMap[$sSp])." genes for species $sSp\n");
}

/*============ Load the ortholog definitions  ======================================*/
$hOrthologList = fopen($sOrthologList , "r");
$arrOrthologs = array();
$arrSpp = array();


$nLn = -1;
echo("Parsing ortholog definitions...\n");
while( ($sLn=fgets($hOrthologList))!==false ) {
        $sLn = trim($sLn);
        if ($sLn == "") {
                continue;
        }
        $nLn++;
        if ($nLn==0) {
                $arrSpp = array_slice( explode("\t" , $sLn) , 3);
		$arrOrthologs = array_combine($arrSpp , array_fill(0, count($arrSpp), array() ) ); 
                if ($arrSpp != array_keys($arrRefGTFs) ) {
                        die("The headers in the ortholog definition file do not match the reference genome names defined in the reference genome file.\n");
                }
                continue; // skip header
        }
        
        $arrFields = explode("\t", $sLn);
        $arrOrIDs =  array_slice($arrFields , 3);
	foreach($arrSpp as $nIdx => $sSp) {
		$arrOrthologs[$sSp][$arrOrIDs[$nIdx]] = $arrFields[0]."\t".$arrFields[1]."\t".$arrFields[2];
	}
}


echo("Loaded ". count($arrOrthologs[$arrSpp[0]]) ." ortholog definitions\n" );
//print_r(array_keys($arrOrthologs));
//die();
// load the counts

foreach($arrSpp as $sSp) {
	$sIn = "$sDir/$sSp/rnaseq_map/mappedPE/featurecounts/counts.txt";
	$sOut = "$sOutPrefix.$sSp.txt";
	if (!file_exists($sIn)) {
		die("$sIn not found!\n");
	}

	$hIn = fopen($sIn , 'r');
	$hOut = fopen($sOut , 'w');
	while( false !== ($sLn = fgets($hIn)) ) {
		$sLn = trim($sLn);
		if ($sLn == '' || $sLn[0] == '#') continue;
		$arrF = explode("\t", $sLn);
		if ($arrF[0] == "Geneid") {
			$arrF[0] = "OrthoID\tGeneSymbol\tGeneName\tSpGeneId\tRNAID";
			fwrite($hOut , implode("\t" , $arrF) . "\n" );
			continue;
		}

		if (!array_key_exists($arrF[0], $arrRNAID2SpGeneIDMap[$sSp])) continue;
		$sSpGeneId = $arrRNAID2SpGeneIDMap[$sSp][$arrF[0]];
		if (!array_key_exists( $sSpGeneId, $arrOrthologs[$sSp])) continue;
		$sGroupId = $arrOrthologs[$sSp][$sSpGeneId];
		$arrF[0] = "$sGroupId\t$sSpGeneId\t$arrF[0]";
		fwrite($hOut , implode("\t" , $arrF) . "\n" );
		
	}
}


?>
