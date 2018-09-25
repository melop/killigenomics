<?php
$sFamF = "../reports/errest_sp_allgenes_noNFZ/cafe_final_summary_fams.txt";
$nFamSizeFilter = 1000000; //only analyze smaller families
$sCountOut = "node_counts.txt";
$hCountOut = fopen($sCountOut , "w");

$hFam = fopen($sFamF , 'r');

$arrCounts = array();
$arrGeneCounts = array();

while( false !== ($sLn = fgets($hFam) ) ) {
	$sLn = trim($sLn);
	if ($sLn[0] == '#') continue;
	if (substr($sLn  ,0, 15) == 'Overall rapid :' ) continue;
	list($sNodeName, $sGeneFamList) = explode("\t" , $sLn);
	$sNodeName = fnNormalizeNodeName($sNodeName);
	$arrCounts[$sNodeName]  = array( 0, 0 , 0 ,0);
	//echo($sNodeName."\n");
	$arrGeneFams = explode(',', $sGeneFamList);


	foreach($arrGeneFams as $sGeneFam) {
		
		$oFam =  fnParseGeneFam($sGeneFam);
		$nFamID = $oFam[1];
		$nNumChange = $oFam[3];
		$nIdx = ($oFam[2]=='+')? 0:1;
		$nIdx2 = ($oFam[2]=='+')? 2:3;
		$arrCounts[$sNodeName][$nIdx] += 1;
		$arrCounts[$sNodeName][$nIdx2] += $nNumChange;

	}

	fwrite($hCountOut , $sNodeName . "\t". implode("\t", $arrCounts[$sNodeName]) ."\n");
}



function fnNormalizeNodeName($s) {
	return str_replace(array('<', '>', ':') , array('Node_' , '','') , $s );
}

function fnParseGeneFam($s) {
	preg_match("/(\d+)\[([+-])(\d+)\S*\]/", $s, $arrR);
	if (count($arrR) != 4 ) {
		die("Cannot parse family description: $s\n");
	}

	return  $arrR;
}

?>
