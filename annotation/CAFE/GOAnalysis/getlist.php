<?php
$sFamF = "../reports/errest_sp_allgenes/cafe_final_summary_fams.txt"; //PATH TO THE CAFE OUTPUT
$sMCLCluster = "../out.seq.mci.I20"; //PATH TO THE INPUT MCI FILE
$nFamSizeFilter = 10000000; //Family size limit.

$hFam = fopen($sFamF , 'r');
list($arrClusterList , $arrExcludeFromBackground) = fnLoadClusters($sMCLCluster);



while( false !== ($sLn = fgets($hFam) ) ) {
	$sLn = trim($sLn);
	if ($sLn[0] == '#') continue;
	if (substr($sLn  ,0, 15) == 'Overall rapid :' ) continue;
	list($sNodeName, $sGeneFamList) = explode("\t" , $sLn);
	$sNodeName = fnNormalizeNodeName($sNodeName);
	//echo($sNodeName."\n");
	$arrGeneFams = explode(',', $sGeneFamList);
	$sExpansion = "$sNodeName.expandedgenes.famsizefilter$nFamSizeFilter.txt";
	$sContraction = "$sNodeName.contractedgenes.famsizefilter$nFamSizeFilter.txt";

	$arrFileHandle = array();
	$arrFileHandle[0] = fopen($sExpansion , "w");
	$arrFileHandle[1] = fopen($sContraction , "w");

	foreach($arrGeneFams as $sGeneFam) {
		
		$oFam =  fnParseGeneFam($sGeneFam);
		$nFamID = $oFam[1];
		$nFileToUse = ($oFam[2]=='+')? 0 : 1;
		if (!array_key_exists( $nFamID , $arrClusterList )) {
			die("FamID not found $nFamID.\n");
		}
		fwrite($arrFileHandle[$nFileToUse] , $arrClusterList[$nFamID] );
		
	}

}

exec('cat '. $sMCLCluster .' | grep -oP  "zebrafish[^\s]+" | grep -oP "ENS[^_]+" > background.zebrafish.genes.txt');
$hBackground = fopen('background.zebrafish.genes.txt' , 'r');
$hBackgroundFilter = fopen("background.zebrafish.genes.famsizefilter$nFamSizeFilter.txt" , 'w');

while( false !== ($sLn = fgets($hBackground) ) ) {
	$sLn = trim($sLn);
	if (!array_key_exists($sLn , $arrExcludeFromBackground) ) {
		fwrite($hBackgroundFilter , $sLn."\n");
	}
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

function fnLoadClusters($sMCLCluster) {
	global $nFamSizeFilter;
	$hF = fopen($sMCLCluster , "r");
	$nFamID = 0;
	$arrRet = array();
	$arrExclude = array();
	while( false !== ($sLn = fgets($hF) ) ) {
		$sLn = trim($sLn);
		if ($sLn =='') continue;
		$nFamID++;
		$arrRet[$nFamID] = "";
		$arrFs = explode("\t" , $sLn);

		$bExcludeFam = (count($arrFs) > $nFamSizeFilter);

		if ($bExcludeFam) {
			echo("family $nFamID size is > $nFamSizeFilter, exclude.\n ");
		}

		foreach($arrFs as $sGeneID) {
			$arrGeneIDF = explode('|' , $sGeneID);
			if (count($arrGeneIDF) != 3 ) continue;
			if ($arrGeneIDF[0] != 'zebrafish') continue; //only consider zebrafish Ensembl IDs

			$arrEnsID = explode('_' , $arrGeneIDF[2]);
			$sEnsID = $arrEnsID[0];

			if (!$bExcludeFam) {

				$arrRet[$nFamID] .= $sEnsID."\n";
			} else {
				$arrExclude[$sEnsID] = true;
			}
		}
	}

	return array($arrRet , $arrExclude);
}
?>
