<?php
$sMCI = "/beegfs/group_dv/home/RCui/killifish_genomes/annotation/UPhO_final/mclclusters/out.seq.mci.I20";
$sOut = "out.seq.mci.I20";

$arrWhiteLists = array('AAU' => "/beegfs/group_dv/home/RCui/killifish_genomes/annotation/maker/AAU_final_1/run_blastmerge_cov_exonerateOn/corrections/cafe_whitelist.txt",
'CTO' => "/beegfs/group_dv/home/RCui/killifish_genomes/annotation/maker/CTO_final_1/run_pilon_exonerateOn/corrections/cafe_whitelist.txt",
'NOR' => "/beegfs/group_dv/home/RCui/killifish_genomes/annotation/maker/NOR_final_1/run_pilon_exonerateOn/corrections/cafe_whitelist.txt",
'PLP' => "/beegfs/group_dv/home/RCui/killifish_genomes/annotation/maker/PLP_final_1/run_pilon_exonerateOn/corrections/cafe_whitelist.txt");

$arrGeneFilters = fnLoadGeneFilters($arrWhiteLists); //only include genes still included in the deduplicated GFF files.

$hMCI = fopen($sMCI , 'r');
$hOut = fopen($sOut, 'w');

$nRemovedCount = 0;

while(false !== ($sLn = fgets($hMCI))) {
	$sLn = trim($sLn);
	$arrF = explode("\t", $sLn);

	$arrPass = array();
	foreach($arrF as $sGeneID) {
		$arrGeneFields = explode('|', $sGeneID);
		if (fnIsGeneIncluded($arrGeneFields[0] , $arrGeneFields[1])) {
			$arrPass[] = $sGeneID;
		} else {
			$nRemovedCount++;
		}
	}

	if (count($arrPass) > 0 ) {
		fwrite($hOut , implode("\t", $arrPass) . "\n" );
	} else {
		echo("Completely removed: $sLn\n");
	}
}

echo("$nRemovedCount genes removed.\n");

function fnLoadGeneFilters($arrWhiteLists) {
	$arrRet = array();

	foreach($arrWhiteLists as $sSp => $sGFF) {
		$arrRet[$sSp] = array();
		$hGFF = fopen($sGFF, 'r');
		while( false !== ($sLn = fgets($hGFF) )) {
			$sLn = trim($sLn);
			$arrF = explode("\t" , $sLn);
			$nGeneID = $arrF[0];
			$arrRet[$sSp][$nGeneID] = true;
		}
	}

	return $arrRet;
}

function fnIsGeneIncluded($sSp, $nGeneID) {
	global $arrGeneFilters;
	if (!array_key_exists($sSp , $arrGeneFilters)) return true;
	return array_key_exists($nGeneID , $arrGeneFilters[$sSp]);
}

?>
