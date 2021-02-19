<?php
$sTarget = "../../plp.nfzref.upstreamcds.50-60kb.gff";
$sBackground = "../../plp.nfzref.upstreamcds.177827.gff";
$sOrtholog = "/beegfs/group_dv/home/RCui/killifish_genomes/annotation/UPhO_fundulus_alm/assigngenesymbol/map2humanensemblID/killi_orthologs_table.ensemblinformed.humanEns.txt";
$sIDMap = "../NFZgeneIDMap.txt";

//hhvm /beegfs/group_dv/home/RCui/killifish_genomes/pathway_analyses/pathway.php relaxed.p${nPCutoff}.txt 0 sum.txt 0 relaxed.p${nPCutoff}.results > relaxed.p${nPCutoff}.log

//hhvm /beegfs/group_dv/home/RCui/killifish_genomes/pathway_analyses/pathway_genesymbol.php relaxed.p${nPCutoff}.txt 3 sum.txt 3 relaxed.p${nPCutoff}.genesymbol.results > relaxed.p${nPCutoff}.genesymbol.log

$sPathwayScript = "hhvm /beegfs/group_dv/home/RCui/killifish_genomes/pathway_analyses/pathway.php";
$sPathwayGeneSymbolScript = "hhvm /beegfs/group_dv/home/RCui/killifish_genomes/pathway_analyses/pathway_genesymbol.php";

$arrIDMap = fnLoadIDMap($sIDMap);
$arrGroupIDMap = fnGroupIDMap( $sOrtholog, 'NFZgenbank');

fnConvert($sTarget, 'target.txt');
fnConvert($sBackground, 'bg.txt');
exec("ln -sf $sOrtholog ortholog.humanEns.txt");
exec("$sPathwayScript target.txt 0 bg.txt 0 results > results.log");
exec("$sPathwayGeneSymbolScript target.txt 3 bg.txt 3 genesymbol.results > genesymbol.results.log");

function fnConvert($s , $sOut) {
	global $arrIDMap, $arrGroupIDMap;
	$h = fopen($s, 'r');
	$hO = fopen($sOut, 'w');
	while(false !== ($sLn = fgets($h) ) ) {
		$sLn = trim($sLn);
		$arrF = explode("\t", $sLn);
		if (count($arrF) <9) continue;
		$arrAnnot = fnStr2Array($arrF[8]);
		$sNCBIid = $arrAnnot['Parent'];
		if (!array_key_exists($sNCBIid, $arrIDMap) ) {
			echo("$sNCBIid not found\n");
			continue;
		}
		$sRNAid = $arrIDMap[$sNCBIid];
		if (!array_key_exists($sRNAid, $arrGroupIDMap) ) {
			echo("$sNCBIid ($sRNAid) not found\n");
			continue;
		}
		fwrite($hO, $arrGroupIDMap[$sRNAid]."\t$sNCBIid\t$sRNAid\t".$arrAnnot['gene']."\t".$arrAnnot['proteinid']."\n");
	}
	
}

function fnLoadIDMap($sIDMap) {
	$h = fopen($sIDMap, 'r');
	$arrMap = array();
	while(false !== ($sLn = fgets($h) ) ) {
		$sLn = trim($sLn);
		$arrF = explode("\t", $sLn);
		if (count($arrF) < 2 ) continue;
		$arrMap[$arrF[1]] = $arrF[0];
	}
	return $arrMap;
}

function fnGroupIDMap($sOrthologList, $sTaxon) {
	$h = fopen($sOrthologList, 'r');
	$arrMap = array();
	$nLn = -1;
	$nCol = 0;
	while(false !== ($sLn = fgets($h) ) ) {
		$sLn = trim($sLn);
		$arrF = explode("\t", $sLn);
		//print_r($arrF);
		if (count($arrF) < 4 ) continue;
		$nLn++;
		if ($nLn == 0) {
			$arrFlip = array_flip($arrF);
			if (!array_key_exists($sTaxon , $arrFlip) ) {
				//print_r($arrFlip);
				die("$sTaxon not found in $sOrthologList\n");
			}
			$nCol = $arrFlip[$sTaxon];
			continue;
		}

		$sGroupID = $arrF[0];
		list($s1, $s2, $sRNAID ) = explode('|', $arrF[$nCol]);
		$arrMap[$sRNAID] = $sGroupID;
	}

	return $arrMap;
}


function fnStr2Array($s) {
	$arrRet = array();
	$arrF1 = explode(';', $s);
	foreach($arrF1 as $s2) {
		$arrF2 = explode('=', $s2);
		if (count($arrF2) == 2) {
			$arrRet[$arrF2[0]] = $arrF2[1];
		}
	}

	return $arrRet;
}
?>
