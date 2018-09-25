<?php
$sOrthologTable = "../killi_orthologs_table.ensemblinformed.txt";
$sOrthoGroups = "../UPhO_orthogroups.genesymbol.txt";
$sFish2HumanDir = "maps";

$sMasterOut = "groupid2human.txt";
$sOrthoTableNew = "killi_orthologs_table.ensemblinformed.humanEns.txt";

$arrGroup2FishEnsembl = fnLoadOrthoGroups();
$arrFishEns2Human = fnLoadFish2HumanMap();

//print_r($arrFishEns2Human);

$arrGroup2Human = array();
$hMasterOut = fopen($sMasterOut , 'w');
$hOrthoTableNew = fopen($sOrthoTableNew , 'w');

foreach($arrGroup2FishEnsembl as $sGroupID => $arrFishEnsIDs) {
	$arrHumanIDs = array();
	foreach($arrFishEnsIDs as $sFishEns) {
		if (array_key_exists( $sFishEns ,$arrFishEns2Human )) {
			$sHumanID = $arrFishEns2Human[$sFishEns];
			if (!array_key_exists($sHumanID , $arrHumanIDs) ) {
				$arrHumanIDs[$sHumanID] = 0;
			}

			$arrHumanIDs[$sHumanID]++;
		}
	}

	if (count($arrHumanIDs) > 1 ) {
		echo("Warning: $sGroupID maps to more than 1 human gene: ". implode(",",array_keys($arrHumanIDs) )."\n");
	}

	fwrite($hMasterOut, "$sGroupID\t".implode(",",array_keys($arrHumanIDs) )."\t".implode(",",array_values($arrHumanIDs) )."\n");

	$arrGroup2Human[$sGroupID] = $arrHumanIDs;
}

$hIn = fopen($sOrthologTable , 'r');
$nLn = 0;
while(false !== ($sLn = fgets($hIn))) {
	$sLn = trim($sLn);
	$nLn++;
	$arrF = explode("\t" , $sLn);
	if ($nLn == 1) { //
		$arrF[] = "HumanEnsemblGeneId";
		fwrite($hOrthoTableNew, implode("\t" , $arrF)."\n" );
		continue;
	}

	$arrGroupIds = fnJoinedGroupId2GroupId($arrF[0]);
	$arrPotentialHumanEnsemblIds = array();

	foreach($arrGroupIds as $sGroupId) {
		if (!array_key_exists($sGroupId , $arrGroup2Human) ||  count($arrGroup2Human[$sGroupId])==0 ) {
			continue;
		}

		foreach($arrGroup2Human[$sGroupId] as $sHumanEnsId => $nCount) {
			if (!array_key_exists($sHumanEnsId  , $arrPotentialHumanEnsemblIds ) ) {
				$arrPotentialHumanEnsemblIds[$sHumanEnsId] = 0;
			}

			$arrPotentialHumanEnsemblIds[$sHumanEnsId] += $nCount;
		}
	}

	if (count($arrPotentialHumanEnsemblIds) <= 1 ) {
		$arrF[] = implode(",",array_keys($arrPotentialHumanEnsemblIds) );
		if (count($arrPotentialHumanEnsemblIds) ==0) {
			echo("Notice: $sGroupId does not map to human ortholog\n");
		}
	} else {
		$nMaxCount = max($arrPotentialHumanEnsemblIds);
		foreach($arrPotentialHumanEnsemblIds as $sHumanId => $nCount) {
			if ($nMaxCount == $nCount) {
				$arrF[] = $sHumanId ;
				echo("Warning: multiple human ensembl ids are mapped to $sGroupId, $sHumanId  is chosen\n");
			}
		}
	}

	fwrite($hOrthoTableNew, implode("\t" , $arrF)."\n" );
	
}

function fnJoinedGroupId2GroupId($sJoinedGroupId) {
	preg_match_all("/Group_\d+_\d+/", $sJoinedGroupId, $arrM);

	return $arrM[0];
}

function fnLoadFish2HumanMap() {
	global $sFish2HumanDir;

	$arrFishEns2Human = array();

	$arrFiles = glob("$sFish2HumanDir/*.txt");
	foreach($arrFiles as $sFile) {
		$arrFishEns2Human = $arrFishEns2Human + fnParseFish2HumanMap($sFile);
	}

	return $arrFishEns2Human;
}

function fnParseFish2HumanMap($sFile) {
	$h = fopen($sFile, 'r');
	$sHumanHeadID = "Human gene stable ID";
	$nLn = 0;
	$nHumanCol = 0;
	$nFishCol = 1;
	$arrRet = array();
	while(false !== ($sLn = fgets($h) )) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		$nLn++;
		$arrF = explode("\t" , $sLn);

		if (count($arrF) != 2) {
			continue;
		}
		if ($nLn == 1) { //parse header
			if ($arrF[1] == $sHumanHeadID) {
				$nHumanCol = 1;
				$nFishCol = 0;
			} else if ($arrF[0] == $sHumanHeadID) {
				$nHumanCol = 0;
				$nFishCol = 1;
			} else {
				die("Failed to find header '$sHumanHeadID' in the ensembl map file $sFile\n");
			}

			continue;
		}

		$arrRet[$arrF[$nFishCol]] = $arrF[$nHumanCol];
	}

	return $arrRet;
}

function fnLoadOrthoGroups() {
	global $sOrthoGroups;
	$hIn = fopen($sOrthoGroups , 'r');
	$arrRet= array();
	while(false !== ($sLn = fgets($hIn))) {
		$sLn = trim($sLn);
		$arrF = explode("\t", $sLn);
		if (count($arrF)!= 4) continue;

		$sGroupID = $arrF[0];
		$sGeneSymbol = $arrF[1];
		$sGeneDesc = $arrF[2];
		$arrRet[$sGroupID] = fnParseEnsGenes($arrF[3]);

	}

	return $arrRet;
}

function fnParseGenes($s) {
	$arrF1 = explode(',' , $s);
	$arrRet = array(); //key is taxon name, value is array with gene ids
	foreach($arrF1 as $sGene) {
		if (substr($sGene, 0, 6) == '#trees') continue;
		$arrF2 = explode("|" , $sGene);
		$sTaxon = trim($arrF2[0]);
		$sGeneId = trim($arrF2[1]);

		if (!array_key_exists($sTaxon, $arrRet)) {
			$arrRet[$sTaxon] = array();
		}
		$arrRet[$sTaxon][] = $sGeneId;
	}

	return $arrRet;
}

function fnParseEnsGenes($s) {
	$arrF1 = explode(',' , $s);
	$arrRet = array(); //key is taxon name, value is array with gene ids
	foreach($arrF1 as $sGene) {
		if (substr($sGene, 0, 6) == '#trees') continue;
		$arrF2 = explode("|" , $sGene);
		$sTaxon = trim($arrF2[0]);
		$sGeneId = trim($arrF2[2]);

		if (substr($sGeneId, 0, 3) == 'ENS') { // this is ensembl gene
			$arrF3 = explode('_', $sGeneId);
			$arrRet[] = $arrF3[0];
		}
	}

	return $arrRet;
}

?>
