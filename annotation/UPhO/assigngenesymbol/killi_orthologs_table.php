<?php
/*
Assign gene symbols to killifish ortholog IDs "Group_0_0 ...."
*/

$sOrthoGroups = "UPhO_orthogroups.genesymbol.txt"; #UPhO ortholog output with gene symbols
$sOut = "killi_orthologs_table.ensemblinformed.txt"; #output

$arrEnsemblOrthologIDMaps = array('ensembl_fish_ortholog_map1.tsv' , 'ensembl_fish_ortholog_map2.tsv');
//ensembl ortholog table. used to inform orthology in case an orthologous group is broken into several parts by UPhO. Download this from Ensembl Biomart


$arrEnsOrthMap = fnReadEnsemblOrthologMap();



/*$arrWhiteLists = array('AAU' => "/beegfs/group_dv/home/RCui/killifish_genomes/annotation/maker/AAU_final_1/run_blastmerge_cov_exonerateOn/corrections/cafe_whitelist.txt",
'CTO' => "/beegfs/group_dv/home/RCui/killifish_genomes/annotation/maker/CTO_final_1/run_pilon_exonerateOn/corrections/cafe_whitelist.txt",
'NOR' => "/beegfs/group_dv/home/RCui/killifish_genomes/annotation/maker/NOR_final_1/run_pilon_exonerateOn/corrections/cafe_whitelist.txt",
'PLP' => "/beegfs/group_dv/home/RCui/killifish_genomes/annotation/maker/PLP_final_1/run_pilon_exonerateOn/corrections/cafe_whitelist.txt");
*/

$arrWhiteLists = array('AAU' => "/beegfs/group_dv/home/RCui/killifish_genomes/annotation/maker/AAU_final_1/run_blastmerge_cov_exonerateOn/corrections/cafe_whitelist.improved.txt",
'CTO' => "/beegfs/group_dv/home/RCui/killifish_genomes/annotation/maker/CTO_final_1/run_pilon_exonerateOn/corrections/cafe_whitelist.improved.txt",
'NOR' => "/beegfs/group_dv/home/RCui/killifish_genomes/annotation/maker/NOR_final_1/run_pilon_exonerateOn/corrections/cafe_whitelist.improved.txt",
'PLP' => "/beegfs/group_dv/home/RCui/killifish_genomes/annotation/maker/PLP_final_1/run_pilon_exonerateOn/corrections/cafe_whitelist.improved.txt");


$arrSpecies = array("AAU" , "CTO" , "NOR", "PLP");

$arrGeneFilters = fnLoadGeneFilters($arrWhiteLists); //only include genes still included in the deduplicated GFF files.

$hIn = fopen($sOrthoGroups , 'r');
$hOut = fopen($sOut , 'w');

$nGroupCount = 0;

$arrIncludedOrtho = array();

$arrPotentialCombinable = array(); //saves the incomplete groups, with fewer than 4 ingroup species; these can potentially be combined by looking into the ensembl ortholog definition. The key is super group name. e.g. Group_0_0 => Group_0 is super group name.

fwrite($hOut, "Group_Id\tGene_Symbol\tGene_Name\t".implode("\t", $arrSpecies  )."\n");
while(false !== ($sLn = fgets($hIn))) {
	$sLn = trim($sLn);
	$arrF = explode("\t", $sLn);
	if (count($arrF)!= 4) continue;

	$sGroupID = $arrF[0];
	$sGeneSymbol = $arrF[1];
	$sGeneDesc = $arrF[2];
	$arrGenes = fnParseGenes($arrF[3]);

	/*if ($sGeneSymbol == 'sirt6') {
		print_r($arrGenes);
		die();
	}*/

	//now check that it's single copy in every target taxon:

	$arrGeneIDs = array();
	$nFoundSingleCopySpecies = count($arrSpecies);
	$bAnyMultipleCopy = false;
	foreach($arrSpecies as $sSp) {
		if (!array_key_exists($sSp , $arrGenes) || (count($arrGenes[$sSp])!=1) ) {
			$nFoundSingleCopySpecies--; //problem occur for this species
			if (array_key_exists($sSp , $arrGenes) && count($arrGenes[$sSp]) > 1) {
				$bAnyMultipleCopy = true;
			}
		} else {
			$arrGeneIDs[] = $arrGenes[$sSp][0];
		}
	}

	if ( $nFoundSingleCopySpecies > 0 && $nFoundSingleCopySpecies < count($arrSpecies) && (!$bAnyMultipleCopy)) { 
                 // if the number of ingroup species found to be single copied is lower than expected, keep it to look for later.
		$arrEnsGenes = fnParseEnsGenes($arrF[3]); // also get the accompanying ensembl gene ids
		if (count($arrEnsGenes) > 0 ) {
			$sSuperGroupName = implode('_', array_slice( explode('_', $sGroupID ), 0, 2 ));
			if (!array_key_exists($sSuperGroupName, $arrPotentialCombinable) ) {
				$arrPotentialCombinable[$sSuperGroupName] = array();
			}

			$arrPotentialCombinable[$sSuperGroupName][$sGroupID] = array( 'ingroup_genes' => $arrGenes , 'ensemblgenes' => $arrEnsGenes, 'genesymbol' => $sGeneSymbol, 'genedesc' => $sGeneDesc);
		}

		continue;
	}

	if ( $nFoundSingleCopySpecies !=  count($arrSpecies)) {
		continue;
	}

	if ($bAnyMultipleCopy) {
		continue;
	}

	$sUID = implode("_", $arrGeneIDs  );
	if (array_key_exists($sUID, $arrIncludedOrtho) ) continue;
	$arrIncludedOrtho[$sUID] = true;
	fwrite($hOut, "$sGroupID\t$sGeneSymbol\t$sGeneDesc\t".implode("\t", $arrGeneIDs  )."\n");
	$nGroupCount++;
}

echo("$nGroupCount orthologs included.\n");

echo("Now try to combine split orthogroups by using ensembl ortholog definitions\n");

//print_r($arrPotentialCombinable);
//die();

//now, go through each super family, and do pairwise search to see if we could combine any
foreach($arrPotentialCombinable as $sSuperGroup => $arrSubGroups) {
	

	foreach($arrSubGroups as $sGroup1 => $arrGroup1) {
		$arrIngroupTaxa1 = array_intersect( array_keys($arrGroup1['ingroup_genes']) , $arrSpecies);
		$arrEnsGenes1 = $arrGroup1['ensemblgenes'];
		$arrInGenes1 = sub_array( $arrGroup1['ingroup_genes'] , $arrIngroupTaxa1);

		foreach($arrSubGroups as $sGroup2 => $arrGroup2) {
			if (strcmp($sGroup1 , $sGroup2) >= 0) continue;
			$arrIngroupTaxa2 = array_intersect( array_keys($arrGroup2['ingroup_genes']) , $arrSpecies);
			$arrInTaxa =array_merge( $arrIngroupTaxa1 , $arrIngroupTaxa2);
			if (count(array_unique($arrInTaxa)) != count($arrSpecies) ) {
				continue;
			}
			// it's possible, check if the two sets of ensemble ids have orthologs
			$arrEnsGenes2 = $arrGroup2['ensemblgenes'];

			$arrInGenes2 = sub_array( $arrGroup2['ingroup_genes'] , $arrIngroupTaxa2);

			$arrInGenes = $arrInGenes1 + $arrInGenes2;

			$arrFoundEnsOrth = fnAnyEnsOrtholog($arrEnsGenes1 , $arrEnsGenes2);
			if ($arrFoundEnsOrth === false) {
				continue;
			}
			//found ortholog
			echo("Found ensembl orthologs linking families $sGroup1 and $sGroup2: " );
			echo(implode(',' , $arrFoundEnsOrth ) . "\n");

			//examine in-paralogs:
			$sGroupID = "$sGroup1-$sGroup2";
			$sGeneSymbol = $arrGroup1['genesymbol'];
			$sGeneDesc = $arrGroup1['genedesc'];

			if (count($arrInTaxa)==count($arrSpecies)) { 
				echo("No in-paralogs found. $sGeneSymbol\n");

				fwrite($hOut, "$sGroupID\t$sGeneSymbol\t$sGeneDesc");
				foreach($arrSpecies as $sSp) {
					fwrite($hOut, "\t".$arrInGenes[$sSp][0]);
				}
				fwrite($hOut, "\n");
				$nGroupCount++;
			} else {
				echo("In paralogs found: ". implode("," , $arrIngroupTaxa1) .";" . implode("," , $arrIngroupTaxa2) . "\n");
				$arrInGenes = array();
				foreach($arrInGenes1 as $sTaxon => $arrGeneId) {
					if (!array_key_exists($sTaxon ,$arrInGenes2 )) {
						$arrInGenes[$sTaxon] = $arrGeneId;
					} else { //compare the two genes
						$sGene1Id = $arrInGenes1[$sTaxon][0];
						$sGene2Id = $arrInGenes2[$sTaxon][0];
						$nCompleteness1 = $arrGeneFilters[$sTaxon][$sGene1Id];
						$nCompleteness2 = $arrGeneFilters[$sTaxon][$sGene2Id];
						echo("$sTaxon : $sGene1Id $sGeneSymbol ($nCompleteness1 complete), $sGene2Id ($nCompleteness2 complete) ");
						if ($nCompleteness1 > $nCompleteness2) {
							$arrInGenes[$sTaxon] = array($sGene1Id);
							echo(" keep $sGene1Id\n");
						} else {
							$arrInGenes[$sTaxon] = array($sGene2Id);
							echo(" keep $sGene2Id\n");
						}
					}
				}

				$arrLeftOverTaxa = array_diff($arrSpecies , array_keys($arrInGenes));
				foreach($arrLeftOverTaxa as $sLeftTaxa ) {
					$arrInGenes[$sLeftTaxa] = $arrInGenes2[$sLeftTaxa];
				}
				fwrite($hOut, "$sGroupID\t$sGeneSymbol\t$sGeneDesc");
				foreach($arrSpecies as $sSp) {
					fwrite($hOut, "\t".$arrInGenes[$sSp][0]);
				}
				fwrite($hOut, "\n");
				$nGroupCount++;


			}

		}
	}

}

echo("In total $nGroupCount groups found after merging\n");
function fnParseGenes($s) {
	$arrF1 = explode(',' , $s);
	$arrRet = array(); //key is taxon name, value is array with gene ids
	foreach($arrF1 as $sGene) {
		if (substr($sGene, 0, 6) == '#trees') continue;
		$arrF2 = explode("|" , $sGene);
		$sTaxon = trim($arrF2[0]);
		$sGeneId = trim($arrF2[1]);
		if ( !fnIsGeneIncluded($sTaxon , $sGeneId) ) continue;

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

function fnLoadGeneFilters($arrWhiteLists) {
	$arrRet = array();

	foreach($arrWhiteLists as $sSp => $sGFF) {
		$arrRet[$sSp] = array();
		$hGFF = fopen($sGFF, 'r');
		while( false !== ($sLn = fgets($hGFF) )) {
			$sLn = trim($sLn);
			$arrF = explode("\t" , $sLn);
			$nGeneID = $arrF[0];
			$arrRet[$sSp][$nGeneID] = floatval($arrF[2]);
		}
	}

	return $arrRet;
}

function fnIsGeneIncluded($sSp, $nGeneID) {
	global $arrGeneFilters, $arrSpecies;
	if (!array_key_exists($sSp , $arrGeneFilters)) return true;
	return array_key_exists($nGeneID , $arrGeneFilters[$sSp]);
}

function fnReadEnsemblOrthologMap() {
	global $arrEnsemblOrthologIDMaps;

	$arrRet = array();

	$nMapFile = 0;
	foreach( $arrEnsemblOrthologIDMaps as $sMapFile ) {
		$nMapFile++;

		$hMap = fopen($sMapFile , 'r');
		$nOrthCount = 0;
		while(false !== ($sLn = fgets($hMap) )) {
			$nOrthCount++;
			$sLn = trim($sLn);
			$arrF = explode("\t", $sLn);

			$sGeneFamName = $nMapFile ."_" . $nOrthCount ;//$arrF[0]; // this is the zebrafish gene name, used as the family name
			foreach( array_slice($arrF, 1 ) as $sEnsID) {
				$arrRet[$sEnsID] = $sGeneFamName;
			}

		}
	}

	return $arrRet;
}

function fnAnyEnsOrtholog($arrEnsGenes1 , $arrEnsGenes2) {
	global $arrEnsOrthMap;

	$arrMap = array();

	foreach($arrEnsGenes1 as $sID) {
		if (array_key_exists($sID , $arrEnsOrthMap) ) {
			$sGeneFam = $arrEnsOrthMap[$sID];
			if (!array_key_exists($sGeneFam , $arrMap) ) {
				$arrMap[$sGeneFam] = array();
			}
			$arrMap[$sGeneFam][] = $sID;
		}
	}

	$bFoundOrth = false;
	$sFoundGeneFam = '';
	foreach($arrEnsGenes2 as $sID) {
		if (array_key_exists($sID , $arrEnsOrthMap) ) {
			$sGeneFam = $arrEnsOrthMap[$sID];
			if (!array_key_exists($sGeneFam , $arrMap) ) {
				continue;
			} else {
				$arrMap[$sGeneFam][] = $sID;
				$bFoundOrth = true;
				$sFoundGeneFam = $sGeneFam;
			}
		}
	}

	if ($bFoundOrth) {
		return $arrMap[$sFoundGeneFam];
	} else {
		return false;
	}

}

function sub_array($haystack, $needle)
{
    return array_intersect_key($haystack, array_flip($needle));
}

?>
