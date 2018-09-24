<?php
//assign supplement current gff with additional gene symbol annotations.
//only changes the gene annotation in case it is "unkown".
$arrOrthoGroups = array("/beegfs/group_dv/home/RCui/killifish_genomes/annotation/UPhO_final/assigngenesymbol/UPhO_orthogroups.genesymbol.txt", //these gene symbols are assigned from many teleosts including zebrafish, gar etc. suffers from LBA sometimes.
                        "/beegfs/group_dv/home/RCui/killifish_genomes/annotation/UPhO_Cyprinodontiformes/assigngenesymbol/UPhO_orthogroups.genesymbol.txt"); //these gene symbols are assigned from Xmac and Molly
$sGFF = "maker.replace_frag_aln.noolverlap.gff";
$sOut = "maker.genesymbol.supple.gff";

$sThisSp = 'NOR';
$sNum2ID = "../../../../UPhO_final/$sThisSp"."_withid.fst";

$hOut = fopen($sOut , 'w');

echo("$sNum2ID \n");

$arrNum2ID  = fnLoadNum2ID($sNum2ID );
$arrID2Num = array_flip($arrNum2ID);

//echo($arrID2Num['NORv12scf1-frag-merged-region2-clust1']);

//print_r($arrID2Num);

$arrGeneSymbolCounts = array();

$arrGene2RNAMap = array();
$hGFF = fopen($sGFF, 'r');

$arrAvailableRNAIDs = array();
$arrRNAIDsWithSymbols = array();
$arrRNAIDsWithSymbols2 = array();

	while( false !== ($sLn = fgets($hGFF ) ) ) {
		if ($sLn == '#') continue;
		$sLn = trim($sLn);
		if ($sLn == '') continue;

		$arrF = explode("\t" , $sLn);

		if (count($arrF) != 9) continue;
		
		if ($arrF[2] == 'mRNA') {
			$arrAnnot = fnParseFields($arrF[8]);
			if (!array_key_exists("ID" , $arrAnnot) ) continue;
			if (!array_key_exists("Parent" , $arrAnnot) ) continue;
			$arrGene2RNAMap[$arrAnnot['Parent'] ] = $arrAnnot['ID'];
			$sRNAID =str_replace(".", "_", $arrAnnot['ID']);
			$arrAvailableRNAIDs[$sRNAID] = true;

			continue;
		}

	}
fclose($hGFF);

$arrID2Symbol = fnLoadOrthoGroups($arrOrthoGroups);

//print_r($arrID2Symbol);


$hGFF = fopen($sGFF, 'r');

	while( false !== ($sLn = fgets($hGFF ) ) ) {
		if ($sLn == '#') {
			fwrite( $hOut , $sLn);
			continue;
		}

		$sLn = trim($sLn);
		if ($sLn == '') continue;

		$arrF = explode("\t" , $sLn);

		if (count($arrF) != 9) {
			fwrite( $hOut , $sLn."\n");
			continue;
		}
		
		if ($arrF[2] == 'gene' || $arrF[2] == 'mRNA') {
			$arrAnnot = fnParseFields($arrF[8]);
			if (!array_key_exists("ID" , $arrAnnot) ) continue;
			if ($arrF[2] == 'gene') {
				if (!array_key_exists($arrAnnot['ID'] , $arrGene2RNAMap) ) {
					fwrite( $hOut , $sLn."\n");
					continue;
				}
			}
			$sRNAID = ($arrF[2] == 'gene')? $arrGene2RNAMap[$arrAnnot['ID']] : $arrAnnot['ID'];
			$sRNAID =str_replace(".", "_", $sRNAID);
		

			$arrAnnot['spgeneid'] = $arrID2Num[$sRNAID];

			if (!array_key_exists($sRNAID , $arrID2Symbol) ) {
				$arrAnnot['gene'] = 'unkown';
				$arrAnnot['description'] = 'unkown';
			} else {
				$arrAnnot['gene'] = $arrID2Symbol[$sRNAID]['genesymbol'];
				$arrAnnot['description'] = $arrID2Symbol[$sRNAID]['description'];
				$arrAnnot['ensemblorthologs'] = $arrID2Symbol[$sRNAID]['ensembl_orthologs'];
				if ($arrF[2] == 'mRNA') echo("Updated $sRNAID to ".$arrAnnot['gene'] ."\n");

			}
			$arrF[8] = fnArr2Annot($arrAnnot);
			fwrite( $hOut , implode("\t", $arrF)."\n");
			continue;
		} else {
			fwrite( $hOut , $sLn."\n");
			continue;
		}



	}
fclose($hGFF);



function fnLoadOrthoGroups($arrOrthoGroups) {
	global $sThisSp,  $arrNum2ID , $arrGeneSymbolCounts, $arrAvailableRNAIDs, $arrRNAIDsWithSymbols, $arrRNAIDsWithSymbols2;

	$arrRet = array(); //key is gene id
	$arrGeneSymbolCounts = array();

	foreach( $arrOrthoGroups as $sOrthoGroup) {
		echo("Loading ortholog definition $sOrthoGroup\n");
		$hF = fopen($sOrthoGroup , 'r');


		while( false !== ($sLn = fgets($hF) )) {
			$sLn = trim($sLn);
			if ($sLn == '') continue;
			list($sGroupID, $sGeneSymbol, $sDescription, $sGenes) = explode("\t", $sLn); 
			$arrGenes = explode(',' , $sGenes);

			$arrEnsemblOrthologs = array();
			$arrAddedGeneIDs = array();
			foreach( $arrGenes as $sGene) {
				if (substr($sGene,0,7) == '#trees_') continue;
				list($sSp, $sNum, $sGeneID) = explode('|' , $sGene);
				//echo("geneid:".$sGeneID."\n");
	 

				if (substr($sGeneID , 0, 3) == 'ENS' ) {
					$arrEnsemblOrthologs[] = preg_replace("/_\d+/" , '', $sGeneID);
				}



				if ($sSp == $sThisSp) {

					if (array_key_exists($sNum, $arrNum2ID )) {
						$sGeneID = $arrNum2ID[$sNum];
						//echo("$sGeneID\n");
					} 

					$sGeneID =str_replace(".", "_", $sGeneID);

					$arrAddedGeneIDs[] = $sGeneID;
					if (!array_key_exists($sGeneID , $arrRet)) {
						$arrRet[$sGeneID] = array();
					}

					if ( (!array_key_exists('genesymbol' , $arrRet[$sGeneID])) || ($arrRet[$sGeneID]['genesymbol'] == 'unknown')  ) {
						$arrRet[$sGeneID]['internalspeciesgeneid'] = $sNum;
						$arrRet[$sGeneID]['genesymbol'] = $sGeneSymbol;
						$arrRet[$sGeneID]['description'] = $sDescription;
						if (!array_key_exists($sGeneSymbol , $arrGeneSymbolCounts)) {
							$arrGeneSymbolCounts[$sGeneSymbol] = 0;
						}
						if (array_key_exists($sGeneID , $arrAvailableRNAIDs) && (!array_key_exists($sGeneID, $arrRNAIDsWithSymbols) )) {
							$arrGeneSymbolCounts[$sGeneSymbol] += 1;
							$arrRNAIDsWithSymbols[$sGeneID] = true;
						}

					}


				}
			}


			foreach($arrAddedGeneIDs as $sGeneID) {
				$arrRet[$sGeneID]['ensembl_orthologs'] = implode(',' , $arrEnsemblOrthologs);
			}
		}


	} //end foreach file

	$arrGeneSymbolCounter = array();
	foreach($arrRet as $sGeneID => &$arrInfo) {
		if ($arrGeneSymbolCounts[$arrInfo['genesymbol']] > 1) { //
			if (!array_key_exists($arrInfo['genesymbol'] , $arrGeneSymbolCounter) ) {
				$arrGeneSymbolCounter[$arrInfo['genesymbol']] = 0;
			}
			if (array_key_exists($sGeneID , $arrAvailableRNAIDs) && (!array_key_exists($sGeneID, $arrRNAIDsWithSymbols2))) {
				$arrGeneSymbolCounter[$arrInfo['genesymbol']]  += 1;
				$arrRNAIDsWithSymbols2[$sGeneID] = true;
			}
			$arrInfo['genesymbol'] = $arrInfo['genesymbol']." (".$arrGeneSymbolCounter[$arrInfo['genesymbol']]." of ".$arrGeneSymbolCounts[$arrInfo['genesymbol']].")";
		}
	}


	return $arrRet;
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

function fnArr2Annot($arr) {
	$s = "";
	
	foreach($arr as $sKey => $sVal) {
		if ($sVal === false) $sVal=0;
		$s .= "$sKey=$sVal;";
	}

	return substr($s, 0, strlen($s)-1);
}

function fnLoadNum2ID($sNum2ID ) {
	$hF = popen("grep '>' $sNum2ID" , 'r');
	$arrRet = array();
	while(false !== ($sLn = fgets($hF) )) {
		$arrF = explode('|', trim($sLn) );
		if (count($arrF)!=3) continue;
		$arrRet[$arrF[1]] = $arrF[2];
	}

	return $arrRet;
}


?>
