<?php
require_once(dirname(__FILE__) . "/lib.php");
//go through the fasta alignments again,
//check the Amino acids at the candidate sites,
//make sure that the mutation is unique in the target genera Callopanchax and Nothobranchius, and it's not shared by the immediate non-annual outgroup
// sharing with semi-annual is allowed.
$nTotalSitesExamined = 0;

$sMonoDesc = "export_protein_multiref/monophyly_force_callo_nothos_monophyly3.txt";
$sPremadeFasta = "export_protein_multiref/ret_improved_cleanDNA.fasta";
$sOutDIR = "annotate_aa_changes_excludeCMD";

$arrOutgroupGenera = array("Scriptaphyosemion", "Aphyosemion", "Fundulopanchax", "Archiaphyosemion", "Epiplatys", "Fundulosoma");
$arrConvergenceGenera = array("Callopanchax", "Nothobranchius");

$arrCountGenera = array("Callopanchax", "Nothobranchius", "Scriptaphyosemion", "Aphyosemion", "Fundulopanchax", "Archiaphyosemion", "Epiplatys");


$nMinTargetTaxaWithDerivedAAPerGenusPerc = 1; // float, 0 - 1, in each convergence genus, at least how many percentage of species need to have this convergent AA?
$nMaxAAStatesInOutgroup = 1; //maximum of 1 AA state in outgroup.

exec("mkdir -p $sOutDIR");
$hPremadeFasta = fopen($sPremadeFasta , "r");
$sLn = "";

$sSeqName = "";
$sSeq = "";
$arrFullCodingSeq = array();
$sCurrentGene = "";
$nGeneCount = -1;

$nThisPart = 0;
$nTotalPart = 1;

$sOutFile = "$sOutDIR/ret_minintaxaPerc$nMinTargetTaxaWithDerivedAAPerGenusPerc.maxAAout$nMaxAAStatesInOutgroup.part$nThisPart.of.$nTotalPart.txt";

//print_r($arrComputedResults);

$hOut = fopen($sOutFile , 'w'); //append to it

do {
	$sLn = fgets($hPremadeFasta);

	if ($sLn==false) {
		fnProcessRecord($sSeqName, $sSeq);
		break;
	}

	$sLn = trim($sLn);
	if (substr($sLn , 0,1) == ">") {
		fnProcessRecord($sSeqName, $sSeq);
		$sSeqName = substr($sLn , 1);
		$sSeq = "";
		
	} else {
		$sSeq .= $sLn;
	}

} while (true);
fnProcessAln($sCurrentGene,$arrFullCodingSeq);

echo("Processed: $nTotalSitesExamined AAs\n");
function fnProcessRecord($sSeqName, $sSeq) {
	global $arrFullCodingSeq , $sCurrentGene;
	if (strpos($sSeqName , "direction") !== false  ) {
		return; //don't process fragments
	}

	preg_match("/Ortho:(\S+);Sp:(\S+);MappedToRef:(\S+);SpGeneId:(\S+);GroupId:(\S+)/",  $sSeqName, $arrParsed);
/*
array(6
0	=>	Ortho:0;Sp:AAU;MappedToRef:AAU;SpGeneId:029799;GroupId:Group_0_0
1	=>	0
2	=>	AAU
3	=>	AAU
4	=>	029799
5	=>	Group_0_0
)

*/

	if (count($arrParsed)!=6 ) {
		return;
	}
	$sGeneName = $arrParsed[5];
	$sSpecies = $arrParsed[2];

	if ($sCurrentGene == $sGeneName ) {
		if (!array_key_exists($sSpecies , $arrFullCodingSeq)) {
			$arrFullCodingSeq[$sSpecies] = "";
		}

		$arrFullCodingSeq[$sSpecies] .= $sSeq;
	} else {
		if (count($arrFullCodingSeq) > 0) {
			fnProcessAln($sCurrentGene,$arrFullCodingSeq);
		}
		$arrFullCodingSeq = array();
		$sCurrentGene = $sGeneName;
		$arrFullCodingSeq[$sSpecies] = $sSeq;
	}

}


function fnProcessAln($sGeneName,$arrFullCodingSeq) {
	global $sOutDIR, $sMonoDesc, $nThisPart, $nTotalPart, $nGeneCount,$hOut, $arrOutgroupGenera , $arrConvergenceGenera, $nMinTargetTaxaWithDerivedAAPerGenusPerc, $nMaxAAStatesInOutgroup,$nTotalSitesExamined, $sMonoDesc, $arrCountGenera;

	$nGeneCount++;
	if ($nGeneCount % $nTotalPart != $nThisPart) return;



		$arrFullCodingSeq = fnExcludeLastStop($arrFullCodingSeq);


		//Check if all sequences are identical, if so then no need to continue with codeml to save time
		$bDiff = false;
		$sTempSeq = -1;
		foreach($arrFullCodingSeq as $sTaxon => $sSeq) {
			if ($sTempSeq == -1) {
				$sTempSeq = $sSeq;
				continue;
			}
			if ($sTempSeq != $sSeq) {
				$bDiff = true;
				break;
			}
		}

		if (!$bDiff) {
			//fwrite($hOutTable, "$sGeneName\tAll sequences identical.".PHP_EOL);
			echo("$sGeneName : all sequences identical \n");
			fwrite($hOut , "$sGeneName\tAllSeqIdentical\n");
			return;
		}
			
		//remove taxa consist of only N's or gaps:
		$arrRemoveTaxa = array();
		foreach($arrFullCodingSeq as $sTaxon => $sSeq) {
			if ( preg_replace("/[N\-n]/", "", $sSeq) == '') {
				$arrRemoveTaxa[] = $sTaxon;
			}
		}

		foreach($arrRemoveTaxa as $sRemoveTaxon) {
			unset($arrFullCodingSeq[$sRemoveTaxon]);
		}


			$oConstraintTree = fnParseMonoDesc($sMonoDesc, array_flip(array_keys($arrFullCodingSeq)) );

			$arrCountGenSpp = array();
			foreach($arrCountGenera as $sGenus) {

				$arrCountGenSpp[$sGenus] = fnGetAllDescendantTips($oConstraintTree, $sGenus) ;
			}

			//print_r(array_keys($arrFullCodingSeq) );
			$arrExcludeSites = array();

			foreach($arrCountGenSpp as $sGenus => $arrSpp) {

				foreach($arrSpp as $nSp1 => $sSp1) {
					$nAlnLen = strlen($arrFullCodingSeq[$sSp1]);
					foreach($arrSpp as $nSp2 => $sSp2) {
						if ($nSp1 == $nSp2) continue;
						//compare two sequences codon by codon
						for($nPos=0;$nPos<$nAlnLen;$nPos+=3) {
							$sCodon1 = substr($arrFullCodingSeq[$sSp1], $nPos, 3);
							$sCodon2 = substr($arrFullCodingSeq[$sSp2], $nPos, 3);
							if (strlen($sCodon1) !=3 || strlen($sCodon2) !=3) continue;

							$nSubType = fnIsMultipleSub($sCodon1, $sCodon2);
							if ($nSubType > 1) {
								$arrExcludeSites[$nPos] = true;
							}
						}
					}

				}

			}

			//now exclude positions by masking them as N
			if (count($arrExcludeSites) > 0) {
				$arrExcludeTheseSites = array_keys($arrExcludeSites);
				sort($arrExcludeTheseSites);
				echo("\nWarning: sites excluded: ".implode(" ", $arrExcludeTheseSites) . "\n");
			} else {
				echo("\nNo site needs exclusion\n");
			}

			foreach($arrExcludeSites as $nPos => $bDump) {
				foreach($arrFullCodingSeq as $sTaxon => &$sSeq) {
					$arrFullCodingSeq[$sTaxon] = substr_replace( $sSeq , 'NNN' ,$nPos, 3); //mask
				}
			}



		$arrAlnRet = fnTranslateAlignment($arrFullCodingSeq , false); //exclude last stop if exists.

		if (!$arrAlnRet["stopcodon"] ) { //there is no stop codon, do raxml
			$arrTaxa = array_keys($arrFullCodingSeq);


			//$sConstraintTree = fnTree2Newick($oConstraintTree , array());
			echo(count($arrTaxa) . " taxa included for $sGeneName \n" );
			fnGetAllDescendantTips($oConstraintTree ,"Nothobranchius" );
			//now go over the alignment.
			$arrAAAln =  $arrAlnRet['proteins'];

			//print_r($arrAAAln);
			$nAlnLen = strlen($arrAAAln[$arrTaxa[0]]);
			//convert genera names to spp names

			$arrConvergentGenSpp = array();
			foreach($arrConvergenceGenera as $sGenus) {

				$arrConvergentGenSpp[$sGenus] = fnGetAllDescendantTips($oConstraintTree, $sGenus) ;
			}

			$arrOutgroupSpp = array();

			foreach($arrOutgroupGenera as $sGenus) {

				$arrOutgroupSpp = array_merge($arrOutgroupSpp, fnGetAllDescendantTips($oConstraintTree, $sGenus)) ;
			}



			for($nPos=0;$nPos<$nAlnLen;$nPos++) {
				$nTotalSitesExamined++;
				//find if there are shared AA between two genera:
				$arrCommonAA = array();
				$nTaxCount = 0;
				$arrTargetGeneraCounts = array();
				foreach($arrConvergenceGenera as $sGenus) {
					$arrAAcounts = fnCountAAForSpp($arrAAAln , $nPos , $arrConvergentGenSpp[$sGenus]);
					$arrTargetGeneraCounts[$sGenus] = $arrAAcounts;
					$nTaxCount++;
					if ($nTaxCount == 1) {
						$arrCommonAA = array_keys($arrAAcounts);
						continue;
					}

					$arrCommonAA = array_values(array_intersect($arrCommonAA  , array_keys($arrAAcounts) ));
				}

				if (count($arrCommonAA) == 0 ) { // no common AA between the convergent taxa, continue
					continue;
				}


				$arrAAOutgroup =  fnCountAAForSpp($arrAAAln , $nPos , $arrOutgroupSpp);

				if (count($arrAAOutgroup) > $nMaxAAStatesInOutgroup) {
					continue;
				}

				//make sure that the common AA are not found in the outgroup:
				foreach($arrCommonAA as $sAA) {
					if (array_key_exists($sAA , $arrAAOutgroup) ) {
						continue; //this aa is also found in the outgroups, skip
					}

					//check $nMinTargetTaxaWithDerivedAAPerGenusPerc
					foreach($arrTargetGeneraCounts as $sGenus => &$arrCounts) {
						if ( ($arrCounts[$sAA] / array_sum($arrCounts)) < $nMinTargetTaxaWithDerivedAAPerGenusPerc) {
							continue 2;
						}
					}

					fwrite($hOut , "$sGeneName\t$nPos\t$sAA\t". fnArray2Str($arrTargetGeneraCounts) ."\t".http_build_query($arrAAOutgroup,'',',')."\n");
					
				}
			}

			//die();
		}
		else {
			echo("$sGeneName : internal stop codon found, skip \n");
			fwrite($hOut , "$sGeneName\tStopFound\n");
			return;
		}
		


	
}

function fnArray2Str($arrTargetGeneraCounts) {
	$s = "";
	foreach($arrTargetGeneraCounts as $sGenus => &$arrCounts) {
		$s .= "$sGenus:".http_build_query( $arrCounts, '',',').";";
	}

	return $s;
}

function fnCountAAForSpp(&$arrProtAln, $nPos, &$arrSpp) {

	$arrAACounts = array();

	foreach($arrSpp as &$sSp) {
		$sAA = $arrProtAln[$sSp][$nPos];
		if ($sAA == '?' || $sAA == 'X') continue; //do not count ambiguous sites
		if (!array_key_exists($sAA, $arrAACounts)) {
			$arrAACounts[$sAA] = 0;
		}

		$arrAACounts[$sAA]++;
	}

	return $arrAACounts;
}


function fnGetAllDescendantTips($oTree, $sTaxon) {

  $arrTips = array();


  foreach( $oTree as $sName => $arrChild) {
    if ($sName == $sTaxon) {
      #cat(sName , " ", sTaxon, "\n");
      if ( $arrChild === true ) {
        return array($sName);
      } else {
        #get all children
        foreach(array_keys($arrChild) as $sChildName ) {
          $arrTips = array_merge($arrTips, fnGetAllDescendantTips($arrChild , $sChildName ));
        }
        return $arrTips;
      }
    } else { #search in lower level
      #get all children
      if ($arrChild !== true) {
      		$arrTips = array_merge( $arrTips , fnGetAllDescendantTips($arrChild , $sTaxon) );
	}
    }
  }
  
  return $arrTips;

}

function fnIsMultipleSub($sCodon1, $sCodon2) {
	
	$sCodon1 = strtoupper($sCodon1);
	$sCodon2 = strtoupper($sCodon2);
	$arrDiffSites = array();
	for($i=0;$i<3;$i++) {
		if ($sCodon1[$i] == 'N' || $sCodon1[$i] == '-' || $sCodon2[$i] == 'N' || $sCodon2[$i] == '-' ) {
			continue;
		}

		if ($sCodon1[$i] != $sCodon2[$i] ) {
			$arrDiffSites[] = $i;
		}
	}

	if (count($arrDiffSites) == 0 ) return 0; 

	if (count($arrDiffSites) == 1 ) return 1; 

	if (count($arrDiffSites) == 3) return 3;

	return ( (abs($arrDiffSites[1] - $arrDiffSites[0]) == 1)? 4:2) ; //is adjacent or not?
}



function fnParseMonoDesc($sIn, $arrIncludeTaxon = true) {
	$hIn = fopen($sIn, 'r');

	$arrMap = array();
	while( false !== ($sLn = fgets($hIn) ) ) {
		$sLn = trim($sLn);
		$arrF = explode("\t" , $sLn);
		if (count($arrF) !=2) continue;
		$arrMap[$arrF[0]] = explode(",", $arrF[1]);
	}

	//print_r($arrMap);
	//now build tree, starting from root:
	$arrTree = fnBuildTree('root', $arrMap , $arrIncludeTaxon);

	return $arrTree;
}

function fnBuildTree($sNodeName, &$arrMap, &$arrIncludeTaxon) {

	$arrTree = array();

	if (!array_key_exists($sNodeName , $arrMap) ) {
		die("Node name $sNodeName undefined in map\n");
	}

	$arrChilds = $arrMap[$sNodeName];

	
	foreach($arrChilds as $sChild) {
		$sChild = trim($sChild);
		if ( $sChild[0] == '*' ) { //expand it
			$sChildNodeName = substr($sChild, 1); //delete * from taxon name
			//echo($sChildNodeName."\n");
			$oChildNode = fnBuildTree($sChildNodeName , $arrMap , $arrIncludeTaxon);
			if ( count($oChildNode) > 0 ) {
				$arrTree[$sChildNodeName] = $oChildNode;
			} 
		} else {
			if ($arrIncludeTaxon === true || array_key_exists($sChild, $arrIncludeTaxon) ) {
				$arrTree[$sChild] = true; //terminal taxon
			}
		}
	}

	return $arrTree;
}

function fnTree2Newick(&$oNode , $arrMarkBranches, $nStyle = 0, $bMark=false, $sMark="") { //arrMarkBranches : branch name => array( bMarkChild? , mark) . style 0 - turn off marking. style 1 - hyphy relax, 2 - codeml, $bMark and $sMark are used inernally, do not set
	$sOut = (count($oNode)>1) ? "(" : "";
	$nCount = -1;

	//print_r($arrMarkBranches) ;
	foreach($oNode as $sNodeName => $arrNodeContent) {
		$nCount++;
		if ($nCount > 0) {
			$sOut .= ",";
		}

		$bMarkChild = false;
		$bMarkNode = (array_key_exists($sNodeName , $arrMarkBranches) || $bMark) ;


		if (!$bMark) { //check if instructed to mark

			if ($bMarkNode) {
				$sMark = ($nStyle==1)? '{'.$arrMarkBranches[$sNodeName][1].'}' : ( ($nStyle == 2)? '#'.$arrMarkBranches[$sNodeName][1] : '');
			} 

			//echo("$sNodeName style $nStyle, mark: $sMark\n");
		} else {
			$bMarkChild = true;
		}

		if (array_key_exists($sNodeName , $arrMarkBranches) && $arrMarkBranches[$sNodeName][0] ) {
			$bMarkChild = true;
		}

		if ($arrNodeContent === true) { // terminal node
			$sOut .= $sNodeName;
		} else {
			$sOut .= fnTree2Newick( $arrNodeContent , $arrMarkBranches, $nStyle, $bMarkChild, $sMark );
		}

		if ($bMarkNode) {
			$sOut .= (count($arrNodeContent)>1 || $arrNodeContent === true)? $sMark:'';
		}

		if ( (!$bMarkNode) && $nStyle == 1) {
			$sOut .= (count($arrNodeContent)>1 || $arrNodeContent === true)? '{R}':'';
		}


	}

	$sOut .= (count($oNode)>1) ? ")" : "";

	return $sOut;
}

?>
