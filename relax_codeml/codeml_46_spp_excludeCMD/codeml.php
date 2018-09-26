<?php
require_once(dirname(__FILE__) . "/lib.php");

/*
This version excludes all codons with multiple substitutions
*/

$sPremadeFasta = "..//export_protein_multiref//ret_improved_cleanDNA.fasta";
$sTreeFolder = "../export_protein_multiref_23_july_2017/Codeml_Nothobranchius";
$sOutDIR = "Codeml_Nothobranchius";
$sGeneSymbolMap = "../../annotation/UPhO_final/assigngenesymbol/killi_orthologs_table.ensemblinformed.txt";

$sMonoDesc = "../export_protein_multiref/monophyly_force_callo_nothos_monophyly3.txt";
$arrCountGenera = array("Callopanchax", "Nothobranchius", "Scriptaphyosemion", "Aphyosemion", "Fundulopanchax", "Archiaphyosemion", "Epiplatys");


$nThisPart = 0;
$nTotalPart = 1;

while(count($argv) > 0) {
    $arg = array_shift($argv);
    switch($arg) {
        case '-p':
            $sPremadeFasta = trim(array_shift($argv));
            break;
        case '-o':
            $sOutDIR  = trim(array_shift($argv));
            break;
        case '-T':
            $sTreeFolder  = trim(array_shift($argv));
            break;
        case '-S':
            $sGeneSymbolMap  = trim(array_shift($argv));
            break;
        case '-N':
            $nTotalPart  = trim(array_shift($argv));
            break;
        case '-f':
            $nThisPart  = trim(array_shift($argv));
            break;

    }
}



exec("mkdir -p $sOutDIR");

$arrGeneSymbolMap = fnLoadGeneSymbolMap($sGeneSymbolMap);

$hPremadeFasta = fopen($sPremadeFasta , "r");
$sLn = "";

$sSeqName = "";
$sSeq = "";
$arrFullCodingSeq = array();
$sCurrentGene = "";
$nGeneCount = -1;

$sOutFile = "$sOutDIR/ret_part$nThisPart.of.$nTotalPart.txt";
$arrComputedResults = fnLoadCurrentRet($sOutFile ); //load previous results, and redo items that failed or not done yet.
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
	global $sOutDIR, $nThisPart, $nTotalPart, $nGeneCount,$hOut,$arrComputedResults,$sTreeFolder, $arrGeneSymbolMap, $sMonoDesc, $arrCountGenera;

	$nGeneCount++;
	if ($nGeneCount % $nTotalPart != $nThisPart) return;

		$sFixTreeFile = $sTreeFolder."/$sGeneName.txt";

		if (!file_exists($sFixTreeFile) ) {
			echo("$sGeneName is filtered out\n");
			return;
		}		


	$arrFullCodingSeq = fnExcludeLastStop($arrFullCodingSeq);
	$arrTaxa = array_keys($arrFullCodingSeq);

	if (array_key_exists($sGeneName , $arrComputedResults)) {
		echo("Prev result available:\n");
		echo(implode("\t", $arrComputedResults[$sGeneName]). "\n");
		echo("Curr: " . strlen(trim($arrFullCodingSeq[$arrTaxa[0]]))/3 . "\n");
	}

	if (array_key_exists($sGeneName , $arrComputedResults) && $arrComputedResults[$sGeneName][1] == 'Success' && trim($arrComputedResults[$sGeneName][5]) != '' && trim($arrComputedResults[$sGeneName][6]) != ''  && ($arrComputedResults[$sGeneName][3]*3) == strlen(trim($arrFullCodingSeq[$arrTaxa[0]])) ) {
		$sCodemlTmp = $arrComputedResults[$sGeneName][7];
		echo("Checking previous run $sCodemlTmp \n");
		$sNex = "$sCodemlTmp/in.fas";
		if (file_exists($sNex)) {
			$sNexContent = file_get_contents($sNex);
			
				$nPrevTax = substr_count($sNexContent , '>'); 
				if ($nPrevTax == count($arrTaxa) ) { //seems promising, check each taxon
					$bAllPresent = true;
					foreach($arrTaxa as $sTaxon) {
						if (strpos($sNexContent, ">".$sTaxon) === false) {
							$bAllPresent = false;
							break;
						}
					}

					if ($bAllPresent) {

						echo("$sGeneName done. Skip.\n");
						fwrite($hOut, implode("\t", $arrComputedResults[$sGeneName]). "\n");
						return;
					}

				}
			
		}

	} 
		$arrAlnRet = fnTranslateAlignment($arrFullCodingSeq , false); //exclude last stop if exists.


		$sFixTree = str_replace( ';','',trim( file_get_contents($sFixTreeFile))); 
		if (!$arrAlnRet["stopcodon"] ) { //there is no stop codon, do codeml

			echo(" ... Doing Codeml... \n");

			
			if ($sFixTree=="") {
				echo("$sGeneName : User needs to specify tree!\n");
				return;
			}

			//exclude codons with multiple substitutions
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
								//echo("$nPos $sSp1 : $sSp2 = $nSubType\n");
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
				echo("\nNo site needs exclusion, skip\n");
				return;
			}

			//print_r($arrFullCodingSeq);

			foreach($arrExcludeSites as $nPos => $bDump) {
				foreach($arrFullCodingSeq as $sTaxon => &$sSeq) {
					$arrFullCodingSeq[$sTaxon] = substr_replace( $sSeq , 'NNN' ,$nPos, 3); //mask
				}
			}

			//print_r($arrFullCodingSeq);
			//die();

			$sTree = $sFixTree;

			$sGeneSymbol = array_key_exists($sGeneName, $arrGeneSymbolMap)? $arrGeneSymbolMap[$sGeneName][0] : "-";
			$sGeneFullName = array_key_exists($sGeneName, $arrGeneSymbolMap)? $arrGeneSymbolMap[$sGeneName][1] : "-";

			$oCodeml = new CodemlUtil();
			$nSeqLen = 0;
			foreach($arrFullCodingSeq as $sKey => $sSeq) {
				$oCodeml->AddSequence($sKey , $sSeq);
				$nSeqLen = strlen( $sSeq);
			}
			$oCodeml->SpecifyTree($sTree);
			$arrCodemlRet = $oCodeml->RunCodemlBranchSite( true); //only need to get dn ds
			if (!$arrCodemlRet["Success"]) {
				fwrite($hOutCodeml, "$sGeneName\tCodeml error:".$arrCodemlRet["tmpfolder"].PHP_EOL);
				continue;
			}
			
			
			fwrite($hOut, "$sGeneName\tSuccess\t$sGeneSymbol\t".($nSeqLen/3)
				."\t".$arrCodemlRet["MA-MAnull"]
				."\t".$arrCodemlRet["MA"]
				."\t".$arrCodemlRet["MA_null"]
				."\t".$arrCodemlRet["tmpfolder"]
				."\t$sGeneFullName"
			.PHP_EOL);
			
			
		}
		else {
			fwrite($hOut, "$sGeneName\tStopCodonFoundIn:".implode(",",$arrAlnRet["stoptaxa"]).PHP_EOL);
		}


	
}




function fnLoadCurrentRet($sOutFile ) {

	$arrRet = array();

	if (!file_exists($sOutFile) ) {
		return $arrRet;
	}

	$hPrevOut = fopen($sOutFile , 'r');

	while( false !== ($sLn = fgets($hPrevOut) ) ) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		$arrF = explode("\t", $sLn);
		$arrRet[$arrF[0]] = $arrF;
	}
	fclose($hPrevOut);

//	$hOut = fopen($sOutFile , 'w');

	foreach($arrRet as $sGroup => $arrLn) {
		if ($arrLn[1] == 'Success' && trim($arrLn[5]) != '' && trim($arrLn[6]) != '' ) {
//			fwrite($hOut , implode("\t" , $arrLn) . "\n"); //only write back good ones.
		} else {
			echo("To be retried: $sGroup \n");
		}
	}

//	fclose($hOut);

	return $arrRet;
}

function fnLoadGeneSymbolMap($sOrthologList) {

	if (!file_exists($sOrthologList) ) {
		die("Gene symbol definition not found $sOrthologList\n");
	}

	$hOrthologList = fopen($sOrthologList , "r");
	$arrOrthologMeta = array();


	$nLn = -1;
	echo("Parsing ortholog definitions...\n");
	while( ($sLn=fgets($hOrthologList))!==false ) {
		$sLn = trim($sLn);
		if ($sLn == "") {
			continue;
		}
		$nLn++;
		if ($nLn==0) {
			continue; // skip header
		}
	
		$arrFields = explode("\t", $sLn);
		$arrOrthologMeta[$arrFields[0]] = array_slice($arrFields , 1, 2);
	}

	echo("Loaded ". count($arrOrthologMeta) ." ortholog definitions\n" );

	return $arrOrthologMeta;
}

function str_lreplace($search, $replace, $subject)
{
    $pos = strrpos($subject, $search);

    if($pos !== false)
    {
        $subject = substr_replace($subject, $replace, $pos, strlen($search));
    }

    return $subject;
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

?>
