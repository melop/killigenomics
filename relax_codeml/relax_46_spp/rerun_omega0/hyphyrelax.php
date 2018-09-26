<?php
require_once(dirname(__FILE__) . "/lib.php");

/*
Test for monophyly constraints present in monophyly_desc.txt
branch marking styles:
hyphy relax: target NodeName{T}, or (A{T},B{T}){T} . reference {R}, do not use {U}
Codeml style: marks all child branches as foreground: $1 (not used in this script) , marks only the one branch as foreground: #1 

This program uses Hyphy 2.3.8, which has a better convergence behavior
It identifies the problematic cases from the Hyphy 2.2.5 runs and rerun them up to 20 times, to see if a solution with 0<omega0<1 with a better ML score can be identified.
*/

$sPremadeFasta = "/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/export_protein_multiref_23_july_2017/ret_improved_cleanDNA.fasta";
$sTreeFolder = "../../../export_protein_multiref_23_july_2017/Relax_Nothobranchius";
$sOutDIR = "Relax_Nothobranchius";

$sGeneSymbolMap = "/beegfs/group_dv/home/RCui/killifish_genomes/annotation/UPhO_final/assigngenesymbol/killi_orthologs_table.ensemblinformed.txt";

$nMaxTries = 20;

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

$sOldOutDIR = "../$sOutDIR";

exec("mkdir -p $sOutDIR");

$arrGeneSymbolMap = fnLoadGeneSymbolMap($sGeneSymbolMap);

$hPremadeFasta = fopen($sPremadeFasta , "r");
$sLn = "";

$sSeqName = "";
$sSeq = "";
$arrFullCodingSeq = array();
$sCurrentGene = "";
$nGeneCount = -1;

$sOutFile = "$sOldOutDIR/ret_part*.txt";
$arrComputedResults = fnLoadCurrentRet($sOutFile ); //load previous results, and redo items that failed or not done yet.
//print_r($arrComputedResults);

$sNewOutFile = "$sOutDIR/ret_part$nThisPart.of.$nTotalPart.txt";
$hOut = fopen($sNewOutFile , 'w'); //append to it

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
	global $sOutDIR, $nThisPart, $nTotalPart, $nGeneCount,$hOut,$arrComputedResults,$sTreeFolder, $arrGeneSymbolMap, $nMaxTries ;

	$nGeneCount++;
	if ($nGeneCount % $nTotalPart != $nThisPart) return;

		$sFixTreeFile = $sTreeFolder."/$sGeneName.txt";

		if (!file_exists($sFixTreeFile) ) {
			echo("$sGeneName is filtered out\n");
			return;
		}	

	$arrFullCodingSeq = fnExcludeLastStop($arrFullCodingSeq);
	$arrTaxa = array_keys($arrFullCodingSeq);

	//echo($arrComputedResults[$sGeneName][3] . " current algn: " .  strlen($arrFullCodingSeq[$arrTaxa[0]]) );
	//die();
	if (array_key_exists($sGeneName , $arrComputedResults)) {
		echo("Prev result available:\n");
		echo(implode("\t", $arrComputedResults[$sGeneName]). "\n");
		echo("Curr: " . strlen(trim($arrFullCodingSeq[$arrTaxa[0]]))/3 . "\n");
	} else {
		return;
	}

	$bProblemCondition = floatval($arrComputedResults[$sGeneName][4]) < 0.05 && 
		floatval($arrComputedResults[$sGeneName][8]) > 1 && 
		floatval($arrComputedResults[$sGeneName][10]) == 0 ; // these are the problematic cases, where omega0==0 and k > 1

	if (array_key_exists($sGeneName , $arrComputedResults) && 
		 $arrComputedResults[$sGeneName][1] == 'Success'  && 
		(!$bProblemCondition) && 
		($arrComputedResults[$sGeneName][3]*3) == strlen(trim($arrFullCodingSeq[$arrTaxa[0]])) ) {
		$sHyphyTmp = $arrComputedResults[$sGeneName][23];
		echo("$sGeneName done. Skip.\n");
		fwrite($hOut, implode("\t", $arrComputedResults[$sGeneName]). "\n");
		return;

	} 

	echo("Gene : $sGeneName is problematic, redo estimation with Hyphy 2.3.8 for up to $nMaxTries tries\n");
	$arrAlnRet = fnTranslateAlignment($arrFullCodingSeq , false); //exclude last stop if exists.

	if (!$arrAlnRet["stopcodon"] ) { //there is no stop codon, do codeml

		$oCurrBest = false;

		for($nTry=1;$nTry<=$nMaxTries;$nTry++) {
			echo("	===try $nTry ====\n");
			$sFixTree = '[&R] ' .   str_replace( ';','',trim( file_get_contents($sFixTreeFile))); 


			echo(" ... Doing RELAX... \n");

			
			if ($sFixTree=="") {
				echo("$sGeneName : User needs to specify tree!\n");
				return;
			}

			$sTree = $sFixTree;

			$oHyphy = new HyPhyUtil();
			$nSeqLen = 0;
			foreach($arrFullCodingSeq as $sKey => $sSeq) {
				$oHyphy->AddSequence($sKey , $sSeq);
				$nSeqLen = strlen( $sSeq);
			}
			$oHyphy->SpecifyTree($sTree);
			$arrHyphyRet = $oHyphy->RunRELAX( true); //only need to get dn ds
			//unset($oHyphy);
			if (!$arrHyphyRet["Success"]) {
				echo("$sGeneName\tRELAX error:".$arrHyphyRet["tmpfolder"].PHP_EOL);
				continue;
			}

			//check if results have converged:
			if ( floatval($arrHyphyRet["R_omega0"]) == 0 || floatval($arrHyphyRet["T_omega0"]) == 0) { //if parameter estimates are 0

				if ($oCurrBest !== false && floatval($arrHyphyRet["Alternative"]) ==  floatval($oCurrBest["Alternative"])  ) {
					echo("likelihood doesn't change between runs, give up\n");
					break;
				}

				if ($oCurrBest === false || floatval($arrHyphyRet["Alternative"]) >  floatval($oCurrBest["Alternative"]) ) {
					$oCurrBest = $arrHyphyRet;
				}


				continue;
			} else if (floatval($arrHyphyRet["R_omega0"]) > 0 && floatval($arrHyphyRet["T_omega0"]) > 0) {
				if ($oCurrBest === false || floatval($arrHyphyRet["Alternative"]) >  floatval($oCurrBest["Alternative"])) {
					$oCurrBest = $arrHyphyRet;
					echo("Better solution found at iteration $nTry\n");
					print_r($oCurrBest);
					break;
				}
			}
						


		} //end for

		if ($oCurrBest === false || (floatval($oCurrBest["R_omega0"])==0 || floatval($oCurrBest["T_omega0"])==0  ) ) {
			if (array_key_exists($sGeneName , $arrComputedResults)) {
				echo("No better solution found\nuse old results\n");
				fwrite($hOut, implode("\t", $arrComputedResults[$sGeneName]). "\n");
				return;
			}


		} else if ($oCurrBest !== false) {
			echo("\nreport the run with the best LR.\n");
		} else {
			echo("No better solution found\nNo old report available, skip.\n");
			return;
		}

			$arrHyphyRet = $oCurrBest;
			$sGeneSymbol = array_key_exists($sGeneName, $arrGeneSymbolMap)? $arrGeneSymbolMap[$sGeneName][0] : "-";
			$sGeneFullName = array_key_exists($sGeneName, $arrGeneSymbolMap)? $arrGeneSymbolMap[$sGeneName][1] : "-";

			fwrite($hOut, "$sGeneName\tSuccess\t$sGeneSymbol\t".($nSeqLen/3)
				."\t".$arrHyphyRet["P"] 
				."\t".$arrHyphyRet["MG94xREV"] 
				."\t".$arrHyphyRet["NULL"] 
				."\t".$arrHyphyRet["Alternative"] 
				."\t".$arrHyphyRet["K"] 
				."\t".$arrHyphyRet["LR"] 
				."\t".$arrHyphyRet["R_omega0"] 
				."\t".$arrHyphyRet["R_omega0_prop"] 
				."\t".$arrHyphyRet["R_omega1"] 
				."\t".$arrHyphyRet["R_omega1_prop"] 
				."\t".$arrHyphyRet["R_omega2"] 
				."\t".$arrHyphyRet["R_omega2_prop"] 
				."\t".$arrHyphyRet["T_omega0"] 
				."\t".$arrHyphyRet["T_omega0_prop"] 
				."\t".$arrHyphyRet["T_omega1"] 
				."\t".$arrHyphyRet["T_omega1_prop"] 
				."\t".$arrHyphyRet["T_omega2"] 
				."\t".$arrHyphyRet["T_omega2_prop"] 
				."\t$sGeneFullName"
				."\t".$arrHyphyRet["tmpfolder"]
			.PHP_EOL);


	}
	else {
			fwrite($hOut, "$sGeneName\tStopCodonFoundIn:".implode(",",$arrAlnRet["stoptaxa"]).PHP_EOL);
	}

	//die();
}




function fnLoadCurrentRet($sOutFile ) {

	$arrRet = array();

	if (count(glob($sOutFile)) == 0 ) {
		echo("Previous output $sOutFile not found, die\n ");
		die();
		return $arrRet;
	}

	$hPrevOut = popen("cat $sOutFile | awk '{if ($2 == \"Success\" ) print $0; }' " , 'r');

	while( false !== ($sLn = fgets($hPrevOut) ) ) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		$arrF = explode("\t", $sLn);
		$arrRet[$arrF[0]] = $arrF;
	}
	fclose($hPrevOut);

//	$hOut = fopen($sOutFile , 'w');

	foreach($arrRet as $sGroup => $arrLn) {
		if ($arrLn[1] == 'Success' && floatval($arrLn[4]) < 0.05 &&  floatval($arrLn[8]) == 0 ) {
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



?>
