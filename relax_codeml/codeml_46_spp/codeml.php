<?php
require_once(dirname(__FILE__) . "/lib.php");

/*
Test for monophyly constraints present in monophyly_desc.txt
branch marking styles:
hyphy relax: target NodeName{T}, or (A{T},B{T}){T} . reference {R}, do not use {U}
Codeml style: marks all child branches as foreground: $1 (not used in this script) , marks only the one branch as foreground: #1 
*/

$sPremadeFasta = "/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/export_protein_multiref_23_july_2017/ret_improved_cleanDNA.fasta";
$sTreeFolder = "../export_protein_multiref_23_july_2017/Codeml_Nothobranchius";
$sOutDIR = "Codeml_Nothobranchius";
$sGeneSymbolMap = "/beegfs/group_dv/home/RCui/killifish_genomes/annotation/UPhO_final/assigngenesymbol/killi_orthologs_table.ensemblinformed.txt";

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
	global $sOutDIR, $nThisPart, $nTotalPart, $nGeneCount,$hOut,$arrComputedResults,$sTreeFolder, $arrGeneSymbolMap;

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

			echo(" ... Doing Codeml... ");

			
			if ($sFixTree=="") {
				echo("$sGeneName : User needs to specify tree!\n");
				return;
			}

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



?>
