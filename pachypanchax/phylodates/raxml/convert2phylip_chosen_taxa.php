<?php
ini_set("memory_limit", "-1");

require_once(dirname(__FILE__) . "/lib.php");

/*
Test for monophyly constraints present in monophyly_desc.txt
branch marking styles:
hyphy relax: target NodeName{T}, or (A{T},B{T}){T} . reference {R}, do not use {U}
Codeml style: marks all child branches as foreground: $1 (not used in this script) , marks only the one branch as foreground: #1 
*/

$sPremadeFasta = "../ret_improved_cleanDNA.fasta"; //"exampleDNA.fasta"
$sOutDIR = "./phylipformat_PML_PLP_APL";
$arrChosenTaxa = array("PML", "PLP"); 

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

    }
}

$arrChosenTaxa = array_flip($arrChosenTaxa);

exec("mkdir -p $sOutDIR");


$hPremadeFasta = fopen($sPremadeFasta , "r");
$sLn = "";

$sSeqName = "";
$sSeq = "";
$arrFullCodingSeq = array();
$sCurrentGene = "";
$nGeneCount = -1;

$sOutFileFull = "$sOutDIR/full.phy";
//$sOutFile4Fold = "$sOutDIR/4fold.phy";
$sOutFile4Fold = "$sOutDIR/syn.phy";
$sOutFileCodon12 = "$sOutDIR/codon12.phy";

$hOutFull = fopen($sOutFileFull , 'w');
$hOut4fold = fopen($sOutFile4Fold , 'w');
$hOutCodon12 = fopen($sOutFileCodon12 , 'w');

$arrProccessedAln = array();
$arrTaxaList = array();

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

echo("Genes for each taxon:\n");
foreach($arrTaxaList as $sTaxon => $nGeneCount) {
	echo("$sTaxon\t$nGeneCount\n");
}

fnWritePhylip($arrProccessedAln['cleancols'], $hOutFull);
fnWritePhylip($arrProccessedAln['pos12'], $hOutCodon12);
fnWritePhylip($arrProccessedAln['4fold'], $hOut4fold);



function fnProcessRecord($sSeqName, $sSeq) {
	global $arrFullCodingSeq , $sCurrentGene, $arrChosenTaxa;
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

	if (!array_key_exists($sSpecies, $arrChosenTaxa) ) {
		return;
	}


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
	global $arrProccessedAln, $arrTaxaList;
	
	$arrTaxa = array_keys($arrFullCodingSeq);
	$nAlnLen = strlen($arrFullCodingSeq[$arrTaxa[0]]);

	if ($nAlnLen % 3 != 0) {
		die("Error: Gene $sGeneName not a multiplication of 3!\n");
	}

	$arrRetCol = array();
	$arrRetPos12 = array();
	$arrRet4fold = array();

	foreach($arrTaxa as $sTaxon) {
		if (!array_key_exists($sTaxon , $arrTaxaList) ) {
			$arrTaxaList[$sTaxon] = 0;
		}
		$arrTaxaList[$sTaxon]++;
		if (strlen($arrFullCodingSeq[$sTaxon]) != $nAlnLen ) {
			die("Error: Gene $sGeneName unequal length at taxon $sTaxon ! check if properly aligned!\n");
		}

		$arrRetCol[$sTaxon] = '';
		$arrRetPos12[$sTaxon] = '';
		$arrRet4fold[$sTaxon] = '';
		
	}

	//go over each alignment column.


	for($nPos=0;$nPos<$nAlnLen; $nPos+=3) {
		$arrCol = array();
		$arrPos12 = array();
		$arrLastPos = array();
		$bAll4fold = true;
		$sFirstTwoPos = "";
		$sAA = "";
		foreach($arrFullCodingSeq as $sTaxon => $sSeq) {
			$sCodon = substr($sSeq, $nPos, 3);
			list($sThisAA, $sTmp) = fnDNA2Peptide($sCodon, true, true);
			
			if ($sAA == '' && $sThisAA != '') {
				$sAA = $sThisAA;
			}

			$bIs4Fold = ($sAA == $sThisAA); //is synonymous?
			//$bIs4Fold = fnIs4FoldDegen($sCodon);
			/*
			if ($bIs4Fold === -1) {
				continue 2; //skip this position due to missing data
			}
			if ($bIs4Fold === true && $sFirstTwoPos=='') {
				$sFirstTwoPos = substr($sCodon, 2);
			}
			if ($bIs4Fold === true && $sFirstTwoPos!='') {
				$bIs4Fold = (substr($sCodon, 2) == $sFirstTwoPos);
			}
			*/ // do not look at first two codons now

			$bAll4fold = $bAll4fold && $bIs4Fold;
			$arrCol[$sTaxon] = $sCodon;
			$arrPos12[$sTaxon] = substr($sCodon, 0,2);
			$arrLastPos[$sTaxon] = $sCodon[2];
		}

		foreach($arrTaxa as $sTaxon) {
			$arrRetCol[$sTaxon] .= $arrCol[$sTaxon];
			$arrRetPos12[$sTaxon] .= $arrPos12[$sTaxon];
			if ($bAll4fold)  {
				$arrRet4fold[$sTaxon] .= $arrLastPos[$sTaxon];
			}
		}
		
	}

	$arrProccessedAln['cleancols'][$sGeneName] = $arrRetCol;
	$arrProccessedAln['pos12'][$sGeneName] = $arrRetPos12;
	$arrProccessedAln['4fold'][$sGeneName] = $arrRet4fold;
	
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

function fnWritePhylip($arrSeqs, $h) {
	global $arrTaxaList;
	$arrOutSeq = array();
	foreach($arrTaxaList as $sTaxon => $n) {
		$arrOutSeq[$sTaxon] = '';
	}

	$nTotalAlnLen = 0;
	foreach($arrSeqs as $sGeneName => $arrGeneSeqs) {
		list($sFirstKey) = array_keys($arrGeneSeqs);
		
		$nAlnLen = strlen($arrGeneSeqs[$sFirstKey]);

		foreach($arrTaxaList as $sTaxon => $n) {
			if (!array_key_exists($sTaxon , $arrGeneSeqs) ) {
				$arrOutSeq[$sTaxon] .= str_repeat('N', $nAlnLen);
				continue;
			}

			$arrOutSeq[$sTaxon] .=  $arrGeneSeqs[$sTaxon];
		}

		$nTotalAlnLen += $nAlnLen;
	}

	fwrite($h , " ".count($arrTaxaList)." $nTotalAlnLen\n");
	foreach($arrOutSeq as $sTaxon => $sSeq) {
		fwrite($h , "$sTaxon     $sSeq\n");
	}
}


?>
