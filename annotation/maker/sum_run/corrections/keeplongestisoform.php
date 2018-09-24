<?php
/*
This script keeps the best scoring isoform for each gene.
*/
require_once(dirname(__FILE__) . "/lib.php");
$sMakerGFF = "maker.refined.gff";
$sGenome = "query.fa";
$sOutGFF = "maker.refined.removeisoform.gff";
$nOverlapA=0.8;
$nOverlapB=0.8;

while(count($argv) > 0) {
    $arg = array_shift($argv);
    switch($arg) {
        case '-R':
            $sGenome = trim(array_shift($argv));
            break;
	case '-M':
	    $sMakerGFF  = trim(array_shift($argv));
	    break;
	case '-o':
	    $sOutGFF  = trim(array_shift($argv));
	    break;
	case '-A':
            $nOverlapA  = floatval( trim(array_shift($argv)));
	    break;
	case '-B':
            $nOverlapB  = floatval( trim(array_shift($argv)));
	    break;
	case '-S':
            $sRule  =  trim(array_shift($argv));
	    break;

    }
}


if (!file_exists($sMakerGFF)) {
	die("Cannot find $sMakerGFF\n");
}


$hOutGFF = fopen($sOutGFF , "w");

$oMakerGFF = new MappedCDSGFF();


$oMakerGFF->LoadGFF($sMakerGFF , $sGenome);


$arrGenes = $oMakerGFF->GetGenes('maker'); // get all maker genes

//collapse overlapping mRNAs under the same "gene"
$sTmpBed = "_tmpbed.bed";
$hTmpBed = fopen($sTmpBed , "w");
$nGeneRank = 0;
foreach($arrGenes as $sGeneID => &$oGene) {
	$nGeneRank++;
	foreach($oGene['mRNAs'] as $smRNAID => &$omRNA) {
		if ( !array_key_exists('score' , $omRNA['annot']) ) {
			if (strpos($smRNAID, "frag-merged-region") !==false ) { //keep the merged fragments!
				$nScore = 10000000;
			} else {
				$nScore = -10000000;
			}

		} else {
			$nScore = $omRNA['annot']['score'];
		}

		fwrite($hTmpBed, $omRNA['scf']."\t".$omRNA['start']."\t".$omRNA['end']."\t$sGeneID|$smRNAID|$nGeneRank\t$nScore\t".$omRNA['strand']."\n");

	}

}

//run bedtools
$sBedOut = "_tmp_intersect.bed.out";
//echo("overlapa: $nOverlapA overlapb: $nOverlapB\n"); die();
exec("bedtools intersect -a $sTmpBed  -b $sTmpBed  -wb -s -f $nOverlapA -F $nOverlapB > $sBedOut; ");

$hBedOut = fopen($sBedOut , "r");
$arrGeneIDReplace = array(); //key is the main gene id, value = array( gene id to delete )
$arrGeneIDReplaceInv = array(); //key is the gene id to delete , value is the main gene id it will be grouped into
while( false !== ($sLn = fgets($hBedOut) ) ) {
	$sLn = trim($sLn );
	if ($sLn == '') continue;
	$arrF = explode("\t" , $sLn);

	if ($arrF[3] == $arrF[9]) {
		continue; //overlap with itself, skip
	}


	//see if already grouped under 1 gene:
	$arrGeneName1 = explode("|", $arrF[3]);
	$arrGeneName2 = explode("|" , $arrF[9]);

	if ($arrGeneName1[0] == $arrGeneName2[0]) {
		continue;
	}

	if ($arrGeneName1[2] > $arrGeneName2[2]) { //gene 1 appears later in the gff file, don't look at it.
		continue;
	}

	//echo($sLn.PHP_EOL);

	if (array_key_exists($arrGeneName1[0] , $arrGeneIDReplace)) { //if 1 is main ID
		if (!array_key_exists($arrGeneName2[0] , $arrGeneIDReplace[$arrGeneName1[0]] ) ) {
			$arrGeneIDReplace[$arrGeneName1[0]][$arrGeneName2[0]] = true;
			$arrGeneIDReplaceInv[$arrGeneName2[0]] = $arrGeneName1[0];
		}
		continue;
	}

	if (array_key_exists($arrGeneName2[0] , $arrGeneIDReplace)) { //if 2 is main ID
		if (!array_key_exists($arrGeneName1[0] , $arrGeneIDReplace[$arrGeneName2[0]] ) ) {
			$arrGeneIDReplace[$arrGeneName2[0]][$arrGeneName1[0]] = true;
			$arrGeneIDReplaceInv[$arrGeneName1[0]] = $arrGeneName2[0];
		}
		continue;
	}

	if (array_key_exists($arrGeneName1[0] , $arrGeneIDReplaceInv)) { //if 1 is to be replaced
		if ($arrGeneIDReplaceInv[$arrGeneName1[0]] != $arrGeneName2[0]) {
			$arrGeneIDReplace[$arrGeneIDReplaceInv[$arrGeneName1[0]]][$arrGeneName2[0]] = true;
			$arrGeneIDReplaceInv[$arrGeneName2[0]] = $arrGeneIDReplaceInv[$arrGeneName1[0]];
		}
		continue;
	}

	if (array_key_exists($arrGeneName2[0] , $arrGeneIDReplaceInv)) { //if 2 is main ID
		if ($arrGeneIDReplaceInv[$arrGeneName2[0]] != $arrGeneName1[0]) {
			$arrGeneIDReplace[$arrGeneIDReplaceInv[$arrGeneName2[0]]][$arrGeneName1[0]] = true;
			$arrGeneIDReplaceInv[$arrGeneName1[0]] = $arrGeneIDReplaceInv[$arrGeneName2[0]];
		}
		continue;
	}

	$arrGeneIDReplace[$arrGeneName1[0]][$arrGeneName2[0]] = true;
	$arrGeneIDReplaceInv[$arrGeneName2[0]] = $arrGeneName1[0];


}



//die();
//now, move the mRNAs into the genes:
$arrGeneIDReplaceInv = array(); //reconstruct array in a clean way;

foreach ($arrGeneIDReplace as $sMainGeneID => &$arrReplaceGeneIDs) {
	foreach($arrReplaceGeneIDs as $sReplaceGeneID => $bDummy) {
		if ($sReplaceGeneID == $sMainGeneID) {
			continue;
		}

		$arrGeneIDReplaceInv[$sReplaceGeneID] = $sMainGeneID;

		foreach($arrGenes[$sReplaceGeneID]['mRNAs'] as $sRNAID => &$oRNA) {
			if (array_key_exists('Parent' , $oRNA['annot'])) {
				$oRNA['annot']['Parent'] = $sMainGeneID;
			}
			$arrGenes[$sMainGeneID]['mRNAs'][$sRNAID] = $oRNA;
		}

		unset($arrGenes[$sReplaceGeneID]); //delete this gene.
	}
}

//print_r($arrGeneIDReplace);

//print_r($arrGeneIDReplaceInv);
//die();

$arrGeneOperations = array();
$arrOperationCounts = array();
$nGeneCount = 0;
$arrKeepRNAID = array();

foreach($arrGenes as $sGeneID => &$oGene) {


	if (substr($sGeneID , 0,4) == 'trna') {
		continue;
	}
	echo("Checking $sGeneID ...\n");

	$arrIsoformScores = array();
	foreach($oGene['mRNAs'] as $smRNAID => &$omRNA) {
		if ( !array_key_exists('score' , $omRNA['annot']) ) {
			if (strpos($smRNAID, "frag-merged-region") !==false ) { //keep the merged fragments!
				$arrIsoformScores[10000000] = $smRNAID; 
			} else {
				$arrIsoformScores[-10000000] = $smRNAID; //if no scores, use a very negative value
			}
			continue;
		}

		$arrIsoformScores[$omRNA['annot']['score']] =  $smRNAID;
	}

	krsort($arrIsoformScores, SORT_NUMERIC); //Sort score in reverse order.

	$arrKeepRNAID[array_shift($arrIsoformScores) ] = $sGeneID; // put the first element. this has the best score.
}


$hMakerGFF = fopen($sMakerGFF , "r");

while( false !== ($sLn = fgets($hMakerGFF) )) {

	$sLn = trim($sLn);
	if ($sLn == "") continue;
	if ($sLn[0] == '#') {
		fwrite($hOutGFF , $sLn."\n");
		continue;
	}

	$arrF = explode("\t" , $sLn);

	if ($arrF[2] == "gene") { //always output the gene.
		preg_match("/ID=([^;]+)/", $arrF[8] , $arrM);
		if (count($arrM) == 2) {
			$sGeneID = trim($arrM[1]);
			if (array_key_exists($sGeneID , $arrGeneIDReplaceInv) ) {
				continue; //this gene id has been replaced.
			}
		}
		fwrite($hOutGFF , $sLn."\n");
		continue;
	} else if ($arrF[2] == "mRNA") {
		preg_match("/ID=([^;]+)/", $arrF[8] , $arrM);
		if (count($arrM) == 2) {
			$sRNAID = trim($arrM[1]);
			if (array_key_exists($sRNAID , $arrKeepRNAID) || substr($sRNAID , 0,4) == 'trna' ) {
				$sGeneID = $arrKeepRNAID[$sRNAID];
				$arrF[8] = preg_replace('/Parent=[^;]+/', "Parent=$sGeneID", $arrF[8]);
				fwrite($hOutGFF , implode("\t", $arrF)."\n");
				continue;
			}
		}

	} else {
		preg_match("/Parent=([^;]+)/", $arrF[8] , $arrM);
		if (count($arrM) == 2) {
			$sRNAID = trim($arrM[1]);
			if (array_key_exists($sRNAID , $arrKeepRNAID) || substr($sRNAID , 0,4) == 'trna' ) {
				fwrite($hOutGFF , $sLn."\n");
				continue;
			}
		}
	}


}


function fnArr2Annot($arr) {
	$s = "";
	
	foreach($arr as $sKey => $sVal) {
		if ($sVal === false) $sVal=0;
		$s .= "$sKey=$sVal;";
	}

	return substr($s, 0, strlen($s)-1);
}

function fnParseAnnoNote($sF,  &$omRNA) {
	global $sGenome;
	$arrAnn = array();

	$hF = fopen($sF, "r");
	$sLn = fgets($hF);
	if ($sLn === false ) return false;
	$arrF = explode("\t" , trim($sLn));
	if (count($arrF) < 2) return false;
	$arrAnn['mRNAID'] = $arrF[0];
	if ($arrAnn['mRNAID'] != $omRNA['annot']['ID'] ) {
		die("mRNA ID not equal!\n");
	}

	$arrAnn['operation'] = $sOperation = $arrF[1];
	$arrAnn['annotation'] = array_key_exists(2 , $arrF)? MappedCDSGFF::fnParseAnnotation($arrF[2]) : '';

	if ($sOperation == 'pass' || $sOperation == 'unknown'  ) { //no need to touch the original model, but pass in the annotations
		$arrAnn['NewRNA'] = $omRNA; //mRNA object unchanged
		if ($sOperation == 'pass') {
			$arrAnnParsed = $arrAnn['annotation'];
			$arrAnn['score'] = $arrAnnParsed['score'];
		} else {
			$arrAnn['score'] = 0;
		}

		return $arrAnn;
	}

	//otherwise, read in the fixed GFF as an object
	$sFixedGFF = dirname($sF)."/fixed.gff";
	if (!file_exists($sFixedGFF)) {
		die("$sFixedGFF not found!\n");
	}


	$oFixedGFF = new MappedCDSGFF();
	$oFixedGFF->LoadGFF($sFixedGFF , $sGenome);
	$arrNewGenes = $oFixedGFF->GetGenes('maker'); // there should be only a single mRNA annotation in here.


	
	foreach($arrNewGenes as $sGeneID => &$oNewGene ) {
		foreach($oNewGene['mRNAs'] as $sNewRNAID => &$oNewRNA) {
			$arrAnn['NewRNA'] = $oNewRNA;
			$arrAnn['score'] = $oNewRNA['annot']['score'];
			$arrAnn['annotation'] = $oNewRNA['annot'];
			return $arrAnn;
			break 2; // only need to look at 1 gene because this is the only one
		}		
	}

}

function fnWriteNewModel($arrNewModel, $oScore, $hFixedGFF, $sPredictor='maker', $bWriteGene=true) { // do not add "gene" line
	if ($arrNewModel['strand'] == '+') {
		ksort($arrNewModel['CDSs']);
	} else {
		krsort($arrNewModel['CDSs']);
	}
/*
NORv12scf258	exonerate	gene	457384	528041	9770	-	.	ID=gene.NORv12scf258:hit:16476:3.10.0.0
NORv12scf258	exonerate	mRNA	457384	528041	9770	-	.	ID=NORv12scf258:hit:16476:3.10.0.0;Name=gi|1007780707|ref|XP_015827581.1|;Parent=gene.NORv12scf258:hit:16476:3.10.0.0
NORv12scf258	exonerate	exon	527937	528041	9770	-	.	ID=NORv12scf258:hsp:20104:3.10.0.0;Parent=NORv12scf258:hit:16476:3.10.0.0;Target=gi|1007780707|ref|XP_015827581.1| 1 35;Length=2344;Gap=M35
NORv12scf258	exonerate	CDS	527937	528041	9770	-	.	ID=NORv12scf258:hsp:20104:3.10.0.0;Parent=NORv12scf258:hit:16476:3.10.0.0;Target=gi|1007780707|ref|XP_015827581.1| 1 35;Length=2344;Gap=M35


Array
(
    [scf] => NORv12scf258
    [start] => 457384
    [end] => 528041
    [strand] => -
    [annot] => Array
        (
            [ID] => SplitChimeric_432.NORv12scf258:hit:22205:4.5.0.0
            [Name] => gi|1007780701|ref|XP_015827578.1|
            [Parent] => SplitChimeric_432.gene.NORv12scf258:hit:22205:4.5.0.0
        )

    [CDSs] => Array
        (
            [457384] => 457765
            [457841] => 458009

*/
	$arrF = array($arrNewModel['scf'] , $sPredictor , 'gene', $arrNewModel['start'] , $arrNewModel['end'],  '.', $arrNewModel['strand'] , '.', "ID=".$arrNewModel['annot']['Parent']);

	if ($bWriteGene) {fwrite($hFixedGFF , implode("\t" , $arrF)."\n");}

	$arrF[2] = 'mRNA';
	$arrF[5] = array_key_exists('score', $oScore)? $oScore['score'] : (array_key_exists('score', $arrNewModel['annot'])? $arrNewModel['annot']['score'] : '.' );
	$arrF[8] = fnArr2Annot($arrNewModel['annot']);
	$arrF[8] .= array_key_exists('score', $oScore)? ";". fnArr2Annot($oScore) : "";
	fwrite($hFixedGFF , implode("\t" , $arrF)."\n");

	$arrF[5] = '.';
	$nCurrLen = 0;
	$nCDScount = 0;
	foreach($arrNewModel['CDSs'] as $nStart => $nEnd) {
		$nCDScount++;
		$nCDSLen = abs($nEnd - $nStart) + 1;
		$nPrevLeftBase = $nCurrLen % 3;
		$nPhase = 3 - $nPrevLeftBase;
		if ($nPhase == 3) $nPhase = 0;
		$nCurrLen += $nCDSLen;
		$arrF[2] = 'exon';
		$arrF[3] = $nStart;
		$arrF[4] = $nEnd;
		$arrF[7] = '.';
		$arrF[8] = "ID=".$arrNewModel['annot']['ID'].".exon$nCDScount;Parent=".$arrNewModel['annot']['ID'];
		fwrite($hFixedGFF , implode("\t" , $arrF)."\n");
		$arrF[2] = 'CDS';
		$arrF[7] = $nPhase;
		$arrF[8] = "ID=".$arrNewModel['annot']['ID'].".cds$nCDScount;Parent=".$arrNewModel['annot']['ID'];
		fwrite($hFixedGFF , implode("\t" , $arrF)."\n");
	}

}

?>
