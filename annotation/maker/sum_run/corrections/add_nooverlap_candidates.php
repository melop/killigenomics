<?php
//this script appends candidate gene models that are not overlapping with any final maker model.
// I found that sometimes maker fails to promote certain candidate models (especially genewise) 
// because sometimes exonerate does not generate a protein alignment to support them
// but genewise is more sensitive
// These models will not be accessed by protein alignments , but selected based on length (longer ones are kept in case they self overlap
require_once(dirname(__FILE__) . "/lib.php");

$sCandidateModels = "split__candidate.gff";
$sMakerModels = "maker.refined.removeisoform.gff";
$sToAppend = "append.candidates.gff";
$sGenome = "query.fa";
$sCombined = "maker.refined.removeisoform.missedaddedback.gff";


exec('grep -P "\tmRNA\t" '. $sCandidateModels .' > candidate.mRNA.gff');
exec('grep -P "\tmRNA\t" ' . $sMakerModels .' | grep -v "trna" > currentmakermodel.mRNA.gff');
exec('bedtools intersect -a candidate.mRNA.gff -b currentmakermodel.mRNA.gff -v -s > nonoverlap_candidate.gff ');
exec('bedtools intersect -a nonoverlap_candidate.gff -b nonoverlap_candidate.gff -r -f 0.5 -s -wb > candidate_selfoverlap.gff ');

//now read in the self overlap candidate file:
$hF = fopen('candidate_selfoverlap.gff' , 'r');

$arrGeneIDKeep = array(); //key is the main gene id, value = array( gene id to delete )
$arrGeneIDDelete = array(); //key is the gene id to delete , value is the main gene id it will be grouped into

while( false !== ( $sLn = fgets($hF ) )) {
	$sLn = trim($sLn);
	if ($sLn == '') continue;
	if ($sLn[0] == '#') continue;

	$arrF = explode("\t", $sLn);
	if (count($arrF) != 18) {
		continue;
	}

	$arrAnnot1 = MappedCDSGFF::fnParseAnnotation($arrF[8]);
	$arrAnnot2 = MappedCDSGFF::fnParseAnnotation($arrF[17]);


	$arrQI1 = explode("|", $arrAnnot1['_QI']);
	$arrQI2 = explode("|", $arrAnnot2['_QI']);

	$nLen1 = array_pop($arrQI1);
	$nLen2 = array_pop($arrQI2);

	if ($arrAnnot1['ID'] == $arrAnnot2['ID']) {
		if (!array_key_exists($arrAnnot1['ID'] , $arrGeneIDDelete)) { // if gene not in delete list
			$arrGeneIDKeep[$arrAnnot1['ID']] = $nLen1; //keep this gene
		}

		continue;
	}


	$sName1 = $arrAnnot1['Name'];
	$sName2 = $arrAnnot2['Name'];

	$nTypeScore1 = fnTypeScore($sName1);
	$nTypeScore2 = fnTypeScore($sName2);
	if ($nTypeScore1 == 0 && $nTypeScore2 == 0) {
		continue; //both are abinitio models, delete.
	}

	$bModel1Better = ($nTypeScore1 == $nTypeScore2)? ($nLen1>$nLen2) : ($nTypeScore1 > $nTypeScore2);

	if ($bModel1Better) { // keep gene 1, discard gene 2
		//keep 1, delete 2
		if (!array_key_exists($arrAnnot1['ID'] , $arrGeneIDDelete)) { // if gene not in delete list
			$arrGeneIDKeep[$arrAnnot1['ID']] = $nLen1; //keep this gene
			$arrGeneIDDelete[$arrAnnot2['ID']] = $nLen2;
			if (array_key_exists($arrAnnot2['ID'] , $arrGeneIDKeep)) {
				unset($arrGeneIDKeep[$arrAnnot2['ID']]);
			}
		}

	} else {
		//keep 1, delete 2
		if (!array_key_exists($arrAnnot2['ID'] , $arrGeneIDDelete)) {
			$arrGeneIDKeep[$arrAnnot2['ID']] = $nLen2; //keep this gene
			$arrGeneIDDelete[$arrAnnot1['ID']] = $nLen1;
			if (array_key_exists($arrAnnot1['ID'] , $arrGeneIDKeep)) {
				unset($arrGeneIDKeep[$arrAnnot1['ID']]);
			}
		}


	}


}



echo(count($arrGeneIDKeep)." non-overlapping candidate models will be added back...\n");

//print_r($arrGeneIDKeep);

$hCandidateGFF = fopen($sCandidateModels , "r");
$hOut = fopen($sToAppend , 'w');


while( false !== ($sLn = fgets($hCandidateGFF) )) {

		if ($sLn == "") continue;
		if ($sLn[0] =='#') continue;

		$arrF = explode("\t", trim( $sLn) );
		$arrF[1] = "maker";

		//print_r($arrF);
		if (count($arrF) < 5) continue;

		if ($arrF[2] == "mRNA") {
			preg_match("/ID=([^;]+);/" , $arrF[8]  , $arrM);
			if (count($arrM) != 2) continue;
			$sID = trim( $arrM[1] );
			$sGeneID = "gene.$sID";
		
			if (array_key_exists($sID  , $arrGeneIDKeep)) {

				$arrF2 = $arrF;
				$arrF2[2] = "gene";
				$arrF2[8] = "ID=$sGeneID";
				$arrF[8] .= ";Parent=$sGeneID";
				fwrite($hOut , implode("\t", $arrF2)."\n" );
				fwrite($hOut , implode("\t", $arrF)."\n" );

			}
			continue;
		}

		if ($arrF[2] == "exon") {
			preg_match("/Parent=([^;]+);/" , $arrF[8]  , $arrM);
			if (count($arrM) != 2) continue;
			$sID = trim( $arrM[1] );
		
			if (array_key_exists($sID  , $arrGeneIDKeep)) {
				fwrite($hOut , implode("\t", $arrF)."\n" );
				$arrF[2] = "CDS";
				fwrite($hOut , implode("\t", $arrF)."\n" );

			}
			continue;
		}

	


}

//add phase information to CDS.

$oFixedGFF = new MappedCDSGFF();
$oFixedGFF->LoadGFF($sToAppend  , $sGenome);
$arrNewGenes = $oFixedGFF->GetGenes('maker'); // load all genes.

fclose($hOut);
$hOut = fopen($sToAppend , 'w');

foreach($arrNewGenes as $sGeneID => &$oNewGene ) {
	foreach($oNewGene['mRNAs'] as $sNewRNAID => &$oNewRNA) {
		fnWriteNewModel($oNewRNA, array(), $hOut, 'maker', true);
	}
}

exec("cat $sMakerModels  $sToAppend > $sCombined ");


function fnTypeScore($sName) {
	$sName = strtolower($sName);
	if (strpos($sName, 'abinit') !== false ) {
		return 0;
	}
	if (strpos($sName, 'evm') !== false ) {
		return 1;
	}
	if (strpos($sName, 'genewise') !== false ) {
		return 2;
	}

	return 0;

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


function fnArr2Annot($arr) {
	$s = "";
	
	foreach($arr as $sKey => $sVal) {
		if ($sVal === false) $sVal=0;
		$s .= "$sKey=$sVal;";
	}

	return substr($s, 0, strlen($s)-1);
}

?>
