<?php
require_once(dirname(__FILE__) . "/lib.php");
$sMakerGFF = "chimericfixed.gff";
$sWorkDir = "./refinement_improved/";
$sGenome = "query.fa";
$sOutGFF = "maker.refined.gff";

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
	    $sWorkDir  = trim(array_shift($argv));
	    break;

    }
}

if (!file_exists($sWorkDir)) {
	die("Cannot find $sWorkDir\n");
}

if (!file_exists($sMakerGFF)) {
	die("Cannot find $sWorkDir\n");
}


$hOutGFF = fopen($sOutGFF , "w");

$oMakerGFF = new MappedCDSGFF();


$oMakerGFF->LoadGFF($sMakerGFF , $sGenome);


$arrGenes = $oMakerGFF->GetGenes('maker'); // get all maker genes

$arrGeneOperations = array();
$arrOperationCounts = array();
$nGeneCount = 0;
foreach($arrGenes as $sGeneID => &$oGene) {

	//if ($nGeneCount++ >= 10) break;
	if (substr($sGeneID , 0,4) == 'trna') {
		continue;
	}
	echo("Checking $sGeneID ...\n");

	$sGeneDIR = "$sWorkDir/".$oGene['scf']."/".md5($sGeneID);
	$bRevised = false;


	$arrOperation = array(); // this saves the operation for each mRNA isoform, key is the score of the model.
	$bPassOrReplacementFound = false;
	$bNewGeneStart = 0;
	$bNewGeneEnd = 0;

	foreach($oGene['mRNAs'] as $smRNAID => &$omRNA) {

		$sRNADIR = $sGeneDIR."/".md5($smRNAID);
		$sFixedGFF = "$sRNADIR/fixed.gff";
		$sAnntNote = "$sRNADIR/annotation.notes";

		if (!file_exists($sAnntNote) ) {
			echo("Warning, refinement output not found for $smRNAID, in $sAnntNote \n");
			continue;
		}
		echo("$sAnntNote\n");
		$oAnnot = fnParseAnnoNote($sAnntNote, $omRNA);

		if ($oAnnot === false) {
			die("Annotation output incomplete, has the program finished?\n");
		}

		if (!array_key_exists($oAnnot['operation'], $arrOperationCounts)) {
			$arrOperationCounts[$oAnnot['operation']] = 0;
		}
		$arrOperationCounts[$oAnnot['operation']] += 1;
		$oAnnot['NewRNA']['annot'] = $oAnnot['annotation'];
		$oAnnot['NewRNA']['annot']['FixResult'] = $oAnnot['operation'];
		$oAnnot['NewRNA']['annot']['ID'] = $smRNAID; //use the original mRNAID
		$oAnnot['NewRNA']['annot']['Parent'] = $sGeneID; //use the original gene id

		if ($bNewGeneStart == 0) { $bNewGeneStart = $oAnnot['NewRNA']['start']; }
		if ($bNewGeneEnd == 0) { $bNewGeneEnd = $oAnnot['NewRNA']['end']; }
		if ($bNewGeneStart > $oAnnot['NewRNA']['start']) { $bNewGeneStart = $oAnnot['NewRNA']['start']; }
		if ($bNewGeneEnd < $oAnnot['NewRNA']['end']) { $bNewGeneEnd = $oAnnot['NewRNA']['end']; }

		$arrOperation[$smRNAID] = $oAnnot; //if they refer to the same model, then they will have the same score.

	}

	ksort($arrOperation);

	if (count($arrOperation) ==0 ) {
		die("Something is wrong with the refinement output. count(arrScores) = 0! \n");
	}

	$arrOperation['_genebound'] = array($bNewGeneStart , $bNewGeneEnd);

	//
	$arrGeneOperations[$sGeneID] = $arrOperation;

	//print_r($arrGeneOperations);

	//die();

}

echo(implode("\t" , array_keys($arrOperationCounts)) . "\n" );
echo(implode("\t" , array_values($arrOperationCounts)) . "\n" );
//die();

//print_r($arrGeneOperations);
//die();

echo("Correcting...\n");
$hMakerGFF = fopen($sMakerGFF , "r");
$arrDeletemRNA = array();
$arrKeepOperations = array('pass' => true , 'unknown' => true, 'notpass_cannotfix' => true); //these are operation types that don't need to change the gene, only need to update the annotation string
$arrSwapOperations = array('replaced_by_exonerate' => true , 'replaced_by_cand' => true); //these are operations types that require deleting the original mRNA.


while( false !== ($sLn = fgets($hMakerGFF) )) {

	$sLn = trim($sLn);
	if ($sLn == "") continue;
	if ($sLn[0] == '#') {
		fwrite($hOutGFF , $sLn."\n");
		continue;
	}

	$arrF = explode("\t" , $sLn);
	$arrOut = $arrF;

	if ($arrF[2] == "gene") {
		preg_match("/ID=([^;]+)/", $arrF[8] , $arrM);
		if (count($arrM) == 2) {
			$sID = trim($arrM[1]);
			if (strpos($sID, 'trna') === false) { //if this is not trnagene 
				echo("Gene $sID\n");
				if (array_key_exists($sID, $arrGeneOperations)) { //check if need to perform operations on gene
					$oOper = $arrGeneOperations[$sID];
					$nNewGeneStart = $oOper['_genebound'][0];
					$nNewGeneEnd = $oOper['_genebound'][1]; //check against the current gene boundary
					$arrOut[3] = min($arrOut[3] , $nNewGeneStart);
	 				$arrOut[4] = max($arrOut[4] , $nNewGeneEnd);
				}
			}
		}
	} else if ($arrF[2] == "mRNA") {
		preg_match("/Parent=([^;]+)/", $arrF[8] , $arrM);
		if (count($arrM) == 2) {
			$sParentID = trim($arrM[1]);
			preg_match("/ID=([^;]+)/", $arrF[8] , $arrM);
			$sID = trim($arrM[1]);
			if (strpos($sParentID, 'trna') === false) { //if this is not trnagene 
				if (array_key_exists($sParentID , $arrGeneOperations)) {
					$oOper = $arrGeneOperations[$sParentID];
					if (array_key_exists($sID , $oOper) ) {
						echo("mRNA $sID\n");
						$mRNAOper = $oOper[$sID];
						if (array_key_exists($mRNAOper['operation'] , $arrKeepOperations) ) { //The operation is to keep the original mRNA
							//update some annotations
							$mRNAOper['NewRNA']['annot']['ID'] = $sID;
							$mRNAOper['NewRNA']['annot']['Parent'] = $sParentID;
							$arrOut[5] = array_key_exists('score', $mRNAOper['NewRNA']['annot'])? $mRNAOper['NewRNA']['annot']['score'] : '.';
							$arrOut[8] = fnArr2Annot($mRNAOper['NewRNA']['annot']);
							echo("Keep mRNA $sID ...\n");
						}
						if (array_key_exists($mRNAOper['operation'] , $arrSwapOperations)) { //If the operation is to replace
							$mRNAOper['NewRNA']['annot']['ID'] = $sID;
							$mRNAOper['NewRNA']['annot']['Parent'] = $sParentID;
							fnWriteNewModel($mRNAOper['NewRNA'], array(), $hOutGFF, 'maker', false) ; 
							$arrDeletemRNA[$sID] = true;
							echo("replaced mRNA $sID with better model...\n");
							continue; //mark this mRNA for deletion.
						}
					}
					else { //this mRNA was discarded.
						$arrDeletemRNA[$sID] = true;
						continue; //mark this mRNA for deletion.
					}
				} //if gene id not found, then output 
			}
		}

	} else {
		preg_match("/Parent=([^;]+)/", $arrF[8] , $arrM);
		if (count($arrM) == 2) {
			$sID = trim($arrM[1]);
			if (array_key_exists($sID , $arrDeletemRNA)) {
				//echo("remove part $sID\n");
				continue;
			}
		}
	}


	fwrite($hOutGFF , implode("\t", $arrOut)."\n");

}

echo("\nFINISHED!\n");

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
