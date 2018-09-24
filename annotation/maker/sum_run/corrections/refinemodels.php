<?php
require_once(dirname(__FILE__) . "/lib.php");
$sMakerGFF = "chimericfixed.gff";
$sGenome = "query.fa";
$sProteinGFF = "split__protein2genome.gff";
$sCandidateGFF = "split__candidate.gff";
$sProteinFasta = "all.protein.evidence.fa";
$sWorkDir = "./refinement_improved/";
$nMinBlastPIdentity = 50; //filter HSP if lower than this.
$nMaxBlastPIdentityDiff = 20; //If identity of an HSP block is 20% smaller than other blocks in the same alignment, then exclude this hsp
$nTreatAsBigGap = 2; // if larger than this, treat as "big gap".
$nTreatAsHugeGap = 10; // if larger than this, treat as "huge gap".
$sDoOnly = ""; //only do this md5 mRNA name
$sForceRedo = false;
$nHardLimitFromNeighbors = -20;

$nThisPart = 0;
$nTotalParts = 1;
//$sTest = "000000000000000111111111111111111110000001111001000000000000000002222222222222222222222200000000001110000000";
//echo($sTest."\n".fnSmoothCov($sTest)."\n");
//die();

while(count($argv) > 0) {
    $arg = array_shift($argv);
    switch($arg) {
        case '-R':
            $sGenome = trim(array_shift($argv));
            break;
	case '-M':
	    $sMakerGFF  = trim(array_shift($argv));
	    break;
	case '-P':
	    $sProteinGFF  = trim(array_shift($argv));
	    break;
	case '-C':
	    $sCandidateGFF  = trim(array_shift($argv));
	    break;
	case '-D':
	    $sDoOnly  = trim(array_shift($argv));
	    break;
	case '-p':
	    $sProteinFasta  = trim(array_shift($argv));
	    break;
	case '-o':
	    $sWorkDir  = trim(array_shift($argv));
	    break;
	case '-N':
	    $nTotalParts  = intval(trim(array_shift($argv)));
	    break;
	case '-f':
	    $nThisPart  = intval(trim(array_shift($argv)));
	    break;
	case '-m':
	    $nMinBlastPIdentity  = floatval(trim(array_shift($argv)));
	    break;
	case '-d':
	    $nMaxBlastPIdentityDiff  = floatval(trim(array_shift($argv)));
	    break;
	case '-g':
	    $nTreatAsBigGap  = intval(trim(array_shift($argv)));
	    break;
	case '-G':
	    $nTreatAsHugeGap  = intval(trim(array_shift($argv)));
	    break;
	case '-F':
	    $sForceRedo = true;

    }
}

exec("mkdir -p $sWorkDir");

$oProteinFasta = new FastaHack();
$oProteinFasta->SetDb($sProteinFasta);

$oMakerGFF = new MappedCDSGFF();
$oMakerGFF->LoadGFF($sMakerGFF , $sGenome);




$arrGenes = $oMakerGFF->GetGenes('maker'); // get all maker genes

$nGeneCount = -1;

foreach($arrGenes as $sGeneID => &$oGene) {

	$nGeneCount++;
	if ( $nGeneCount % $nTotalParts != $nThisPart) continue;

	echo("Looking at $sGeneID ...\n");

	$sGeneDIR = "$sWorkDir/".$oGene['scf']."/".md5($sGeneID);
	$oNextGeneCoord = $oMakerGFF->GetNeighborGeneOnScf($oGene['scf'] , $oGene['start'], $oGene['strand'], $nOffSet = 1);
	$oPrevGeneCoord = $oMakerGFF->GetNeighborGeneOnScf($oGene['scf'] , $oGene['start'], $oGene['strand'], $nOffSet = -1);

	$nLeftHardLimit = ($oPrevGeneCoord === false)? 0 : ($oPrevGeneCoord['end']+$nHardLimitFromNeighbors);
	$nRightHardLimit = ($oNextGeneCoord === false)? PHP_INT_MAX : ($oNextGeneCoord['start']-$nHardLimitFromNeighbors);

	foreach($oGene['mRNAs'] as $smRNAID => &$omRNA) {

		$smRNAIDMD5= md5($smRNAID);

		if ($sDoOnly != '') {
			if ($sDoOnly != $smRNAIDMD5) {
				continue;
			}
		}

		$sRNADIR = $sGeneDIR."/".$smRNAIDMD5;
		$sScoreNotes = "$sRNADIR/annotation.notes";
		$sFixedGFF = "$sRNADIR/fixed.gff";

		if (!$sForceRedo) {
			if ( file_exists($sScoreNotes) ) { //if there are already results
				if (filesize($sScoreNotes) > 0) {
					echo("computed, skip.\n");
					continue;
				}
			}
		}

		exec("mkdir -p $sRNADIR");

		$hScoreNotes = fopen($sScoreNotes , 'w');

		echo(" -- mRNA  $smRNAID in dir $sRNADIR ...\n");
		$arrRNARange = array( array($omRNA['scf'] , $omRNA['start'], $omRNA['end'] , $omRNA['strand']) );

		$sLocalProteinGFF = fnGetOverlapGFF($sRNADIR , $sProteinGFF, $arrRNARange, 0.30, $nPercOverlapB=0.10);
		$oLocalProteinGFF = new MappedCDSGFF();
		$oLocalProteinGFF->LoadGFF($sLocalProteinGFF , $sGenome);

		if ($oLocalProteinGFF === false) {
			echo("Cannot read in evidence protein GFF $sLocalProteinGFF for $smRNAID. Leave as is.\n");
			fwrite( $hScoreNotes , "$smRNAID\tunknown\tquality=unkown\n");
			continue;

		}


		$sLocalCandidateModelGFF = fnGetOverlapGFF($sRNADIR , $sCandidateGFF, $arrRNARange, 0.30, 0.10);
		$oLocalCandidateModelGFF = new MappedCDSGFF();
		$oLocalCandidateModelGFF->LoadGFF($sLocalCandidateModelGFF , $sGenome);


		$arrExtracted = $oMakerGFF->ExtractmRNASequence('maker', $smRNAID);

		// now score this main model:
		$sMainModelProt = $arrExtracted['AA'];
		$oMainModelScore = fnGetBestScore($sMainModelProt , $oLocalProteinGFF , $sRNADIR );


		if ($oMainModelScore['supports'] == 0) {
			echo("Model $smRNAID has no protein match, cannot evaluate. Leave as is.\n");
			fwrite( $hScoreNotes , "$smRNAID\tunknown\tquality=unkown\n");
			continue;
		}

		if ($oMainModelScore['pass'] ) {
			echo("Model $smRNAID passed examination\n");
			fwrite( $hScoreNotes , "$smRNAID\tpass\t".fnArr2Annot($oMainModelScore['bestscore'])."\n");
			continue;
		}

		echo("Model $smRNAID did not pass examination, best score = ".  $oMainModelScore['notpassbestscore']['score']."\n ". fnArr2Annot($oMainModelScore['notpassbestscore']) ."\n try to find a better model...\n");
		if ($oMainModelScore['notpassbestscore']['hicov'] > 0) {
			echo("Warning: Repeated protein detected! Could be misassembly!\n");
			echo(">pred\n".$sMainModelProt."\n");

		}

		if ($oLocalCandidateModelGFF !== false && ($oMainModelScore['notpassbestscore']['hicov'] == 0) ) { //if the hicov is not 0, then only look at exonerate models!
			//try to look into the candidate models to see if can find a good one:
			$arrGoodCandModels = array();
			foreach($oLocalCandidateModelGFF->GenePredictors() as $sPredictor) {
				$arrCandidateModels = $oLocalCandidateModelGFF->GetGenes($sPredictor);
				foreach( $arrCandidateModels  as $sCandidateModelID => &$arrCandidateModelRNA ) {
					foreach($arrCandidateModelRNA['mRNAs'] as $sCandidateRNAID => &$oCandidateModel) {
						if ($oCandidateModel['start'] < $nLeftHardLimit  || $oCandidateModel['end'] > $nRightHardLimit ) {
							echo("\t\tCandidate model $sCandidateRNAID overlaps with neighboring genes, ignore.\n");
							continue;
						}

						$arrExtracted = $oLocalCandidateModelGFF->ExtractmRNASequence($sPredictor , $sCandidateRNAID);
						$oCandModelScore = fnGetBestScore($arrExtracted['AA'], $oLocalProteinGFF , $sRNADIR , max($oMainModelScore['notpassbestscore']['completeness'] , 0.95) );
						if ($oCandModelScore['pass']) {
							$arrGoodCandModels[$oCandModelScore['bestscore']['score']] = array($oCandidateModel, $oCandModelScore['bestscore'] );
						} 
					}
				}
			}
		

		
			$oBestCandModel = false;

			if (count($arrGoodCandModels) > 0 ) {
				//a passing candidate model found
				krsort($arrGoodCandModels , SORT_NUMERIC); //sort by score
				$arrScores = array_keys($arrGoodCandModels);
				$oBestCandModel = $arrGoodCandModels[$arrScores[0]];
				if ($arrScores[0] > $oMainModelScore['notpassbestscore']['score']) { // if the score is better than the non-passing score of the main model
					fnWriteNewModel( $oBestCandModel[0] , $oBestCandModel[1],  fopen($sFixedGFF , 'w') );
					fwrite( $hScoreNotes , "$smRNAID\treplaced_by_cand\t".$oBestCandModel[0]['annot']['ID']."\n");
					echo ("Candidate model used in place of the original maker model: ".$oBestCandModel[0]['annot']['ID'] .  "\n" . fnArr2Annot($oBestCandModel[1]) . "\n" );
					continue;
				}
			}

			echo("No candidate model pass examination try to find in exonerate protein models...\n");

		} else {
			echo("Candidate GFF contains some errors, skip\n");
		}

		// otherwise, put the best protein alignment as the final model
		$arrGoodExonerateModels = array();
		foreach($oLocalProteinGFF->GenePredictors() as $sPredictor) {
			$arrLocalProteinGenes = $oLocalProteinGFF->GetGenes($sPredictor);
			foreach( $arrLocalProteinGenes  as $sProteinGFFID => &$arrProteinModels ) {
				foreach($arrProteinModels['mRNAs'] as $sProtID => &$oProteinModel) {
					if ($oProteinModel['start'] < $nLeftHardLimit  || $oProteinModel['end'] > $nRightHardLimit ) {
						echo("\t\Exonerate model $sProtID overlaps with neighboring genes, ignore.\n");
						continue;
					}

					$arrExtracted = $oLocalProteinGFF->ExtractmRNASequence($sPredictor , $sProtID);
					$sProteinQueryName = $oProteinModel['annot']['Name'];
					$sProteinSeq = $oProteinFasta->GetContig($sProteinQueryName);
					//echo("$sProteinSeq\n");
					$oScore = fnScoreExonerateModel($arrExtracted['AA'] , $sProteinSeq , $sRNADIR ); //get a score for this underlying protein
					if ( (!$oScore['internalstop']) && ($oScore['score'] > 0) && ($oScore['hicov'] == 0)) {
						$arrGoodExonerateModels[$oScore['score']] = array($oProteinModel , $oScore , $sProteinSeq);
					} 
				}
			}
		}

		if (count($arrGoodExonerateModels) > 0 ) {
			//a passing candidate model found
			krsort($arrGoodExonerateModels , SORT_NUMERIC); //sort by score
			$arrScores = array_keys($arrGoodExonerateModels);
			if ( ($arrScores[0] > $oMainModelScore['notpassbestscore']['score'])  || // if the score is better than the non-passing score of the main model
                               ($oMainModelScore['notpassbestscore']['hicov'] > 0  && $arrGoodExonerateModels[$arrScores[0]][1]['hicov'] == 0) ) // or main model is repeated but this model is not
			{ 
				fnWriteNewModel( $arrGoodExonerateModels[$arrScores[0]][0] , $arrGoodExonerateModels[$arrScores[0]][1],  fopen($sFixedGFF , 'w') );
				fwrite( $hScoreNotes , "$smRNAID\treplaced_by_exonerate\t".$arrGoodExonerateModels[$arrScores[0]][0]['annot']['ID']."\n");
				echo ("Exonerate model used in place of the original maker model: ".$arrGoodExonerateModels[$arrScores[0]][0]['annot']['ID'] . "\n" . fnArr2Annot($arrGoodExonerateModels[$arrScores[0]][1]) . "\n");
				echo(">pred\n". ($oLocalProteinGFF->ExtractmRNASequenceByModel($arrGoodExonerateModels[$arrScores[0]][0]))['AA'] ."\n>evidence\n".$arrGoodExonerateModels[$arrScores[0]][2]."\n");
				continue;
			}
		}

		//cannot really find a better option...
		fnWriteNewModel( $omRNA  , $oMainModelScore['notpassbestscore'] ,  fopen($sFixedGFF , 'w') );
		fwrite( $hScoreNotes , "$smRNAID\tnotpass_cannotfix\t\n");
		echo ("No better model is found to replace the Maker model: ".$smRNAID . "\n");


		
	}
}

function fnGetBestScore(&$sMainModelProt, &$oLocalProteinGFF, $sRNADIR, $nMinCompleteness=0.95 , $MaxHugeInsertionCount=0, $MaxHugeDeletionCount=0, $MaxHiCov = 0 , $MaxBigInsertionLen=5, $MaxBigDeletionLen=5) {
	global $oProteinFasta;
		$nSupports = 0;
		$oBestScore = false;
		$oNotPassBestScore =  false; //saves best scores for the models that do not pass cutoffs
		$bPass = false;

		foreach($oLocalProteinGFF->GenePredictors() as $sPredictor) {
			$arrLocalProteinGenes = $oLocalProteinGFF->GetGenes($sPredictor);
			foreach( $arrLocalProteinGenes  as $sProteinGFFID => &$arrProteinModels ) {
				foreach($arrProteinModels['mRNAs'] as &$oProteinModel) {
					$sProteinQueryName = $oProteinModel['annot']['Name'];
					$sProteinID = fnSafeID( $oProteinModel['annot']['ID'] );//needed for fastahack
					//echo("$sProteinQueryName\n");
					$sProteinSeq = $oProteinFasta->GetContig($sProteinQueryName);
					//echo("$sProteinSeq\n");
					$nSupports++;
					$oScore = fnScoreExonerateModel($sMainModelProt , $sProteinSeq , $sRNADIR ); //get a score for this underlying protein

/*
	score
	5primemiss
	5primemissperc
	3primemiss
	3primemissperc
	5primeuncertain
	3primeuncertain
	completeness
	hicov
	smallinsertionlen
	smallinsertioncount
	biginsertionlen
	biginsertioncount
	hugeinsertionlen
	hugeinsertioncount

	smalldeletionlen
	smalldeletioncount
	bigdeletionlen
	bigdeletioncount
	hugedeletionlen
	hugedeletioncount
*/
					if ($oScore['score'] == 0 || $oScore['internalstop'] ||  
							 $oScore['completeness'] < $nMinCompleteness  || $oScore['hugeinsertioncount'] > $MaxHugeInsertionCount || 
							 $oScore['hugedeletioncount'] > $MaxHugeDeletionCount || $oScore['hicov'] > $MaxHiCov ||
                                                         $oScore['biginsertionlen'] > $MaxBigInsertionLen || $oScore['bigdeletionlen'] > $MaxBigDeletionLen ) {
						//check for underlying models
						if ($oNotPassBestScore === false) {
							$oNotPassBestScore = $oScore;
						} else if ($oNotPassBestScore['score'] < $oScore['score']) {
							$oNotPassBestScore = $oScore;
						}
						continue; //don't pass
					} else {
						$bPass = true;
						if ($oBestScore === false) {
							$oBestScore = $oScore;
						} else if ($oBestScore['score'] < $oScore['score']) {
							$oBestScore = $oScore;
						}
					}
				}
			}
		}

	return array('pass'=>$bPass , 'supports' => $nSupports , 'bestscore' => $oBestScore , 'notpassbestscore' => $oNotPassBestScore);
}

function fnGetOverlapGFF($sTmpDir , $sGFF, $arrRanges, $nPercOverlapA=0.01, $nPercOverlapB=0.01) {

	$sFilter = "$sTmpDir/overlapped.genes.txt";
	$sOut = "$sTmpDir/proteinsubset.gff";

	$sBed = "$sTmpDir/ranges.bed";
	$hBed = fopen($sBed, "w");

	exec("grep -P \"\tmRNA\t\" $sGFF > $sTmpDir/mRNA.gff");


	foreach($arrRanges as $oRange) {
		fwrite($hBed , $oRange[0]."\t".($oRange[1]-1)."\t$oRange[2]\tA\t1000\t$oRange[3]\n"); 
	}

	exec("bedtools intersect -a $sBed -b $sTmpDir/mRNA.gff -wb -s -f $nPercOverlapA -F $nPercOverlapB > $sFilter; rm $sTmpDir/mRNA.gff");



	$hFilter = fopen($sFilter, "r");
	$hOut = fopen($sOut , "w");

	$arrCanGenes = array();
	while(false!==($sLn = fgets($hFilter) )) {
		if ($sLn == "") continue;
		$arrF = explode("\t", $sLn);

		//print_r($arrF);
		if (count($arrF) !=15) continue;

		$sCandidateAln = $arrF[14];
		preg_match("/ID=([^;]+);/" , $sCandidateAln , $arrM);
		if (count($arrM) != 2) continue;

		$sCandidateGeneID = trim($arrM[1]);
		//echo($sCandidateGeneID . "\n");
		$arrCanGenes[$sCandidateGeneID] = true;
	}

	//print_r($arrCanGenes);
	$hGFF = fopen($sGFF , "r");

	while(false!==($sLn = fgets($hGFF) )) {
		if ($sLn == "") continue;
		if ($sLn[0] =='#') continue;

		$arrF = explode("\t", trim( $sLn) );
		$arrF[1] = "exonerate";

		//print_r($arrF);
		if (count($arrF) < 5) continue;

		if ($arrF[2] == "mRNA") {
			preg_match("/ID=([^;]+);/" , $arrF[8]  , $arrM);
			if (count($arrM) != 2) continue;
			$sID = trim( $arrM[1] );
			$sGeneID = "gene.$sID";
		
			if (array_key_exists($sID  , $arrCanGenes)) {

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
		
			if (array_key_exists($sID  , $arrCanGenes)) {
				fwrite($hOut , implode("\t", $arrF)."\n" );
				$arrF[2] = "CDS";
				fwrite($hOut , implode("\t", $arrF)."\n" );

			}
			continue;
		}

	
	}

	return $sOut;

}

function fnMergeOverlap(&$arrCoord, $nExtend=0) {


		
		$nRegions = -1;
		//echo("======== $sContig ==========\n");
		while(true) {
			$arrCoord = fnMergeOverlapContig($arrCoord, $nExtend);
			//print_r($arrRegions[$sContig]);
			
			if ($nRegions == count($arrCoord)) {
				break;
			}
			else {
				$nRegions = count($arrCoord);
			}
		}
		

	

}

function fnMergeOverlapContig(&$arrCoord, $nExtend) {

		$arrMergedCoord = array();

		reset($arrCoord);
		while($nCurrEnd = current($arrCoord)) {
			$nCurrStart = key($arrCoord);
			$nNextEnd = next($arrCoord);
			$nNextStart = key($arrCoord);
			
			
			$nMergedStart = $nCurrStart;
			
			if ($nNextEnd === false ) {
				$nMergedEnd = $nCurrEnd;
				
			}
			else {
			
				if ($nCurrEnd + $nExtend < $nNextStart) { // current end before next start
					$nMergedEnd = $nCurrEnd; // keep this
					
				}
				else  {
					$nMergedEnd = max($nNextEnd ,  $nCurrEnd);
					next($arrCoord);
				}
				
			}
			
			$arrMergedCoord[$nMergedStart] = $nMergedEnd ;

			
		}
		
		return $arrMergedCoord;
}

function fnScoreExonerateModel($sPredictedProt , $sEvidenceProt , $sDIR ) {
	global $nMinBlastPIdentity , $nMaxBlastPIdentityDiff , $nTreatAsBigGap , $nTreatAsHugeGap; 

	//echo("$sPredictedProt\n\n$sEvidenceProt\n");
	$arrScores = array('score' => 0, 'internalstop' => false);

	if (substr($sPredictedProt , -1,1) == '*' ) {
		$sPredictedProt = substr($sPredictedProt , 0, strlen($sPredictedProt)-1 );
	}

	if (substr($sEvidenceProt , -1,1) == '*' ) {
		$sEvidenceProt = substr($sEvidenceProt , 0, strlen($sEvidenceProt)-1 );
	}


	if (strpos($sPredictedProt, '*') !== false ) {
		$arrScores['internalstop'] = true;
		//echo("stop found, discard\n");
		return $arrScores; //bad model, has internal stops
	}

	//echo("$sPredictedProt\n\n$sEvidenceProt\n");
	//echo("hi");

	$nPredLen = strlen($sPredictedProt);
	$nEvidLen = strlen($sEvidenceProt);

	$sPredFa = "$sDIR/_predprot.fa";
	$sEvidFa = "$sDIR/_evidprot.fa";

	$hPredFa = fopen($sPredFa , "w");
	$hEvidFa = fopen($sEvidFa , "w");

	fwrite($hPredFa , ">pred\n$sPredictedProt");
	fwrite($hEvidFa , ">evid\n$sEvidenceProt");

	$sBlastOut = "$sDIR/_scoremodel.blast.out";

	exec("makeblastdb -in $sPredFa -dbtype prot; blastp -task blastp -db $sPredFa -query $sEvidFa  -evalue 1e-10 -out $sBlastOut -outfmt '7 qseqid sseqid qstart qend sstart send evalue bitscore length pident mismatch gapopen btop'");

 	//now examine the alignment:
	$hBlastOut = fopen($sBlastOut, "r");
	$arrPredBlastRet = array(); //key is protein evidence ID, value is array of hits
	$arrEvidBlastRet = array(); //key is protein evidence ID, value is array of hits
	$arrPercIden = array();
	while(false !== ($sLn = fgets($hBlastOut) )) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		if ($sLn[0] == '#') continue;
		$arrF = explode("\t" , $sLn);
		if (count($arrF) != 13) continue;

		$nQStart = $arrF[2];
		$nQEnd = $arrF[3];
		$nHStart = $arrF[4];
		$nHEnd = $arrF[5];
		$nPercIden = $arrF[9];
		if ($nQStart > $nQEnd || $nHStart > $nHEnd ) {
			continue;
		}

		if ($nPercIden < $nMinBlastPIdentity) continue;


		$arrPredBlastRet[$nHStart] = array('qstart' => $nQStart, 'qend' =>$nQEnd, 'hstart' => $nHStart, 'hend' =>  $nHEnd, 'percid' => $nPercIden, 'btop' => $arrF[12]);
		$arrEvidBlastRet[$nQStart] = array('qstart' => $nQStart, 'qend' =>$nQEnd, 'hstart' => $nHStart, 'hend' =>  $nHEnd, 'percid' => $nPercIden, 'btop' => $arrF[12]);
		$arrPercIden[] = $nPercIden;
	}

	if (count($arrPercIden) ==0 ) {
		return $arrScores; // no blast results
	}

	$nMinPercIden = max($arrPercIden) - $nMaxBlastPIdentityDiff;
	ksort($arrPredBlastRet);
	ksort($arrEvidBlastRet);

	$arrPredScores = array_fill(1 , $nPredLen , 0); //initialize scores 
	$arrPredCov = array_fill(1 , $nPredLen , 0); //initialize coverage, this is supposed to be 1 

	foreach($arrPredBlastRet as $nAlnStart => &$oAlnInfo ) {
		if ($oAlnInfo['percid'] < $nMinPercIden) {
			continue;
		}

		preg_match_all("/\d+|[A-Z-][A-Z-]/", $oAlnInfo['btop'] , $arrBTOP);

		$nCurrAlnPos = $nAlnStart;

		foreach($arrBTOP[0] as &$sBTOP) {
			if (intval($sBTOP)==0) { //substitution or gap
				if ($sBTOP[0] == '-') { //insertion
					$arrPredScores[$nCurrAlnPos] -= 1;
					$nCurrAlnPos++;
					continue;
				}
				if ($sBTOP[1] == '-') { //deletion do nothing
					continue;
				}
					//substitution, no score
				if ($arrPredCov[$nCurrAlnPos] !== '-' && $arrPredCov[$nCurrAlnPos] > 0) {
					$arrPredCov[$nCurrAlnPos] += 1;//count as coverage, but do not add score.
				} else if ($arrPredCov[$nCurrAlnPos] === 0) {
					$arrPredCov[$nCurrAlnPos] = '-'; // count as gap
				}
 
				$nCurrAlnPos++;

				continue;
			} else {
				$nMatch = intval($sBTOP);
				$nMatchEnd = $nCurrAlnPos+$nMatch-1;
				for( ; $nCurrAlnPos<=$nMatchEnd; $nCurrAlnPos++) {
					$arrPredScores[$nCurrAlnPos] += 1;
					if ($arrPredCov[$nCurrAlnPos] === '-') {$arrPredCov[$nCurrAlnPos] = 1;}
					$arrPredCov[$nCurrAlnPos] += ($arrPredCov[$nCurrAlnPos] <9)? 1:0;
				}
			}
		}
	}

	//print_r($arrPredScores);
	//print_r($arrPredCov);


	$arrEvidScores = array_fill(1 , $nEvidLen , 0); //initialize scores 
	$arrEvidCov = array_fill(1 , $nEvidLen , 0); //initialize coverage, this is supposed to be 1 

	foreach($arrEvidBlastRet as $nAlnStart => &$oAlnInfo ) {
		if ($oAlnInfo['percid'] < $nMinPercIden) {
			continue;
		}

		preg_match_all("/\d+|[A-Z-][A-Z-]/", $oAlnInfo['btop'] , $arrBTOP);

		$nCurrAlnPos = $nAlnStart;
		foreach($arrBTOP[0] as &$sBTOP) {
			if (intval($sBTOP)==0) { //substitution or gap
				$sBTOP = strrev($sBTOP);
				if ($sBTOP[0] == '-') { //insertion
					$arrEvidScores[$nCurrAlnPos] -= 1;
					$nCurrAlnPos++;
					continue;
				}
				if ($sBTOP[1] == '-') { //deletion do nothing
					continue;
				}

					//substitution, no score
				if ( $arrEvidCov[$nCurrAlnPos] !== '-' && $arrEvidCov[$nCurrAlnPos] > 0) {
					$arrEvidCov[$nCurrAlnPos] += 1;//count as coverage, but do not add score.
				} else if ($arrEvidCov[$nCurrAlnPos] === 0) {
					$arrEvidCov[$nCurrAlnPos] = '-'; // count as mismatch
				}
 
				$nCurrAlnPos++;
				continue;
			} else {
				$nMatch = intval($sBTOP);
				$nMatchEnd = $nCurrAlnPos+$nMatch-1;
				for( ; $nCurrAlnPos<=$nMatchEnd; $nCurrAlnPos++) {
					$arrEvidScores[$nCurrAlnPos] += 1;
					if ($arrEvidCov[$nCurrAlnPos] === '-') { $arrEvidCov[$nCurrAlnPos] = 1;}
					$arrEvidCov[$nCurrAlnPos] += ($arrEvidCov[$nCurrAlnPos] <9)? 1:0; 
				}
			}
		}
	}

	//print_r($arrEvidScores);
	//print_r($arrEvidCov);

	$nEvidScore = array_sum($arrEvidScores);
	$nPredScore = array_sum($arrPredScores);

	echo("\t\t\t\tEvidence score: $nEvidScore ; Prediction score: $nPredScore\n");

	$oEvidCovScore = fnCheckCov($arrEvidCov);
	$oPredCovScore = fnCheckCov($arrPredCov);

	//print_r($oEvidCovScore);
	//print_r($oPredCovScore);
	$arrScores['score'] = $nEvidScore + $nPredScore + $oEvidCovScore['gappenalty'] + $oPredCovScore['gappenalty'];
	$arrScores['5primemiss'] = $oEvidCovScore['5primemiss'];
	$arrScores['5primemissperc'] = $oEvidCovScore['5primemissperc'];
	$arrScores['3primemiss'] = $oEvidCovScore['3primemiss'];
	$arrScores['3primemissperc'] = $oEvidCovScore['3primemissperc'];
	$arrScores['5primeuncertain'] = $oPredCovScore['5primemiss'];
	$arrScores['3primeuncertain'] = $oPredCovScore['3primemiss'];
	$arrScores['completeness'] = $oEvidCovScore['singlecoverage'];
	$arrScores['hicov'] = $oEvidCovScore['hicoverage'];

	$arrScores['smalldeletionlen'] = $oEvidCovScore['smallgaplen'];
	$arrScores['smalldeletioncount'] = $oEvidCovScore['smallgapcount'];
	$arrScores['bigdeletionlen'] = $oEvidCovScore['biggaplen'];
	$arrScores['bigdeletioncount'] = $oEvidCovScore['biggapcount'];
	$arrScores['hugedeletionlen'] = $oEvidCovScore['hugegaplen'];
	$arrScores['hugedeletioncount'] = $oEvidCovScore['hugegapcount'];


	$arrScores['smallinsertionlen'] = $oPredCovScore['smallgaplen'];
	$arrScores['smallinsertioncount'] = $oPredCovScore['smallgapcount'];
	$arrScores['biginsertionlen'] = $oPredCovScore['biggaplen'];
	$arrScores['biginsertioncount'] = $oPredCovScore['biggapcount'];
	$arrScores['biginsertionlist'] = $oPredCovScore['biggaps'];
	$arrScores['hugeinsertionlen'] = $oPredCovScore['hugegaplen'];
	$arrScores['hugeinsertioncount'] = $oPredCovScore['hugegapcount'];
	$arrScores['hugeinsertionlist'] = $oPredCovScore['hugegaps'];
	//$arrScores['btop'] = $oAlnInfo['btop'];

	return $arrScores;
}

function fnCheckCov(&$arrEvidCov) {
	global $nTreatAsBigGap , $nTreatAsHugeGap; 

	$nFullLen = count($arrEvidCov);
	$sEvidCov = implode('' , $arrEvidCov);
	//echo("Before: $sEvidCov\n");
	$sEvidCov = fnSmoothCov($sEvidCov);
	//echo("After: $sEvidCov\n");
	$n5PrimeMiss = 0;
	$n3PrimeMiss = 0;
	$arrInternalSmallGapLen = array();
	$arrInternalLargeGapLen = array();
	$arrInternalHugeGapLen = array();
	$arrHiCovLen = array();

	$arrInternalHugeGapList = array();
	$arrInternalLargeGapList = array();


	preg_match("/^0+/", $sEvidCov, $arrM); //check coverage at 5'
	if (count($arrM) > 0) {
		$n5PrimeMiss = strlen($arrM[0]);
	}

	preg_match("/0+$/", $sEvidCov, $arrM); //check coverage at 3'
	if (count($arrM) > 0) {
		$n3PrimeMiss = strlen($arrM[0]);
	}

	preg_match_all("/0+/", $sEvidCov, $arrM, PREG_OFFSET_CAPTURE); //find internal gaps (cov = 0)
	if (count($arrM[0]) > 0) {
		foreach($arrM[0] as &$arrGap) {
			$nGapLen = strlen($arrGap[0]);
			$nGapPos = $arrGap[1];
			if ($nGapPos == 0 || ($nGapPos+$nGapLen) >= $nFullLen ) continue;

			if ($nGapLen > $nTreatAsHugeGap) {
				$arrInternalHugeGapLen[] = $nGapLen;
				$arrInternalHugeGapList[] = ($nGapPos+1) . "-". ($nGapPos + $nGapLen);
			} else if ($nGapLen > $nTreatAsBigGap) {
				$arrInternalLargeGapLen[] = $nGapLen;
				$arrInternalLargeGapList[] = ($nGapPos+1) . "-". ($nGapPos + $nGapLen);
			} else {
				$arrInternalSmallGapLen[] = $nGapLen;
			}
		}
	}

	preg_match_all("/[^01]+/", $sEvidCov, $arrM); //find internal gaps (cov = 0)
	if (count($arrM[0]) > 0) {
		foreach($arrM[0] as &$sHiCov) {
			$arrHiCovLen[] = strlen($sHiCov);			
		}
	}

	$nSmallGapLen = array_sum($arrInternalSmallGapLen);
	$nHugeGapLen = array_sum($arrInternalHugeGapLen);
	$nLargeGapLen = array_sum($arrInternalLargeGapLen);

	$nSmallGapCount = count($arrInternalSmallGapLen);
	$nHugeGapCount = count($arrInternalHugeGapLen);
	$nLargeGapCount = count($arrInternalLargeGapLen);

	$nHiCovLen = array_sum($arrHiCovLen);

	$nGapPenalty = 0;
	$nGapPenalty = $nHugeGapLen * 3 + $nLargeGapLen * 2; //extra penalty for gaps

	return array('5primemiss' => $n5PrimeMiss , '5primemissperc' => $n5PrimeMiss / $nFullLen, 
		     '3primemiss' => $n3PrimeMiss, '3primemissperc' => $n3PrimeMiss / $nFullLen, 
		     'gappenalty' => -$nGapPenalty , 
		     'smallgaplen' => $nSmallGapLen , 'smallgapcount' => $nSmallGapCount, 
		     'biggaplen' => $nLargeGapLen , 'biggapcount' => $nLargeGapCount, 
		     'hugegaplen' => $nHugeGapLen , 'hugegapcount' => $nHugeGapCount,
		     'singlecoverage' => ($nFullLen - $n5PrimeMiss - $n3PrimeMiss - $nLargeGapLen - $nHugeGapLen - $nHiCovLen) / $nFullLen, 
		     'hicoverage' => $nHiCovLen / $nFullLen,
		     'hugegaps' => implode(',' , $arrInternalHugeGapList),
		     'biggaps' => implode(',' , $arrInternalLargeGapList)
		);
	
}

function fnSmoothCov($s) {

	$sRet = $s;
	// first, replace all mismatches around gaps as gaps ---0000000--- into 00000000000
	preg_match_all('/-+0+-+/', $sRet, $arrM, PREG_OFFSET_CAPTURE);

	foreach($arrM[0] as $oMatch ) {
		$nMLen = strlen($oMatch[0]);
		$nOffSet = $oMatch[1];
		$sRet = substr_replace($sRet , str_repeat('0' , $nMLen) , $nOffSet ,  $nMLen);
	}

	//now use a window to scan
	$nWinSize = 10;
	$arrLeftEdge = array();
	$nEnd = strlen($sRet)-$nWinSize;
	for($i=0;$i<$nEnd;$i++) {
		$sWin = substr($sRet , $i, $nWinSize);
		$nBad = substr_count($sWin , '0');
		$nBad += substr_count($sWin , '-') * 0.5;
		$arrLeftEdge[$i] = 1- ($nBad / $nWinSize);
	}

	$arrRightEdge = array();
	$nEnd = strlen($sRet)-1;
	for($i=$nEnd;$i>=($nWinSize-1);$i--) {
		$sWin = substr($sRet , $i - $nWinSize + 1, $nWinSize);
		$nBad = substr_count($sWin , '0');
		$nBad += substr_count($sWin , '-') * 0.5;
		$arrRightEdge[$i] = 1- ($nBad / $nWinSize);
	}

	//print_r($arrLeftEdge);
	//print_r($arrRightEdge); die();

	for($i=0;$i<strlen($sRet);$i++) {
		if ($sRet[$i]=='0') continue;
		$nLeftEdgeScore = array_key_exists($i, $arrLeftEdge)? $arrLeftEdge[$i]: 0;
		$nRightEdgeScore = array_key_exists($i, $arrRightEdge)? $arrRightEdge[$i]: 0;
		$nScore = min($nLeftEdgeScore , $nRightEdgeScore);
		if ( intval($sRet[$i]) >= 1) {
			$nScore = max($nLeftEdgeScore , $nRightEdgeScore);
		}
		if ($nScore >= 0.7) {
			$sRet[$i]  = $sRet[$i];
		} else {
			$sRet[$i]  = '0';
		}
	}

	return  str_replace('-','0', $sRet);
}

function fnWriteNewModel($arrNewModel, $oScore, $hFixedGFF, $sPredictor='maker') {
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
	fwrite($hFixedGFF , implode("\t" , $arrF)."\n");
	$arrF[2] = 'mRNA';
	$arrF[5] = $oScore['score'];
	$arrF[8] = fnArr2Annot($arrNewModel['annot']);
	$arrF[8] .= ";". fnArr2Annot($oScore);
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

function fnSafeID($s) {
	return urlencode($s);
}

function fnUnSafeID($s) {
	return urldecode($s);
}

?>
