<?php
require_once(dirname(__FILE__) . "/lib.php");

///////////////////////////////////////

$sMakerGFF = "split__main.gff";
//$sMakerGFF = "problem_model.gff";
$sGenome = "query.fa";
$sProteinGFF = "split__protein2genome.gff";
$sProteinFasta = "all.protein.evidence.fa";
$sCandidiateGFF = "split__candidate.gff"; // this is the candidate gene models; the genewise "pred_gff_genewise" models will be used. other models ignored.
//$sWorkDir = "./splitchimeric/";
$sWorkDir = "./splitchimeric_improved/";
$nMinBlastPIdentity = 50; //filter HSP if lower than this.
$nMaxBlastPIdentityDiff = 20; //If identity of an HSP block is 20% smaller than other blocks in the same alignment, then exclude this hsp
$nTreatAsBigGap = 2; // if larger than this, treat as "big gap".
$nTreatAsHugeGap = 10; // if larger than this, treat as "huge gap".

$nThisPart = 0;
$nTotalParts = 1;

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
	case '-p':
	    $sProteinFasta  = trim(array_shift($argv));
	    break;
	case '-C':
	    $sCandidiateGFF  = trim(array_shift($argv));
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

	//if (md5($sGeneID) != '94c9d1b9e532eb714497555a9515ac52') continue;
	//if (md5($sGeneID) != 'c0fda4466793b5b2589689892e2dfb84') continue;



	echo("Looking at $sGeneID ...\n");

	$sGeneDIR = "$sWorkDir/".$oGene['scf']."/".md5($sGeneID);

	$oGeneRevised = $oGene;
	$bRevised = false; // flag for whether a revision has been made

	foreach($oGene['mRNAs'] as $smRNAID => &$omRNA) {

		$sRNADIR = $sGeneDIR."/".md5($smRNAID);
		exec("mkdir -p $sRNADIR");
		echo(" -- mRNA  $smRNAID in dir $sRNADIR ...\n");
		$arrRNARange = array( array($omRNA['scf'] , $omRNA['start'], $omRNA['end'] , $omRNA['strand']) );
		$sLocalProteinGFF = fnGetOverlapGFF($sRNADIR , $sProteinGFF, $sCandidiateGFF, $arrRNARange, 0.01, $nPercOverlapB=0.01);
		$oLocalProteinGFF = new MappedCDSGFF();
		$oLocalProteinGFF->LoadGFF($sLocalProteinGFF , $sGenome);

		$sProteinSeqFa = "$sRNADIR/prot.evidence.fa"; //this is the original NCBI proteins used as evidence
		$arrProtLength = array(); //this writes down the length of the proteins
		$hProteinSeqFa = fopen($sProteinSeqFa , 'w');

		foreach($oLocalProteinGFF->GenePredictors() as $sPredictor) {
			$arrLocalProteinGenes = $oLocalProteinGFF->GetGenes($sPredictor);
			foreach( $arrLocalProteinGenes  as $sProteinGFFID => &$arrProteinModels ) {
				foreach($arrProteinModels['mRNAs'] as &$oProteinModel) {
					$sProteinQueryName = $oProteinModel['annot']['Name'];
					$sProteinID = fnSafeID( $oProteinModel['annot']['ID'] );//needed for fastahack
					if ( $sPredictor != 'genewise' ) { //if it is not a genewise alignment, use the original protein
						//echo("$sProteinQueryName\n");
						$sProteinSeq = $oProteinFasta->GetContig($sProteinQueryName);
						//echo("$sProteinSeq\n");
						$arrProtLength[$sProteinID] = strlen($sProteinSeq);
						fwrite($hProteinSeqFa , ">$sProteinID\n$sProteinSeq\n");
					} else {
						$oGenewiseSeq = $oLocalProteinGFF->ExtractmRNASequenceByModel( $oProteinModel);
						$arrProtLength[$sProteinID] = strlen($oGenewiseSeq['AA']);
						
							//print_r($oProteinModel);
							//print_r($oGenewiseSeq);

						fwrite($hProteinSeqFa , ">$sProteinID\n".$oGenewiseSeq['AA']."\n");
					}
				}
			}
		}

		if ( count($arrProtLength) == 0) {
			// if nothing aligns
			echo(" --- no evidence found, skip.\n");
			continue;
		}

		$arrExtracted = $oMakerGFF->ExtractmRNASequence('maker', $smRNAID);
		$sGenemodelFasta = "$sRNADIR/predicted.fa";
		$hGenemodelFasta = fopen($sGenemodelFasta , 'w');
		fwrite($hGenemodelFasta , ">$smRNAID\n".$arrExtracted['AA']."\n");

		// now blast each evidence protein against the predicted protein:
		$sBlastOut = "$sRNADIR/blast.out";
		exec("if [ ! -e $sBlastOut.ok ]; then makeblastdb -in $sGenemodelFasta -dbtype prot; blastp -task blastp -db $sGenemodelFasta -query $sProteinSeqFa  -evalue 1e-10 -out $sBlastOut -outfmt '7 qseqid sseqid qstart qend sstart send evalue bitscore length pident mismatch gapopen btop' && touch $sBlastOut.ok; fi");
		$arrClassified = fnClassifyGeneModel($sBlastOut , $arrProtLength);

		if ($arrClassified === false ) {
			echo("\tcannot classify. skip\n");
			continue;
		}
		//print_r($oLocalProteinGFF);

		if (!$arrClassified['isfused']) {
			echo("\tmRNA ID $smRNAID is ok...\n");
		}

		if ($arrClassified['isfused']) {
			echo("\tmRNA ID $smRNAID is a chimeric model from ".$arrClassified['pieces']." genes\n\t\tSelect from the best exonerate models to promote...\n");
			//print_r($arrClassified);
			//scan the exonerate alignments for the best model to promote:
			$oEvidenceProtFasta = new FastaHack();
			$oEvidenceProtFasta->SetDb($sProteinSeqFa);

			$hFixedGFF = fopen("$sRNADIR/fixed.gff", 'w');

			//go over sub fragment by sub fragment:
			foreach($arrClassified['partitioned_proteins'] as $nParStart => &$oPartition) {
				echo("\t\tPartition starting at: $nParStart \n");

				$nCurrentHighScore = 0;
				$oCurrentBestModel = array();
				foreach($oPartition as &$sEvidenceID) {
					//examine evidence:

					echo("\t\t\t $sEvidenceID \n");
					if ($oLocalProteinGFF->mRNAIDExistsInPredictor('exonerate', fnUnSafeID($sEvidenceID))) {
						$arrPredictionExtracted = $oLocalProteinGFF->ExtractmRNASequence('exonerate', fnUnSafeID($sEvidenceID) );
					
						$oScore = fnScoreExonerateModel($arrPredictionExtracted['AA'] , $oEvidenceProtFasta->GetContig($sEvidenceID) , $sRNADIR );
						//print_r($oScore);
						//die();
						if ($oScore['internalstop']) continue;
						if ($oScore['score'] > $nCurrentHighScore) {
							$nCurrentHighScore = $oScore['score'];
							$oCurrentBestModel =  $oScore;
							$oCurrentBestModel['id'] = $sEvidenceID;
						}
					} else if ($oLocalProteinGFF->mRNAIDExistsInPredictor('genewise', fnUnSafeID($sEvidenceID))) {
						$arrPredictionExtracted = $oLocalProteinGFF->ExtractmRNASequence('genewise', fnUnSafeID($sEvidenceID) );
						$oScore = fnScoreExonerateModel($arrPredictionExtracted['AA'] , $arrPredictionExtracted['AA'] , $sRNADIR ); //aligns to itself, scores not meaningful.
						if ($oScore['internalstop']) continue;
						if ($oScore['score'] > $nCurrentHighScore) {
							$nCurrentHighScore = $oScore['score'];
							$oCurrentBestModel =  $oScore;
							$oCurrentBestModel['id'] = $sEvidenceID;
						}
					}

				}

				if ($nCurrentHighScore == 0) {
					echo("\t\t\t\tNo good model found for this sub fragment.\n");
					continue;
				}

				echo("\t\t\t\tPromote model " . $oCurrentBestModel['id'] ." ; score = ".$oCurrentBestModel['score']."\n");
				$bRevised = true;
				unset($oGeneRevised['mRNAs'][$smRNAID]); //delete the current mRNA.
				$sNewmGeneID = "SplitChimeric_$nParStart.$sGeneID";
				$sNewmRNAID = "SplitChimeric_$nParStart.$smRNAID";
				$arrNewModel = ($oLocalProteinGFF->mRNAIDExistsInPredictor('exonerate', fnUnSafeID($oCurrentBestModel['id'])))? 
						$oLocalProteinGFF->GetmRNA('exonerate' , fnUnSafeID($oCurrentBestModel['id'])) : 
						$oLocalProteinGFF->GetmRNA('genewise' , fnUnSafeID($oCurrentBestModel['id']))
						;
				$arrNewModel['annot']['Parent'] = $sNewmGeneID; //set it to the current gene model.
				$arrNewModel['annot']['ID'] = $sNewmRNAID; //set it to the current gene model.
				//print_r($arrNewModel);
				fnWriteNewModel($arrNewModel, $oCurrentBestModel, $hFixedGFF);
			}
		}
	}
}

function fnGetOverlapGFF($sTmpDir , $sGFF, $sCandidateGFF, $arrRanges, $nPercOverlapA=0.01, $nPercOverlapB=0.01) {

	$sFilter = "$sTmpDir/overlapped.genes.txt";
	$sFilterGenewise = "$sTmpDir/overlapped.genewise.txt";

	$sOut = "$sTmpDir/proteinsubset.gff";

	$sBed = "$sTmpDir/ranges.bed";
	$hBed = fopen($sBed, "w");

	foreach($arrRanges as $oRange) {
		fwrite($hBed , $oRange[0]."\t".($oRange[1]-1)."\t$oRange[2]\tA\t1000\t$oRange[3]\n"); 
	}

	exec("grep -P \"\tmRNA\t\" $sGFF > $sTmpDir/mRNA.gff");
	exec("bedtools intersect -a $sBed -b $sTmpDir/mRNA.gff -wb -s -f $nPercOverlapA -F $nPercOverlapB > $sFilter; rm $sTmpDir/mRNA.gff");

	exec("grep -P \"\tmRNA\t\" $sCandidateGFF > $sTmpDir/mRNA.gff");
	exec("bedtools intersect -a $sBed -b $sTmpDir/mRNA.gff -wb -s -f $nPercOverlapA -F $nPercOverlapB | grep pred_gff_genewise > $sFilterGenewise; rm $sTmpDir/mRNA.gff");


	$hFilter = popen("cat $sFilter $sFilterGenewise", "r");
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
	$hGFF = popen( "cat $sGFF $sCandidateGFF", "r");

	while(false!==($sLn = fgets($hGFF) )) {
		if ($sLn == "") continue;
		if ($sLn[0] =='#') continue;

		$arrF = explode("\t", trim( $sLn) );
		$arrF[1] = (strpos($arrF[1],'genewise') ===false )? "exonerate":"genewise";

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

function fnClassifyGeneModel($sBlastRet, $arrEvidenceProtLength) { //based on blast results, check if Gene model is likely a fused gene. if so, where is the break?
	global $nMinBlastPIdentity;
	$sWD = dirname($sBlastRet);
	$sCoordTab = "$sWD/_coord.tab";
	$sClusterTab = "$sWD/_clusters.tab";
	$sClusterMedTab = "$sWD/_cluster_median.tab";

	$hCoordTab = fopen($sCoordTab, "w");
	$hBlastRet = fopen($sBlastRet , "r");
	$arrBlastRet = array(); //key is protein evidence ID, value is array of hits
	while(false !== ($sLn = fgets($hBlastRet) )) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		if ($sLn[0] == '#') continue;
		$arrF = explode("\t" , $sLn);
		if (count($arrF) != 13) continue;

		$sQuery = $arrF[0];
		$nQStart = $arrF[2];
		$nQEnd = $arrF[3];
		$nHStart = $arrF[4];
		$nHEnd = $arrF[5];
		$nPercIden = $arrF[9];
		if ($nQStart > $nQEnd || $nHStart > $nHEnd ) {
			continue;
		}

		if ($nPercIden < $nMinBlastPIdentity) continue;

		if (!array_key_exists($sQuery , $arrBlastRet)) {
			$arrBlastRet[$sQuery] = array();
		}

		$arrBlastRet[$sQuery][$nHStart] = array('qstart' => $nQStart, 'qend' =>$nQEnd, 'hstart' => $nHStart, 'hend' =>  $nHEnd, 'percid' => $nPercIden);
	}

	if (count($arrBlastRet) == 0) { //no blast hits
		return false;
	}

	//print_r($arrBlastRet);

	$arrHitPoints = array();


	$arrHitStartEnd = array();
	$arrBlastRetRange = array(); //within each query, contains the min and max hit coordinates

	fwrite($hCoordTab , "Start\tEnd\tID\n");
	foreach($arrBlastRet as $sQuery => &$arrBlastRecord) {
		ksort($arrBlastRecord ,SORT_NUMERIC);
		$nMinStart = -1;
		$nMaxEnd = -1;
		$arrBlastRetRange[$sQuery] = array(0,0);

		foreach($arrBlastRecord as $nHStart => $arrInfo) {
			if ($nMinStart == -1) {
				$nMinStart = $nHStart ; //the first one is the smallest
				$arrBlastRetRange[$sQuery][0] = $nHStart ;
			}
			if ($nMaxEnd == -1) {
				$nMaxEnd = $arrInfo['hend'];
				$arrBlastRetRange[$sQuery][1] = $arrInfo['hend'];
			} else if ($nMaxEnd < $arrInfo['hend']) {
				$nMaxEnd = $arrInfo['hend'];
				$arrBlastRetRange[$sQuery][1] = $arrInfo['hend'];
			}
		} 

		$arrHitStartEnd[] = array($nMinStart, $nMaxEnd);
		fwrite($hCoordTab , "$nMinStart\t$nMaxEnd\t$sQuery\n");

	}

	fclose($hCoordTab);

	exec("Rscript coordcluster.R $sCoordTab $sClusterTab $sClusterMedTab");

	$hClusterTab = fopen($sClusterTab , 'r');
	$hClusterMedTab = fopen($sClusterMedTab , 'r'); //calculated median start and end positions of the cluster

	$arrClusters = array();
	$arrPartitionedGeneID = array();
	$arrHitStartEnd = array();

	while( false !== ($sLn = fgets($hClusterMedTab) ) ) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		$arrF = explode("\t" , $sLn);
		if ($arrF[0] == 'clu') continue;
		$arrClusters[$arrF[0]] = array($arrF[1], $arrF[2]);
		$arrPartitionedGeneID[$arrF[1]] = array();
		$arrHitStartEnd[$arrF[1]] = $arrF[2];
	}

	$nSubFragments = count($arrClusters);

	$bFused = ($nSubFragments > 1);

	while( false !== ($sLn = fgets($hClusterTab) ) ) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		$arrF = explode("\t" , $sLn);
		if ($arrF[0] == 'Start') continue;
		if ($arrF[4] != 1) continue ; //doesn't have the "keep" flag
		$nClusterStart = $arrClusters[$arrF[3]][0];
		$arrPartitionedGeneID[$nClusterStart][] = $arrF[2];
	}
/*
	fnMergeOverlap($arrHitStartEnd , -10);
	//print_r($arrHitStartEnd);

	$bFused = false;
	$nSubFragments = count($arrHitStartEnd);

	if ($nSubFragments > 1 ) {
		$bFused = true;
	}

	//calculate the proteins belonging to each sub fragment

	$arrPartitionedGeneID = array();

	foreach($arrHitStartEnd as $nStart => $nEnd) {
		$arrPartitionedGeneID[$nStart] = array();
		foreach($arrBlastRetRange as $sQuery => $arrRange) {
			if ($arrRange[0] >= $nStart && $arrRange[1] <= $nEnd ) {
				$arrPartitionedGeneID[$nStart][] = $sQuery;
			}
		}
	}
*/

	return array('isfused' => $bFused , 'pieces'=> $nSubFragments,  'breaks' => $arrHitStartEnd, 'partitioned_proteins' => $arrPartitionedGeneID);
	
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

				$arrPredCov[$nCurrAlnPos] += 1;//count as coverage, but do not add score.
				$nCurrAlnPos++;

				continue;
			} else {
				$nMatch = intval($sBTOP);
				$nMatchEnd = $nCurrAlnPos+$nMatch-1;
				for( ; $nCurrAlnPos<=$nMatchEnd; $nCurrAlnPos++) {
					$arrPredScores[$nCurrAlnPos] += 1;
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
				$arrEvidCov[$nCurrAlnPos] += 1;//but count as coverage
				$nCurrAlnPos++;
				continue;
			} else {
				$nMatch = intval($sBTOP);
				$nMatchEnd = $nCurrAlnPos+$nMatch-1;
				for( ; $nCurrAlnPos<=$nMatchEnd; $nCurrAlnPos++) {
					$arrEvidScores[$nCurrAlnPos] += 1;
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
	$arrScores['smallinsertionlen'] = $oEvidCovScore['smallgaplen'];
	$arrScores['smallinsertioncount'] = $oEvidCovScore['smallgapcount'];
	$arrScores['biginsertionlen'] = $oEvidCovScore['biggaplen'];
	$arrScores['biginsertioncount'] = $oEvidCovScore['biggapcount'];
	$arrScores['hugeinsertionlen'] = $oEvidCovScore['hugegaplen'];
	$arrScores['hugeinsertioncount'] = $oEvidCovScore['hugegapcount'];

	$arrScores['smalldeletionlen'] = $oPredCovScore['smallgaplen'];
	$arrScores['smalldeletioncount'] = $oPredCovScore['smallgapcount'];
	$arrScores['bigdeletionlen'] = $oPredCovScore['biggaplen'];
	$arrScores['bigdeletioncount'] = $oPredCovScore['biggapcount'];
	$arrScores['hugedeletionlen'] = $oPredCovScore['hugegaplen'];
	$arrScores['hugedeletioncount'] = $oPredCovScore['hugegapcount'];

	return $arrScores;
}

function fnCheckCov(&$arrEvidCov) {
	global $nTreatAsBigGap , $nTreatAsHugeGap; 

	$nFullLen = count($arrEvidCov);
	$sEvidCov = implode('' , $arrEvidCov);
	$n5PrimeMiss = 0;
	$n3PrimeMiss = 0;
	$arrInternalSmallGapLen = array();
	$arrInternalLargeGapLen = array();
	$arrInternalHugeGapLen = array();
	$arrHiCovLen = array();


	preg_match("/^0+/", $sEvidCov, $arrM); //check coverage at 5'
	if (count($arrM) > 0) {
		$n5PrimeMiss = strlen($arrM[0]);
	}

	preg_match("/0+$/", $sEvidCov, $arrM); //check coverage at 5'
	if (count($arrM) > 0) {
		$n3PrimeMiss = strlen($arrM[0]);
	}

	preg_match_all("/[^0](0+)[^0]/", $sEvidCov, $arrM); //find internal gaps (cov = 0)
	if (count($arrM[1]) > 0) {
		foreach($arrM[1] as &$sGap) {
			$nGapLen = strlen($sGap);

			if ($nGapLen > $nTreatAsHugeGap) {
				$arrInternalHugeGapLen[] = $nGapLen;
			} else if ($nGapLen > $nTreatAsBigGap) {
				$arrInternalLargeGapLen[] = $nGapLen;
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
		     'hicoverage' => $nHiCovLen / $nFullLen
		);
	
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
