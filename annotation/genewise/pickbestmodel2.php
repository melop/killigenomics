<?php
//this version considers all fragments identified for the gene alignment, instead of only taking the longest one.
$arrOpts = getopt("g:f:o:k:s:S:");


$sGeneFamily = $arrOpts["g"]; // gene name"v2r";
//$sGeneWiseGTF = "out/Xmac/splitted/all.ret.fasta";
//$sTranslatedGTF = "sGeneFamily_Genewise_translated.gtf";
$sGeneWiseGTF = $arrOpts["f"]; //"out/medaka/splitted/all.ret.fasta";
$sSummaryFile = $arrOpts["s"]; //summary file, containing all the contiguous hits
$sTranslatedGTF =  $arrOpts["o"];//"sGeneFamily_Genewise_translated_medaka.gtf";

$bRemoveOverlap = true;

$hGeneWiseGTF = fopen($sGeneWiseGTF , "r");
$hSummaryFile = fopen($sSummaryFile , "r");
$hTranslatedGTF = fopen($sTranslatedGTF , "w");
$nGeneCount = 0;
$nExonCount = 0;
$nTotalExonCount = 0;
$arrAllGenes = array(); //key is chromosome name, secondary array are genes one chromosome

$nKeepHowMany = intval($arrOpts["k"]);//1;
$nKeepMaxLogRatio = $arrOpts["S"];// keep the hit if the score is this percentage of the best hit. say 0.5 means that it needs to at least have 50% of the best score
$nMinBitScore = -200;
//echo($nKeepScorePercent);

$arrTopScores = array(); //key-score, value - unused.

echo("Loading summary file...\n");

while(false!==($sLn=fgets($hSummaryFile))) {
	$sLn = trim($sLn);
	if ($sLn == "") {
		continue;
	}

	if (substr($sLn , 0, 4) == "Bits" ) {
		continue;
	}

	$reMatch = "/(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)/"; 
	$arrMatches = array();
	preg_match($reMatch, $sLn, $arrMatches);

	if (count( $arrMatches ) != 10) {
		continue;
	}


	list($sChr, $nStart, $nEnd) = fnParseSegmentCoord($arrMatches[5]);
	$nStart = min($nStart , $nEnd);
	$nLocalStart = min($arrMatches[6], $arrMatches[7]);
	$nLocalEnd   = max($arrMatches[6], $arrMatches[7]);
	$nTranslatedStart = $nStart -1 + $nLocalStart;
	$nTranslatedEnd	  = $nStart -1 + $nLocalEnd;


	$oGeneHit = array('Contig' => $sChr , 'Start' => $nTranslatedStart , 'End' => $nTranslatedEnd );
	if (count($arrTopScores) <  $nKeepHowMany) {
		$arrTopScores[$arrMatches[1]] = $oGeneHit ; // write down score.
	} else {

		$nToUnSet = false;
		foreach($arrTopScores as $nScore => $nval) {
			if ($nScore < $arrMatches[1]) {
				$nToUnSet = $nScore; 
				break;
			}
		}

		if ($nToUnSet !== false) {
			unset($arrTopScores[$nToUnSet]);
			$arrTopScores[$arrMatches[1]] = $oGeneHit; // write down score.
			ksort($arrTopScores , SORT_NUMERIC);
		}
	}
  
}

$arrTopScoresFiltered = array();

$nMaxScore = max(array_keys($arrTopScores));
if ( $nMaxScore >= $nMinBitScore ) {

	foreach($arrTopScores as $nScore => $nval) {
		//echo($nScore . " ".$nMaxScore);
		$nLogRatio =  $nMaxScore - $nScore ;
		if ( ($nLogRatio  <= $nKeepMaxLogRatio ) && ($nScore >= $nMinBitScore) ) {
			$arrTopScoresFiltered[$nScore] = $nval;
		}
	}
}

$arrTopScores = $arrTopScoresFiltered;

echo("Loading genewise GTF file...\n");
while(false!==($sLn=fgets($hGeneWiseGTF))) {
	$sLn = trim($sLn);
	if ($sLn == "") {
		continue;
	}
	
	$arrFields = explode("\t", $sLn);

	
	list($sChr, $nStart, $nEnd) = fnParseSegmentCoord($arrFields[0]);
	$nStart = min($nStart , $nEnd);
	$nLocalStart = min($arrFields[3], $arrFields[4]);
	$nLocalEnd   = max($arrFields[3], $arrFields[4]);
	$nTranslatedStart = $nStart -1 + $nLocalStart;
	$nTranslatedEnd	  = $nStart -1 + $nLocalEnd;
	$sOrientation = $arrFields[6];
	$nPhase = $arrFields[7];
	$nExonCount++;
	$nTotalExonCount++;
	
	
	if ($arrFields[2] == "match") {
		$nGeneCount++;
		$nExonCount = 0;
		
		if (!array_key_exists($sChr, $arrAllGenes)) {
			$arrAllGenes[$sChr] = array();
		}
		
		$arrAllGenes[$sChr][] = array(); // this array contains each exon
		$arrAllGenes[$sChr][count($arrAllGenes[$sChr])-1][0] = array("Tstart" => $nTranslatedStart, "Tend" => $nTranslatedEnd, "Orient" => $sOrientation, "Score" => $arrFields[5]); //exon 0 is the complete gene region including all exons. real exon starts from exon 1.
		
		continue;
	}
	
	if ($arrFields[2] != "cds") {
		continue; //only look at "cds" features
	}
	
	$arrAllGenes[$sChr][count($arrAllGenes[$sChr])-1][] = array("Tstart" => $nTranslatedStart, "Tend" => $nTranslatedEnd, "Orient" => $sOrientation, "Phase" => $nPhase); 
	//add exon info
	
	
	//now write to output
	//JH556661.1	protein_coding	exon	195633	195699	.	-	.	 gene_id "ENSXMAG00000014710"; transcript_id "ENSXMAT00000014803"; exon_number "2"; gene_name "TCF7L2"; gene_biotype "protein_coding"; transcript_name "TCF7L2-201"; exon_id "ENSXMAE00000168230";
	/*
	$sOut = $sChr."\t";
	$sOut .= "protein_coding\texon\t";
	$sOut .= $nTranslatedStart."\t";
	$sOut .= $nTranslatedEnd."\t";
	$sOut .= ".\t";
	$sOut .= $sOrientation."\t.\t";
	$sGeneId = "MANUALG_".$sGeneFamily."_".$nGeneCount;
	$sOut .= "gene_id \"$sGeneId\"; transcript_id \"$sGeneId\"; exon_number \"$nExonCount \"; gene_name \"$sGeneId\"; gene_biotype \"protein_coding\"; transcript_name \"$sGeneId\"; exon_id \"MANUAL_".$sGeneFamily."_EX_".$nTotalExonCount."\";" .PHP_EOL;
	
	fwrite($hTranslatedGTF , $sOut);
	*/
}
//print_r($arrAllGenes);

/*
if ($bRemoveOverlap) {
echo("Picking best records on each chromosome...\n");
foreach($arrAllGenes as $sChr => $arrGenes) {
	echo("\tContig: $sChr \n");
	
	$arrGenesToDelete = array();
	for($i=0;$i<count($arrGenes);$i++) {
		for($j=$i+1;$j<count($arrGenes);$j++) {
			$oGene1 = $arrGenes[$i];
			$oGene2 = $arrGenes[$j];

			if ($oGene1[0]['Score'] > $oGene2[0]['Score']) {
				$arrGenesToDelete[] = $j; //delete $j
			} else {

				$arrGenesToDelete[] = $i; //delete $j
			}
			
			
			if (fnIsOverlap($oGene1[0]['Tstart'], $oGene1[0]['Tend'], $oGene2[0]['Tstart'], $oGene2[0]['Tend']) && $oGene1[0]['Orient'] == $oGene2[0]['Orient']) { // there's overlap
				//decide which one to keep
				$nGene1Len = $oGene1[0]['Tend'] - $oGene1[0]['Tstart'];
				$nGene2Len = $oGene2[0]['Tend'] - $oGene2[0]['Tstart'];
				if ($nGene1Len >= $nGene2Len) { //keep gene 1
					echo("Keep: $sChr:".$oGene1[0]['Tstart']."-".$oGene1[0]['Tend']." Remove: ".$oGene2[0]['Tstart']."-".$oGene2[0]['Tend'].PHP_EOL);
					$arrGenesToDelete[] = $j; //delete $j
				}
				else {
					echo("Keep: $sChr:".$oGene2[0]['Tstart']."-".$oGene2[0]['Tend']." Remove: ".$oGene1[0]['Tstart']."-".$oGene1[0]['Tend'].PHP_EOL);
					$arrGenesToDelete[] = $i;  //delete $i
				}
			}
			
		}
	}
	
	$arrGenesToDelete = array_unique($arrGenesToDelete);
	foreach($arrGenesToDelete as $nInxToDel) {
		unset($arrGenes[$nInxToDel]);
	}
	$arrAllGenes[$sChr] = $arrGenes;
}
}
*/
//print_r($arrAllGenes);
echo("Writing GTF ... \n");
$nGeneCount = 0;
$nTotalExonCount = 0;

foreach($arrAllGenes as $sChr => $arrGenes) {

	$nCurrGoodScore = -1;
	$nGeneInMatch = -1;
	$nExonCount = -1;	
	foreach($arrGenes as $oGene) {



		//check if region is found in the top scored genes:

		$bFoundInGoodHit = false;
		foreach( $arrTopScores as $nScore => $oGeneHit) {
			if ($oGeneHit['Contig'] == $sChr && $oGeneHit['Start'] <= $oGene[0]['Tstart'] && $oGeneHit['End'] >= $oGene[0]['Tend']) { //this match is part of a good hit
				$bFoundInGoodHit = true;
				if ($nCurrGoodScore != $nScore) {
					$nGeneInMatch = 1;
					$nCurrGoodScore = $nScore;
					$nGeneCount++;
					$nExonCount = 0;

					$sOut = $sChr."\t";
					$sOut .= "protein_coding\tgene\t";
					$sOut .= $oGeneHit['Start']."\t";
					$sOut .= $oGeneHit['End']."\t";
					$sOut .= $nScore."\t";
					$sOut .= $oGene[0]['Orient']."\t.\t";
					$sGeneId = "GENEWISE_".$sGeneFamily."_".$nGeneCount;
					$sOut .= "gene_id \"$sGeneId\"; transcript_id \"$sGeneId\"; gene_name \"$sGeneId\"; gene_biotype \"protein_coding\";" .PHP_EOL;
					fwrite($hTranslatedGTF , $sOut);


				} else {
					$nGeneInMatch++;
				}
				
				break;
			}
		}

		if (!$bFoundInGoodHit) {
			continue;
		}
 
		foreach($oGene as $nEx => $oExon) {
			
			if ($nEx == 0) {


				continue; // skip the record for the whole gene region.
			}
			$nExonCount++;	
			$nTotalExonCount++;
			
			$sOut = $sChr."\t";
			$sOut .= "protein_coding\texon\t";
			$sOut .= $oExon['Tstart']."\t";
			$sOut .= $oExon['Tend']."\t";
			$sOut .= ".\t";
			$sOut .= $oExon['Orient']."\t". $oExon['Phase']."\t";
			$sGeneId = "GENEWISE_".$sGeneFamily."_".$nGeneCount;
			$sOut .= "gene_id \"$sGeneId\"; transcript_id \"$sGeneId\"; exon_number \"$nExonCount\"; gene_name \"$sGeneId\"; gene_biotype \"protein_coding\"; transcript_name \"$sGeneId\"; exon_id \"GENEWISE_".$sGeneFamily."_EX_".$nTotalExonCount."\";" .PHP_EOL;
			fwrite($hTranslatedGTF , $sOut);
		}
		
	}

}

function fnParseSegmentCoord($sSegmentStr) { //group18:9215509..9226706
	$arr1 = explode(":", $sSegmentStr);
	$arr2 = explode("..", $arr1[1]);
	return array($arr1[0], $arr2[0], $arr2[1]);
}

function fnIsOverlap($nStart1, $nEnd1, $nStart2, $nEnd2) {
	if ($nStart1 > $nEnd1 || $nStart2 > $nEnd2) {
		die("fnIsOverlap error: $nStart1, $nEnd1, $nStart2, $nEnd2");
	}
	
	return !($nEnd1 < $nStart2 || $nEnd2 < $nStart1);
	
}


?>
