<?php
//this version considers all fragments identified for the gene alignment, instead of only taking the longest one.
$arrOpts = getopt("f:o:k:");


//$sGeneFamily = $arrOpts["g"]; // gene name"v2r";
//$sGeneWiseGTF = "out/Xmac/splitted/all.ret.fasta";
//$sTranslatedGTF = "sGeneFamily_Genewise_translated.gtf";
$sGeneWiseGTF = $arrOpts["f"]; //"out/medaka/splitted/all.ret.fasta";
$sTranslatedGTF =  $arrOpts["o"];//"sGeneFamily_Genewise_translated_medaka.gtf";
$nMaxAllowedExonInsertionPerc = 0.01; //no insertion or deletions in the gene model, otherwise alignment was spurious
$nMaxAllowedExonDeletionPerc = 0.01; //no insertion or deletions in the gene model, otherwise alignment was spurious
$bRemoveOverlap = true;

$hGeneWiseGTF = fopen($sGeneWiseGTF , "r");
$hTranslatedGTF = fopen($sTranslatedGTF , "w");
$nGeneCount = 0;
$nExonCount = 0;
$nTotalExonCount = 0;
$arrAllGenes = array(); //key is chromosome name, secondary array are genes one chromosome

$nKeepHowMany = intval($arrOpts["k"]);//1;
$arrTopScores = array(); //key-score, value - unused.

$sGeneFamily = "Gene".rand();

echo("Loading exonerate GTF file...\n");
$bInsideGFFBlock = false;

$sGeneID = "";
$arrExcluded = array(); //gene ids to exclude because they don't align well enough

while(false!==($sLn=fgets($hGeneWiseGTF))) {
	$sLn = trim($sLn);
	if ($sLn == "") {
		continue;
	}

	if (strpos($sLn , "START OF GFF DUMP") !== false) {
		$bInsideGFFBlock = true;
		continue;
	}

	if (strpos($sLn , "END OF GFF DUMP") !== false) {
		$bInsideGFFBlock = false;
		continue;
	}

	if (!$bInsideGFFBlock)  {
		continue;
	}

	if (substr($sLn, 0,1) == '#') {
		continue;
	}
	
	$arrFields = explode("\t", $sLn);

	
	list($sChr, $nStart, $nEnd) = fnParseSegmentCoord($arrFields[0]);
	$nStart = min($nStart , $nEnd);
	$nLocalStart = min($arrFields[3], $arrFields[4]);
	$nLocalEnd   = max($arrFields[3], $arrFields[4]);
	$nTranslatedStart = $nStart -1 + $nLocalStart;
	$nTranslatedEnd	  = $nStart -1 + $nLocalEnd;
	$nFeatureLen = abs($nTranslatedEnd - $nTranslatedStart ) + 1;
	$sOrientation = $arrFields[6];
	$nPhase = $arrFields[7];
	$sComments = count($arrFields)<9? "." : $arrFields[8];
	$nAlnScore = $arrFields[5];
	$nExonCount++;
	$nTotalExonCount++;

	
	if ($arrFields[2] == "gene") {
		$nGeneCount++;
		$nExonCount = 0;
		
		if (!array_key_exists($sChr, $arrAllGenes)) {
			$arrAllGenes[$sChr] = array();
		}

		preg_match("/sequence (\S+)/", $sComments, $arrMatches);

		if (count($arrMatches) !=2) {
			$sGeneID = "Gene".rand();
		} else {
			$sGeneID = $arrMatches[1]."_hit_".$nGeneCount;
		}
		
		$arrAllGenes[$sChr][] = array(); // this array contains each exon
		$arrAllGenes[$sChr][count($arrAllGenes[$sChr])-1][0] = array("Tstart" => $nTranslatedStart, "Tend" => $nTranslatedEnd, "Orient" => $sOrientation, "Score" => $nAlnScore, "CDSLen" => 0, "ExonLen" => 0, "ExonDel" => 0, "ExonInsert" => 0); //exon 0 is the complete gene region including all exons. real exon starts from exon 1.
		//check if is top score:

		$oGeneHit = array('Contig' => $sChr, 'Start' => $nTranslatedStart , 'End' => $nTranslatedEnd, 'Comment' => preg_replace("/gene_id[^;]+;/", "", $sComments), 'GeneID' => $sGeneID );
		if (count($arrTopScores) <  $nKeepHowMany) {
			$arrTopScores[$nAlnScore] = $oGeneHit ; // write down score.
		} else {	
	
			$nToUnSet = false;
			foreach($arrTopScores as $nScore => $nval) {
				if ($nScore <  $nAlnScore) {
					$nToUnSet = $nScore; 
					break;
				}
			}
	
			if ($nToUnSet !== false) {
				unset($arrTopScores[$nToUnSet]);
				$arrTopScores[ $nAlnScore] = $oGeneHit; // write down score.
				ksort($arrTopScores , SORT_NUMERIC);
			}
		}
		continue;
	}
	
	if ($arrFields[2] == "cds") {
		$arrAllGenes[$sChr][count($arrAllGenes[$sChr])-1][0]['CDSLen'] += $nFeatureLen;
		$arrFields[2] = "CDS";
	}

	if ($arrFields[2] == "exon") {
		$arrAllGenes[$sChr][count($arrAllGenes[$sChr])-1][0]['ExonLen'] += $nFeatureLen;
	}

	preg_match("/insertions\s+(\d+).+deletions\s(\d+)/", $sComments, $arrMatches);

	if (count($arrMatches) ==3) {
		if ($arrFields[2] == "exon") {
			$arrAllGenes[$sChr][count($arrAllGenes[$sChr])-1][0]['ExonInsert'] += $arrMatches[1];
			$arrAllGenes[$sChr][count($arrAllGenes[$sChr])-1][0]['ExonDel'] += $arrMatches[2];
		}
	}

	$arrAllGenes[$sChr][count($arrAllGenes[$sChr])-1][] = array("Tstart" => $nTranslatedStart, "Tend" => $nTranslatedEnd, "Orient" => $sOrientation, "Phase" => $nPhase, "Feature" => $arrFields[2],  "Comment" => $sComments, 'GeneID' => $sGeneID); 
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

		//check if there are too many insertions or deletions in this gene
		if ($oGene[0]['ExonInsert'] / $oGene[0]['ExonLen'] > $nMaxAllowedExonInsertionPerc || $oGene[0]['ExonDel'] / $oGene[0]['ExonLen'] > $nMaxAllowedExonDeletionPerc) {

			continue;
		}

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
					$sOut .= "exonerate:cdna2genome\tgene\t";
					$sOut .= $oGeneHit['Start']."\t";
					$sOut .= $oGeneHit['End']."\t";
					$sOut .= $nScore."\t";
					$sOut .= $oGene[0]['Orient']."\t.\t";
					$sOut .= "gene_id \"". $oGeneHit['GeneID'] . "\"; ". "transcript_id \"". $oGeneHit['GeneID'] . "\"; ";
					$sOut .= "cds_len " . $oGene[0]['CDSLen'] . "; exon_len " . $oGene[0]['ExonLen'] . "; exon_insertions " . $oGene[0]['ExonInsert'] . "; exon_deletions " . $oGene[0]['ExonDel'] . "; " . $oGeneHit['Comment'].PHP_EOL;
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
			$sOut .= "exonerate:cdna2genome\t". $oExon['Feature'] ."\t";
			$sOut .= $oExon['Tstart']."\t";
			$sOut .= $oExon['Tend']."\t";
			$sOut .= ".\t";
			$sOut .= $oExon['Orient']."\t". $oExon['Phase']."\t";
			$sOut .= "gene_id \"". $oGeneHit['GeneID'] . "\"; ". "transcript_id \"". $oGeneHit['GeneID'] . "\"; ". $oExon['Comment'].PHP_EOL; 
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
