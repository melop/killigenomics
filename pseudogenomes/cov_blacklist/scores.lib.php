<?php
$nMinBlastPIdentity = 50; //filter HSP if lower than this.
$nMaxBlastPIdentityDiff = 20; //If identity of an HSP block is 20% smaller than other blocks in the same alignment, then exclude this hsp
$nTreatAsBigGap = 2; // if larger than this, treat as "big gap".
$nTreatAsHugeGap = 10; // if larger than this, treat as "huge gap".

function fnExec($sCmd, &$arrRet = false, &$retVal=false) {
	echo("\t\t\t > ".$sCmd.PHP_EOL);
	if ( $arrRet !== false ) {
		if ($retVal===false) {
			exec($sCmd , $arrRet);
		} else {
			exec($sCmd , $arrRet, $retVal);
		}
	} else {
		exec($sCmd);
	}
}

function fnScoreExonerateModel($sPredictedProt , $sEvidenceProt , $sDIR, $bIgnoreInternalStop = false ) {
	global $nMinBlastPIdentity , $nMaxBlastPIdentityDiff , $nTreatAsBigGap , $nTreatAsHugeGap; 

	//echo("$sPredictedProt\n\n$sEvidenceProt\n");
	$arrScores = array('score' => 0, 'internalstop' => false);

	if (substr($sPredictedProt , -1,1) == '*' ) {
		$sPredictedProt = substr($sPredictedProt , 0, strlen($sPredictedProt)-1 );
	}

	if (substr($sEvidenceProt , -1,1) == '*' ) {
		$sEvidenceProt = substr($sEvidenceProt , 0, strlen($sEvidenceProt)-1 );
	}

	if ($bIgnoreInternalStop) {
		$sPredictedProt = str_replace('*' , 'X', $sPredictedProt);
		if (strpos($sPredictedProt, '*') !== false ) {
			$arrScores['internalstop'] = true;
		}
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

function fnArr2Annot($arr) {
	$s = "";
	
	foreach($arr as $sKey => $sVal) {
		if ($sVal === false) $sVal=0;
		$s .= "$sKey=$sVal;";
	}

	return substr($s, 0, strlen($s)-1);
}

function fnParseAnnotation($s) {
		$arrMap = array();
		$arrF = explode(';' , $s);
		foreach($arrF as $v) {
			$arrPair = explode('=' , $v);
			if (count($arrPair) !=2) continue;
			$sKey = trim($arrPair[0]);
			$sValue = preg_replace('/^"|"$/', '', trim($arrPair[1]));
			$arrMap[$sKey] = $sValue;
		}
		return $arrMap;
}

?>
