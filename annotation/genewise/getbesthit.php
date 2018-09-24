<?php
ini_set('memory_limit',-1);

$arrOpts = getopt("i:o:c:x:s:f:e:");

$sIn = $arrOpts["i"];
$sOut = $arrOpts["o"];
$sStrandOut = $arrOpts["c"];
$nExtend = intval($arrOpts["x"]) ; //3kb
$nExtendEnds = intval($arrOpts["s"]) ; //3kb
$sQueryFasta = $arrOpts["f"] ; //the fasta file used for query
$nEValueCutoff = 1e-5; // this is now used differently. the cutoff applies for weighted evalues for each alignment block.

$arrRegionsByQuery = array(); //regions, first level key - query, second level key, chromosome, third level - hit start, value - array[hit end, query start, query end, hit length, evalue]

//get the length of each query fasta sequence:	
$hBlastRet = fopen($sIn, "r");
//fnLoadRegions($hBlastRet , $arrRegions);
fnLoadRegionsByQuery($hBlastRet , $arrRegionsByQuery);

//print_r($arrRegionsByQuery);

$arrTopBlocks = array();
foreach($arrRegionsByQuery as $sQueryName => $arrRegions) {

	//echo(" $sQueryName Positive strand matches");
	//print_r(fnFindHitBlocks($arrRegionsByQuery[$sQueryName], true)); //look on positive strand

	$arrTopBlocks = $arrTopBlocks + fnGetTopBlocks(fnFindHitBlocks($arrRegionsByQuery[$sQueryName], true) , fnFindHitBlocks($arrRegionsByQuery[$sQueryName], false), 4); //get the best 2 blocks 
	//echo(" $sQueryName Negative strand matches");
	//print_r(fnFindHitBlocks($arrRegionsByQuery[$sQueryName], false)); //look on negative strand

	//fnMergeOverlap($arrRegionsByQuery[$sQueryName]);
}

//echo("Selected:\n");
//print_r($arrTopBlocks);

//now reformat the regions for fnMergeOverlap
$arrRegions = array(); //regions, first level key - chromosome name, second level : Start => End, sorted by Start.
foreach($arrTopBlocks as $nScore => $arrBlockDat) {
	$sContig = $arrBlockDat["Contig"];
	if (!array_key_exists( $sContig , $arrRegions) ) {
		$arrRegions[$sContig] = array();
	}

	$nDatBlockStart = $arrBlockDat["Start"];
	$nDatBlockEnd  = $arrBlockDat["data"]["End"];
	$nDatBlockStrand = $arrBlockDat["data"]["Strand"];
	if (!array_key_exists( $nDatBlockStart , $arrRegions[$sContig]) ) {
		$arrRegions[$sContig][$nDatBlockStart] = array( $nDatBlockEnd , $nDatBlockStrand);
	} else {
		$arrRegions[$sContig][$nDatBlockStart] = $nDatBlockEnd > $arrRegions[$sContig][$nDatBlockStart][0] ? 
                                                               array( $nDatBlockEnd , $nDatBlockStrand) : array( $arrRegions[$sContig][$nDatBlockStart][0]  ,  $arrRegions[$sContig][$nDatBlockStart][1]);
	}
}

//echo("reformated:\n");
//print_r($arrRegions);

fnMergeOverlap($arrRegions);

//echo("Merged\n");
//print_r($arrRegions);

//print_r($arrRegions);

$hOut = fopen($sOut , "w");
$hStrandOut = fopen($sStrandOut , "w");
//fwrite($hOut , "Contig\tStart\tEnd".PHP_EOL);

foreach($arrRegions as $sContig => $arrCoord) {
	foreach($arrCoord as $nStart => $arrDat ) {
		list($nEnd, $bStrand) = $arrDat;
		$nStartNew = ($nStart - $nExtendEnds>0)? $nStart - $nExtendEnds : 1;
		$nEndNew = $nEnd + $nExtendEnds;
		fwrite($hOut , "$sContig:$nStartNew..$nEndNew".PHP_EOL);
		fwrite($hStrandOut , ($bStrand?'+':'-').PHP_EOL);
	}
}


function fnLoadRegionsByQuery($hBlastRet, &$arrRegionsByQuery) { 
	global $nEValueCutoff;
	while( false!==($sLn=fgets($hBlastRet)) ) {
		$sLn = trim($sLn);
		if (substr($sLn,0,1) == "#") {
			continue;
		}
		
		$arrFields = explode("	", $sLn);

		$sQueryName = trim($arrFields[0]);
		$sContig = trim($arrFields[1]);
		$nHitLen = intval($arrFields[3]);
		$nMisMatches = intval($arrFields[4]);
		$bStrand = ($arrFields[8] < $arrFields[9]); //true = positive strand, false = negative strand
		$nQStart = intval($arrFields[6]);
		$nQEnd = intval($arrFields[7]);
		$nHStart = intval( $bStrand? $arrFields[8]:$arrFields[9]);
		$nHEnd = intval( $bStrand? $arrFields[9]:$arrFields[8] );
		$nEvalue = $arrFields[10];
		if ($nEValueCutoff < $nEvalue) {
			continue;
		}

 //regions, first level key - query, second level key, chromosome, third level - hit start, value - array[hit end, query start, query end, hit length, evalue]
		if(!array_key_exists($sQueryName, $arrRegionsByQuery)) {
			$arrRegionsByQuery[$sQueryName] = array();
		}

		if(!array_key_exists($sContig, $arrRegionsByQuery[$sQueryName])) {
			$arrRegionsByQuery[$sQueryName][$sContig] = array();
		}
		
		if (!array_key_exists( $nHStart , $arrRegionsByQuery[$sQueryName][$sContig] ) ) {
			//preserve the longer end
			$arrRegionsByQuery[$sQueryName][$sContig][$nHStart] = array("HEnd" => $nHEnd, "QStart" => $nQStart, "QEnd"=>$nQEnd, "HitLen"=>$nHitLen, "EValue" =>$nEvalue, "Strand" => $bStrand , "Mismatches" => $nMisMatches);
		}
		
	}
	
	
	foreach($arrRegionsByQuery as $sQueryName => $arrRecord) {
		
		foreach($arrRecord as $sContig => $arrReg) {
			ksort($arrRegionsByQuery[$sQueryName][$sContig] , SORT_NUMERIC );
		}

		ksort($arrRegionsByQuery[$sQueryName] , SORT_STRING);

	}

	ksort($arrRegionsByQuery, SORT_STRING);
	

}



function fnMergeOverlap(&$arrRegions) {

	foreach($arrRegions as $sContig => $arrCoord) {
		
		$nRegions = -1;
		//echo("======== $sContig ==========\n");
		while(true) {
			$arrRegions[$sContig] = fnMergeOverlapContig($arrRegions[$sContig]);
			//print_r($arrRegions[$sContig]);
			
			if ($nRegions == count($arrRegions[$sContig])) {
				break;
			}
			else {
				$nRegions = count($arrRegions[$sContig]);
			}
		}
		
	}
	

}

function fnFindHitBlocks($arrRegions, $bLookAtStrand=true) { //which strand to look at 

	$arrRetRegions = array();
	foreach($arrRegions as $sContig => $arrCoord) {
		//echo("Doing $sContig\n");
		//echo("arrCoord\n");
		//print_r($arrCoord);
		if ($bLookAtStrand) {
			ksort($arrCoord , SORT_NUMERIC);
		} else {
			krsort($arrCoord , SORT_NUMERIC);
		}
		reset($arrCoord);
		//print_r(current($arrCoord));
		$arrBlocksOnContig = array(); //blocks on contig, first level key: start, value array("End" , "Weighted Evalue");
		$nCurrStrand = true;
		$nCurrBlockStart = -1;
		$nCurrBlockEnd = -1;
		$nPrevQStart = -1;
		$nCurrBlockHitLen = -1; //accumulative hit length
		$nAccEValue = 0; //accumulative e value, multiplied by each hit block's length
		$nAccMismatches = 0;
		while($arrCurrHit = current($arrCoord)) {

			$nCurrHitStart = key($arrCoord);
			//echo("key: " . key($arrCoord) . PHP_EOL);
			$nCurrStrand = $arrCurrHit["Strand"];
			next($arrCoord);

			if ( $nCurrStrand != $bLookAtStrand) {
				continue; //only look at one strand.
			}
		
			$nCurrQStart = $nCurrStrand? $arrCurrHit["QStart"] : $arrCurrHit["QEnd"] ;
			//see if the query start is monotonous
			if ( ( $nCurrStrand?  $nCurrQStart <= $nPrevQStart : $nCurrQStart >= $nPrevQStart ) || $nCurrBlockStart == -1) { // not monotonous, break detected.

				if ($nCurrBlockStart !== -1) {
					$arrBlocksOnContig[$nCurrBlockStart] = array("End" => $nCurrBlockEnd , "WEvalue" => $nAccEValue / $nCurrBlockHitLen, "PropMismatch" => 1e14 * $nAccMismatches / $nCurrBlockHitLen , "HitLen" => $nCurrBlockHitLen, "Strand" => $bLookAtStrand);
				}

				$nCurrBlockStart = $nCurrHitStart;
				$nCurrBlockEnd = $arrCurrHit["HEnd"];
				$nPrevQStart = $nCurrStrand? $arrCurrHit["QStart"] : $arrCurrHit["QEnd"] ;
				$nCurrBlockHitLen =  $arrCurrHit["HitLen"];
				$nAccEValue = $arrCurrHit["EValue"] * $arrCurrHit["HitLen"];
				$nAccMismatches = $arrCurrHit["Mismatches"] ;
				continue;
			} else { // monotonous
				$nCurrBlockEnd = $arrCurrHit["HEnd"];
				$nPrevQStart = $nCurrStrand? $arrCurrHit["QStart"] : $arrCurrHit["QEnd"] ;
				$nCurrBlockHitLen  +=  $arrCurrHit["HitLen"];
				$nAccEValue += ($arrCurrHit["EValue"] * $arrCurrHit["HitLen"]);
				$nAccMismatches += $arrCurrHit["Mismatches"] ;
			}
		}

		if ($nCurrBlockStart !== -1) {
			$arrBlocksOnContig[$nCurrBlockStart] = array("End" => $nCurrBlockEnd , "WEvalue" => $nAccEValue / $nCurrBlockHitLen, "PropMismatch" => 1e14 * $nAccMismatches / $nCurrBlockHitLen , "HitLen" => $nCurrBlockHitLen, "Strand" => $bLookAtStrand);
		}

		if (count($arrBlocksOnContig) > 0 ) {
			$arrRetRegions[$sContig] = $arrBlocksOnContig;
		}
		
	}

	//print_r($arrRetRegions);
	return $arrRetRegions;
}

function fnGetTopBlocks($arrPosStrandBlocks, $arrNegStrandBlocks, $nTopCount) {

	$arrTopBlocks = array(); //key is score, value is the block itself
	fnGetTopBlockSub($arrPosStrandBlocks,  $arrTopBlocks, $nTopCount);
	fnGetTopBlockSub($arrNegStrandBlocks,  $arrTopBlocks, $nTopCount);
	return $arrTopBlocks;
}

function fnGetTopBlockSub($arrBlocks,  &$arrCurrTops, $nTopCount) { //subfunction.

	foreach($arrBlocks as $sContig => $arrBlockOnContig) {
		foreach($arrBlockOnContig as $nBlockStart => $arrBlockData) {
			if (count($arrCurrTops) < $nTopCount) { //if top blocks haven't been filled 
				$arrCurrTops[$arrBlockData["PropMismatch"]] = array("Contig" => $sContig , "Start" => $nBlockStart, "data" => $arrBlockData);
			} else {

				$nDeleteThis = -1;
				foreach($arrCurrTops as $nScore => $arrBlock) {
					if ($nScore > $arrBlockData["PropMismatch"]) { // more mismatches than current, replace it with current
						$nDeleteThis = $nScore;
						break;
					}
				}

				if ($nDeleteThis !==-1) {
					unset($arrCurrTops[$nDeleteThis]);
					$arrCurrTops[$arrBlockData["PropMismatch"]] = array("Contig" => $sContig , "Start" => $nBlockStart, "data" => $arrBlockData);
					ksort($arrCurrTops, SORT_NUMERIC);
				}
			}
		}
	}
}

function fnMergeOverlapContig(&$arrCoord) {
		global $nExtend;
		$arrMergedCoord = array();

		reset($arrCoord);
		while(list($nCurrEnd , $bCurrStrand) = current($arrCoord)) {
			$nCurrStart = key($arrCoord);
			list($nNextEnd, $bNextStrand) = next($arrCoord);
			$nNextStart = key($arrCoord);
			
			
			$nMergedStart = $nCurrStart;
			
			if ($nNextEnd === false ) {
				$nMergedEnd = $nCurrEnd;
				
			}
			else {
			
				if ( ($nCurrEnd + $nExtend < $nNextStart) || ($bCurrStrand != $bNextStrand) ) { // current end before next start
					$nMergedEnd = $nCurrEnd; // keep this
					
				}
				else  {
					$nMergedEnd = max($nNextEnd ,  $nCurrEnd);
					next($arrCoord);
				}
				
			}
			
			$arrMergedCoord[$nMergedStart] = array( $nMergedEnd , $bCurrStrand);

			
		}
		
		return $arrMergedCoord;
}

?>
