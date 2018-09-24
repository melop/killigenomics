<?php
//This script processes the megablast output, then output a list of likely contaminated contigs
$sContigLen = "original_congitlen.sorted.txt"; //the line number correspond to the contig number in the blast result
$sMegaBlastRet = "contig.megablast.txt"; //"contig.megablast.debug.txt"; // 
$nContamCutoff = 99; //need to be higher than 99% to be considered contamination
$nAtLeastLen = 0.2; // at least 50% are identical 
$nExtend = 0 ; //merging tolerance between blocks
$hOutContam = fopen("contaminated.list.txt", "w");


$hContigLen = fopen($sContigLen , "r");

$arrContigLen = array();

$nContigNum = 1;
while(false !== ($sLn = fgets($hContigLen)) ) {
	list($sContig, $nLen) = explode("\t" , trim($sLn) );
	if ($nLen == 0) {die("Zero length contig? $sContig\n");}
	$arrContigLen[$nContigNum] = array($sContig , $nLen);
	$nContigNum++;
}

$hMegaBlastRet = fopen($sMegaBlastRet , "r");

$nCurrContigNum = -1;
$arrCurrContig = array(); // this array contains ranges

while(false !== ($sLn = fgets($hMegaBlastRet)) ) {
	$sLn = trim($sLn);
	if ($sLn[0] == "#" || $sLn =="") continue;

	$arrF = explode("\t" , $sLn );
	$nContigNum = $arrF[0];
	$nHitIdentity = $arrF[2];
	$nHitStart = $arrF[6];
	$nHitEnd = $arrF[7];

	if ($nHitStart > $nHitEnd) {die("file error: hit start larger than hit end:\n$sLn\n");}

	if ($nContigNum != $nCurrContigNum) {
		fnProcess($nCurrContigNum , $arrCurrContig);
		$nCurrContigNum = $nContigNum;
		$arrCurrContig = array();
	}

	if ($nHitIdentity < $nContamCutoff) continue; //not similar enough 

	if (array_key_exists($nHitStart , $arrCurrContig)) {
		if ($nHitEnd > $arrCurrContig[$nHitStart] ) {
			$arrCurrContig[$nHitStart] = $nHitEnd;
		}
	} else {
			$arrCurrContig[$nHitStart] = $nHitEnd;
	}


	
}

fnProcess($nCurrContigNum , $arrCurrContig);

function fnProcess($nCurrContigNum , $arrCurrContig) {
	global $arrContigLen, $nAtLeastLen, $hOutContam;
	if ($nCurrContigNum == -1 || count($arrCurrContig) == 0 ) {
		return;
	}

	if (!array_key_exists($nCurrContigNum  , $arrContigLen) ) {
		die("Contig $nCurrContigNum not found in length description.\n");
	}

	//print_r($arrCurrContig);
	ksort($arrCurrContig , SORT_NUMERIC);
	//print_r($arrCurrContig);

	$arrMergedBlocks = $arrCurrContig;

	$nRegions = -1;
		//echo("======== $sContig ==========\n");
	while(true) {
		$arrMergedBlocks = fnMergeOverlapContig($arrMergedBlocks);
		//print_r($arrRegions[$sContig]);
			
		if ($nRegions == count($arrMergedBlocks)) {
			break;
		}
		else {
			$nRegions = count($arrMergedBlocks);
		}
	}


	//print_r($arrMergedBlocks);

	$nTotalSimLen = 0;
	foreach($arrMergedBlocks as $nStart => $nEnd) {
		$nBlockLen = $nEnd - $nStart + 1;
		$nTotalSimLen += $nBlockLen;
	}

	$nPercentLen = $nTotalSimLen/$arrContigLen[$nCurrContigNum][1];
	echo($nCurrContigNum."\t".$arrContigLen[$nCurrContigNum][0]."\t".$arrContigLen[$nCurrContigNum][1]."\t$nPercentLen\n");
	if ( $nPercentLen >= $nAtLeastLen ) {
		fwrite($hOutContam , $nCurrContigNum."\t".$arrContigLen[$nCurrContigNum][0]."\t".$arrContigLen[$nCurrContigNum][1]."\t$nPercentLen\n");
	}
}

function fnMergeOverlapContig(&$arrCoord) {
		global $nExtend;
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


?>
