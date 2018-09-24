<?php
ini_set('memory_limit',-1);

$arrOpts = getopt("i:o:x:s:e:");

$sIn = $arrOpts["i"];
$sOut = $arrOpts["o"];
$nExtend = intval($arrOpts["x"]) ; //3kb
$nExtendEnds = intval($arrOpts["s"]) ; //3kb
$nEValueCutoff = 1e-20;
$arrRegions = array(); //regions, first level key - chromosome name, second level : Start => End, sorted by Start.


		
$hBlastRet = fopen($sIn, "r");
fnLoadRegions($hBlastRet , $arrRegions);

//print_r($arrRegions);

fnMergeOverlap($arrRegions);

//print_r($arrRegions);

$hOut = fopen($sOut , "w");
//fwrite($hOut , "Contig\tStart\tEnd".PHP_EOL);
foreach($arrRegions as $sContig => $arrCoord) {
	foreach($arrCoord as $nStart => $nEnd) {
		$nStartNew = ($nStart - $nExtendEnds>0)? $nStart - $nExtendEnds : 1;
		$nEndNew = $nEnd + $nExtendEnds;
		fwrite($hOut , "$sContig:$nStartNew..$nEndNew".PHP_EOL);
	}
}

function fnLoadRegions($hBlastRet, &$arrRegions) {
	global $nEValueCutoff;
	while( false!==($sLn=fgets($hBlastRet)) ) {
		$sLn = trim($sLn);
		if (substr($sLn,0,1) == "#") {
			continue;
		}
		
		$arrFields = explode("	", $sLn);
		if ($arrFields[10] > $nEValueCutoff ) {
			continue;
		}
		$sContig = trim($arrFields[1]);
		$nStart = intval( ($arrFields[8] < $arrFields[9])? $arrFields[8]:$arrFields[9]);
		$nEnd = intval(($arrFields[8] < $arrFields[9])? $arrFields[9]:$arrFields[8] );
		if(!array_key_exists($sContig, $arrRegions)) {
			$arrRegions[$sContig] = array();
		}
		
		if (array_key_exists( $nStart , $arrRegions[$sContig]) ) {
			//preserve the longer end
			$arrRegions[$sContig][$nStart] = max($nEnd , $arrRegions[$sContig][$nStart]);
		}
		else {
			$arrRegions[$sContig][$nStart] = $nEnd;
		}
		
	}
	
	foreach($arrRegions as $sContig => $arrCoord) {
	
		ksort($arrRegions[$sContig] , SORT_NUMERIC );
	}
	
	ksort($arrRegions , SORT_STRING );
	
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