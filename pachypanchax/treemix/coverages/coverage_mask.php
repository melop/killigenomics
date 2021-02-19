<?php
//prepare a coverage mask file, read in bcf files, and find the peak coverage. exclude regions with abnormal sequence coverage.

$sBCFTOOLS = "/beegfs/group_dv/software/source/bcftools-1.6/bcftools";
$sSampleDef = "samples.txt";
$nSampleLines = 0.2e9;
$nTotalParts = intval($argv[1]);
$nThisPart = intval($argv[2]);

$arrSamples = fnLoadSamples($sSampleDef);

print_r($arrSamples);

$nRunCount = -1;
foreach($arrSamples['ref'] as $sPop => $arrSamples) {
	$nRunCount++;
	if ($nRunCount % $nTotalParts != $nThisPart) continue;
	fnProcess("ref_$sPop", $arrSamples);
}

/*
foreach($arrSamples['samples'] as $sPop => $arrSamples) {
	foreach($arrSamples as $sSampleID => $sBCF) {
		fnProcess("ref_$sPop"."_$sSampleID", array($sBCF) );
	}
}
*/

function fnProcess($sLabel, $arrSamples) {
	global $sBCFTOOLS, $nSampleLines;
	$sCmd = '';
	if (count($arrSamples) > 1) {
		$sCmd .= "$sBCFTOOLS merge ";
	} else {
		$sCmd .= "$sBCFTOOLS view ";
	}

	$sCmd .= implode(' ', $arrSamples);
	//$sCmd .= " | head -n " . ($nSampleLines + 10000);
	echo($sCmd."\n");

	$descriptorspec = array(
	   0 => array("pipe", "r"),  // stdin is a pipe that the child will read from
	   1 => array("pipe", "w"),  // stdout is a pipe that the child will write to
	   2 => array("pipe", "w") // stderr is a file to write to
	);

	$hProc = proc_open($sCmd, $descriptorspec, $pipes);


	$arrDepthCount = array();
	$nLnCount = 0;
	while( false !== ($sLn = fgets($pipes[1]) ) ) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		if ($sLn[0] == '#') continue;
		$arrF = explode("\t", $sLn);
		if (count($arrF) < 8) continue;

		preg_match('/DP=([0-9]+)/', $arrF[7], $arrM);
		if (count($arrM) != 2) continue;
		$nDP = intval($arrM[1]);
		if (!array_key_exists($nDP , $arrDepthCount)) {
			$arrDepthCount[$nDP] = 0;
		}

		$arrDepthCount[$nDP]+=1;

		//echo($nLnCount."\n");
		if ($nLnCount % 1e6 == 0 ) {
			echo("Proccessed ".($nLnCount / 1e6). "M positions\n");
		}

		if ( $nSampleLines < $nLnCount++) {
			//print_r($arrDepthCount);
			break;
		}
	}

	echo("finished.\n");
	$s = proc_get_status($hProc);
	posix_kill($s['pid'], SIGKILL);
	proc_close($hProc);


	echo("Sort\n");
	ksort($arrDepthCount, SORT_NUMERIC);
	echo("End Sort\n");
	//print_r($arrDepthCount);

	//die("Stop.\n");
	$hOut = fopen("cov_$sLabel", 'w');
	fwrite($hOut, "Depth\tCount\n");
	foreach($arrDepthCount as $nDP => $nCount) {
		fwrite($hOut, "$nDP\t$nCount\n");
	}
}

function fnLoadSamples($sSampleDef) {
	$h = fopen($sSampleDef, 'r');

	$arrRet = array('ref' => array(), 'samples' => array());
	while(false !== ($sLn =fgets($h))) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		if ($sLn[0] == '#') continue;

		list($sRef, $sPop, $sSampleID, $sBCF) = explode("\t", $sLn );

		$sType = ($sRef == 'Ref')? 'ref':'samples'; 
		if (!array_key_exists($sPop, $arrRet[$sType]) ) {
			$arrRet[$sType][$sPop] = array();
		}

		$arrRet[$sType][$sPop][] = $sBCF;		
	}

	return $arrRet;
}

?>
