<?php
//during initial submission, NCBI identifies some adapter contamination which either needs to be excluded or masked.

//$sSp = "PLP";

$sSp = $argv[1]; //"PLP";

if (trim($sSp) == '') {
	die("Usage: xxx.php spname\n");
}

$sOut = "$sSp/softmasked.fa";
list($sRMOut) = glob("$sSp/*out*");
list($sFa) = glob("$sSp/*.fna*");

$arrMasks = fnReadMask($sRMOut);
//print_r($arrMasks);
//die();
$hIn = popen( "zcat -f $sFa" , 'r');
$hOut = fopen($sOut, 'w');
$sSeq = '';
$sSeqName = '';
$nLineWidth = -1;

while(true) {
	$sLn = fgets($hIn);
	if ($sLn === false || ($sLn !== '' && $sLn[0]== '>') ) {
		if ($sSeqName != '' && $sSeq != '') {
			fnProcessSeq($sSeqName, $sSeq);
		}

		if ($sLn === false) {
			break;
		}

		if ($sLn[0]== '>') {
			list($sSeqName) = preg_split('/\s+/', substr(trim($sLn) , 1));
			$sSeq = '';
		}

		continue;
	}	

	if ($sLn == '') continue;

	if ($nLineWidth == -1) {
		$nLineWidth = strlen(trim($sLn));
	}

	$sSeq .= strtoupper(trim($sLn));

}

function fnProcessSeq($sSeqName, $sSeq) {
	global $arrMasks, $hOut, $nLineWidth;

	$nBeforeLen = strlen($sSeq);

	echo("Masking $sSeqName...\n");
	if (array_key_exists($sSeqName, $arrMasks) ) {
		foreach( $arrMasks[$sSeqName] as $nStart => $nEnd ) {
			$nLen = abs($nEnd - $nStart) + 1;
			$sN = strtolower(substr($sSeq , $nStart, $nLen));
			$sSeq = substr_replace($sSeq , $sN, $nStart - 1, strlen($sN));
			//echo("Soft Masked $sSeqName:$nStart..$nEnd, len $nLen\n");
		}
	}


	$nAfterLen = strlen($sSeq);

	if ($nBeforeLen != $nAfterLen) {
		echo("$sSeqName, before $nBeforeLen, after $nAfterLen bp\n");
	}

	fwrite($hOut , ">$sSeqName\n".wordwrap($sSeq, $nLineWidth, "\n", true)."\n" );

}

function fnReadMask($sMask) {
	$h = popen("zcat -f $sMask", 'r');

	$arrMask = array();
	$nCount = 0;
	while(false !== ($sLn = fgets($h)) ) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;

		$arrF = preg_split('/\s+/', $sLn);

		$nStart = intval($arrF[5]);
		$nEnd = intval($arrF[6]);

		if ($nStart == 0 || $nEnd==0 ) {
			continue;
		}

		if ($nStart > $nEnd) {
			die("Error $nStart > $nEnd \n$sLn\n");
		}

		$nCount++;

		if (!array_key_exists($arrF[4], $arrMask)) {
			$arrMask[$arrF[4]] = array();
		}

		$arrMask[$arrF[4]][$nStart] = $nEnd;

		

	}

	echo("Loaded $nCount repeats\n");
	return $arrMask;
}
?>
