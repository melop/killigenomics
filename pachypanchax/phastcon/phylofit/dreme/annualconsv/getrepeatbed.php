<?php

$sRMOut = "../../../../cactus/NFZ2.0/NFZ2.0.repeatmask.out";
$sOut = "nfz.repeats.bed";
$hO = fopen($sOut, 'w');

$arrMasks = fnReadMask($sRMOut);
ksort($arrMasks);

foreach($arrMasks as $sChr => $arrCoord) {
	ksort($arrMasks[$sChr], SORT_NUMERIC);
	foreach($arrCoord as $nStart => $nEnd) {
		fwrite($hO , "$sChr\t".($nStart-1)."\t".$nEnd."\n");
	}
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
