<?php
$sIn = "../PLP.gphocs.1000.50000.in";
$sOutStem = "PLP.gphocs.1000.50000.";
$nPart = 4;

$h = fopen($sIn, 'r');
$nTotal = intval(fgets($h));
$nCountPerPart = intval($nTotal / $nPart);

$arrOut = array();
for($i = 0 ;$i<$nPart; $i++) {
	$arrOut[$i]= fopen($sOutStem.$i.".in" , 'w');
	fwrite($arrOut[$i], $nCountPerPart."\n\n");
}

$nCurr   = -1;
while( false !== ($sLn = fgets($h)) ) {
	$sLn = trim($sLn);
	if ($sLn == '') continue;
	if (count(explode(' ', $sLn) )==3) {
		if ($nCurr!=-1) fwrite($arrOut[$nCurr], "\n");
		$nCurr++;
		if ($nCurr >= $nPart) $nCurr =    0;
	}	

	fwrite($arrOut[$nCurr], $sLn . "\n");
}
?>
