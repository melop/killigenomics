<?php
$sGFF = "genes.gff";
$sIn = "readcounts.txt";
$sOut = "readcounts_spgeneid.txt";


$hGFF = fopen($sGFF, 'r');

$arrID2SpGeneId = array();

while( false !== ($sLn = fgets($hGFF)) ) {
	$sLn = trim($sLn);
	if ($sLn == '') continue;
	$arrF = explode("\t" , $sLn);
	if ($arrF[2] != 'mRNA') {
		continue;
	}

	preg_match("/spgeneid=([^;]+)/", $arrF[8], $arrM);
	if (count($arrM) != 2) {
		echo("Warning: spgeneid not found for line\n$sLn\n");
		continue;
	}

	$nSpGeneId = $arrM[1];

	preg_match("/ID=([^;]+)/", $arrF[8], $arrM);
	if (count($arrM) != 2) {
		echo("Warning: ID not found for line\n$sLn\n");
		continue;
	}

	$sGeneId = $arrM[1];

	$arrID2SpGeneId[$sGeneId] = $nSpGeneId;

}

$hIn = fopen($sIn, 'r');
$hOut = fopen($sOut, 'w');

while( false !== ($sLn = fgets($hIn)) ) {
	$sLn = trim($sLn);
	if ($sLn == '') continue;
	if ($sLn[0] == '#') {
		fwrite($hOut, $sLn ."\n");
		continue;
	}

	$arrF = explode("\t" , $sLn);
	if (array_key_exists($arrF[0] , $arrID2SpGeneId)) {
		$arrF[0] = $arrID2SpGeneId[$arrF[0]];
	} else {
		echo("Warning: $arrF[0] cannot be mapped to spgeneid\n");
	}

	fwrite($hOut, implode("\t", $arrF) ."\n");

}

?>
