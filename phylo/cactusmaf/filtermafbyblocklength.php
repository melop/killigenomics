<?php

$arrOpts = getopt("i:o:l:");

$sIn = $arrOpts['i'];
$sOut = $arrOpts['o'];

$hIn = fopen($sIn , "r");
$hOut = fopen($sOut , "w");

$nMinLen = $arrOpts['l']; //min block is 15

$arrAln = array();
$bProcess = false;
do {

	$sLn = fgets($hIn);

	if ($sLn===false || substr($sLn , 0, 1) =='a' ) {
		if ($bProcess) 	{
			fnProcess($arrAln);
		}
		$arrAln = array();
		$bProcess = false;
		if ($sLn===false) break;
	} else if (substr($sLn , 0, 1) =='s') {
		$arrFields = explode("\t", $sLn);
		if ($arrFields[3] >= $nMinLen) {
			$bProcess = true;
			$arrAln[] = $sLn;
		}
	} else if (trim($sLn)!="" ) {
		fwrite($hOut , $sLn);
	}

} while (true);

function fnProcess($arrAln) {
	global $hOut;
	if (count($arrAln) <=1 ) {
		return;
	}
	fwrite($hOut , "\na\tscore=1000\n");
	foreach($arrAln as $sLn) {
		$arrFields = explode("\t" , $sLn);
		$nLast = count($arrFields) - 1;
		$arrFields[$nLast] = str_replace('N' , 'n', $arrFields[$nLast]);
		fwrite($hOut , implode("\t",$arrFields));
	}
}

?>
