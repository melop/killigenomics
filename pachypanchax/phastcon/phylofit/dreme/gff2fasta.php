<?php
$sGFF = $argv[1];
$sFas = $argv[2];

$nMinLen = 25;

$hG = fopen($sGFF, 'r');
$arrFas = fnParseFasta($sFas);
//echo(count($arrFas));
//die();

while( false !== ($sLn = fgets($hG) ) ) {
	$sLn = trim($sLn);

	$arrF = explode("\t", $sLn);

	$sChr = $arrF[0];
	$nStart = $arrF[3];
	$nEnd = $arrF[4];
	$nLen = $nEnd - $nStart + 1;
	if ($nLen <= 0) {
		die("Error coord. $sLn\n");
	}

	if ($nLen < $nMinLen) {
		continue;
	}

	$sID = $arrF[2]."|$sChr|$nStart|$nEnd";
	echo(">$sID\n");
	echo(substr($arrFas[$sChr], $nStart-1, $nLen));
	echo("\n");	
}

function fnParseFasta($sFas) {
	$arrRet = array();
	$sSeq = '';
	$sName = '';
	$h = fopen($sFas, 'r');
	while(true) {
		$sLn = fgets($h);
		if ($sLn === false || (trim($sLn)!='' && $sLn[0] == '>') ) {
			if ($sName != '') {
				if (array_key_exists($sName, $arrRet) ) {
					die("$sName is repeated in fasta\n");
				}
				$arrRet[$sName] = $sSeq;
				$sSeq = '';
				$sName = '';
			}

			if ($sLn === false) break;
			$sName = substr(trim($sLn),1);
			continue;
		}

		$sSeq .= trim($sLn);
	}

	return $arrRet;
}

?>
