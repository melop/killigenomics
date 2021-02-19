<?php
//make a map of rna id rnaX to XM_XXXX
//because in the UPhO analysis the rnaX ids were used.

$sGFF = "/beegfs_old/group_dv/home/RCui/notfurgenome/Jena/ref_Nfu_20140520_top_level.gff3";
$sOut = "NFZgeneIDMap.txt";

$h = fopen($sGFF, 'r');
$hO = fopen($sOut, 'w');

while(false !== ($sLn = fgets($h) ) ) {
	$sLn = trim($sLn);
	$arrF = explode("\t", $sLn);
	if (count($arrF) < 9) continue;
	if ($arrF[2] != 'mRNA') continue;
	$arrAnnot = fnStr2Array($arrF[8]);
	if (array_key_exists('ID', $arrAnnot) && array_key_exists('Name', $arrAnnot) ) {
		$sID = $arrAnnot['ID'];
		list($sName, $s) = explode('.', $arrAnnot['Name']);
		fwrite($hO, "$sID\t$sName\n");
	}
}

function fnStr2Array($s) {
	$arrRet = array();
	$arrF1 = explode(';', $s);
	foreach($arrF1 as $s2) {
		$arrF2 = explode('=', $s2);
		if (count($arrF2) == 2) {
			$arrRet[$arrF2[0]] = $arrF2[1];
		}
	}

	return $arrRet;
}
?>
