<?php
$sBUSCOList = "full_table_v2NOR_v1.0.tsv"; #this is the output full table from BUSCO
$sGFFDir = "augustus_output/gffs"; #where are the gff outputs from busco?
$sOut = "singlecopy_coord.bed"; #the output coordinates of the single copy buscos

$hOut = fopen($sOut, "w");

$hBUSCOList = fopen($sBUSCOList , 'r');
$arrSingleCopyBUSCO = array();
while( false!== ($sLn = fgets($hBUSCOList)) ) {
	$sLn = trim($sLn);

	if ($sLn == "") continue;
	if ($sLn[0] == "#") continue;
	$arrF = explode("\t", $sLn);
	if ($arrF[1] != "Complete") continue;
	$arrSingleCopyBUSCO[] = $arrF[0];
}

$arrCoord = array();
foreach($arrSingleCopyBUSCO as $sSingleCopyBUSCO) {
	$sGFF = $sGFFDir."/".$sSingleCopyBUSCO.".gff";
	if (!file_exists($sGFF) ) {
		echo(" warning : $sGFF not found!\n");
		continue;
	}
	$hGFF = fopen($sGFF  , "r");
	while( false!== ($sLn = fgets($hGFF)) ) {
		$sLn = trim($sLn);
		if ($sLn == "") continue;
		if ($sLn[0] == "#") continue;
		$arrF = explode("\t", $sLn);
		if ($arrF[2] !== "CDS") continue;
		if (!array_key_exists($arrF[0], $arrCoord)) {
			$arrCoord[$arrF[0]] = array();
		}

		$arrCoord[$arrF[0]][$arrF[3]] = $arrF[4];
	}
}

ksort($arrCoord , SORT_STRING);
//sort 
foreach($arrCoord as $sScfld => $arrB) {
	ksort($arrCoord[$sScfld] , SORT_NUMERIC);
	foreach($arrCoord[$sScfld] as $nStart => $nEnd) {
		fwrite($hOut , "$sScfld\t".($nStart-1)."\t$nEnd\n");
	}
}

?>
