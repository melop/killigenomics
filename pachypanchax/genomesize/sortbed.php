<?php
$sScfOrder = $argv[1]; //"scf_order.txt"; // sort the bed file according to this order
$sInBed = $argv[2]; //"/beegfs/group_dv/home/RCui/killifish_genomes/denovo/discovardenovo/NOR/genomesize_mapping/singlecopy_coord.bed";
$sOutBed = $argv[3]; //"singlecopy_busco.sorted.bed";

$hScfOrder = fopen($sScfOrder , 'r');
$arrSorted = array();

while(false !== ($sLn = fgets($hScfOrder) ) ) {
	$sLn = trim($sLn);

	if ($sLn == '') continue;

	list($sScfName, $nLen) = explode("\t" , $sLn);
	$arrSorted[$sScfName] = array();
	
}

$hIn = fopen($sInBed , 'r');
$hOut = fopen($sOutBed , 'w');

while(false !== ($sLn = fgets($hIn) ) ) {
	$sLn = trim($sLn);

	if ($sLn == '') continue;

	list($sScfName, $nStart, $nEnd) = explode("\t" , $sLn);
	if (!array_key_exists( $sScfName, $arrSorted)) {
		die("$sScfName 's order is not given!\n");
	}

	if (array_key_exists($nStart , $arrSorted[$sScfName]) ) {
		die("$sScfName : $nStart already defined previously!\n");
	}

	
	$arrSorted[$sScfName][$nStart] = $nEnd;
}

foreach($arrSorted as $sScfName => &$arrBedOnScf) {
	ksort($arrBedOnScf , SORT_NUMERIC);
	foreach($arrBedOnScf as $nStart => $nEnd) {
		fwrite($hOut , "$sScfName\t$nStart\t$nEnd\n");
	} 
}

?>
