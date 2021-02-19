<?php
$sIn = "fortreemix.txt.gz";
$arrUsePops = array('B', 'C', 'A', 'F', 'D'); // pool alleles from these pops together for calculations.

$h = popen("zcat -f $sIn", 'r');

$sHeader = trim(fgets($h));
$arrPopIDs = explode("\t", $sHeader);
$arrPopIDs = array_flip($arrPopIDs);
$arrUseCols = array_values(array_intersect_key($arrPopIDs , array_flip($arrUsePops) ));

$nPos = 0;
$arrDevCounts = array(
					0 => array(),
					1 => array(),
					2 => array()
				);
				
$nTotalAllele = 0;
while(false !== ($sLn = fgets($h) ) ) {
	$sLn = trim($sLn);
	$nPos++;
	$nCodonPos = $nPos % 3;
	$arrF = explode("\t", $sLn);
	$nTotalAnc = 0;
	$nTotalDev = 0;
	foreach($arrUseCols as $nCol) {
		$sF = $arrF[$nCol];
		list($nAncAlleleCount, $nDevAlleleCount) = explode(',', $sF);
		$nAncAlleleCount = intval($nAncAlleleCount);
		$nDevAlleleCount = intval($nDevAlleleCount);
		$nTotalAnc += $nAncAlleleCount;
		$nTotalDev += $nDevAlleleCount;
	}
	
	if (!array_key_exists($nTotalDev, $arrDevCounts[$nCodonPos])) {
		$arrDevCounts[$nCodonPos][$nTotalDev] = 0;
	}
	
	$arrDevCounts[$nCodonPos][$nTotalDev] += 1;
	$nTotalAllele = max($nTotalAllele,$nTotalAnc+$nTotalDev );
}

echo("Max allele count : $nTotalAllele\n");
$dN = $arrDevCounts[1][$nTotalAllele] + $arrDevCounts[2][$nTotalAllele];
$dS = $arrDevCounts[0][$nTotalAllele];

echo("dN : $dN, dS : $dS\n");
echo("Derive_al_count\tpN\tpS\n");
for($nDevACount=0; $nDevACount<$nTotalAllele; $nDevACount++) {
	$pN = fnExist0( $nDevACount, $arrDevCounts[1]) + fnExist0( $nDevACount, $arrDevCounts[2]);
	$pS = fnExist0( $nDevACount, $arrDevCounts[0]);
	echo("$nDevACount\t$pN\t$pS\n");
}

function fnExist0($sKey, $arr) {
	return (array_key_exists($sKey, $arr)? $arr[$sKey]:0 );
}

?>