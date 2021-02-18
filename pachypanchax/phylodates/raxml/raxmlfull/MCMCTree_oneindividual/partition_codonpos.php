<?php
//convert the input phylip file into MCMCTree partition format, where codon positions 1, 2 and 3 are in separated blocks.

$sIn = "one_ind_per_sp.phy";
$sOut = "one_ind_for_mcmctree.phy";

$h = fopen($sIn , 'r');
$hOut = fopen($sOut , 'w');

$nLn = 0;
$nLen = -1;

$arrPartitions = array();

while(false !== ($sLn = fgets($h) ) ) {
	$nLn++;
	if ($nLn == 1) continue;
	$sLn = trim($sLn);
	if ($sLn == '') continue;

	$arrRet = array("", "", "");
	list($sTaxon, $sSeq) = preg_split('/\s+/', $sLn);
	if ($nLen == -1) {
		$nLen = strlen($sSeq);
		if ($nLen % 3 != 0) {
			die("Error: input sequence not a multiplication of 3\n");
		}
	}

	if ($nLen != strlen($sSeq)) {
		die("Taxon $sTaxon has an abnormal sequence length!\n");
	}

	for($nPos=0;$nPos < $nLen; $nPos++) {
		$arrRet[$nPos % 3] .= $sSeq[$nPos];
	}

	$arrPartitions[$sTaxon] = $arrRet;
	
}

$nTaxonCount = count($arrPartitions);
$nPerPartitionLen = $nLen / 3;
for($nCodonPos=0; $nCodonPos<3; $nCodonPos++) {
	fwrite($hOut, " $nTaxonCount $nPerPartitionLen\n");
	foreach($arrPartitions as $sTaxon => $arrSeqs) {
		fwrite($hOut, "$sTaxon     ");
		fwrite($hOut, $arrSeqs[$nCodonPos]."\n");
	}

	fwrite($hOut, "\n");
}

?>
