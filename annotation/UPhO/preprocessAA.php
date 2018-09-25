<?php
//this script reads in deduplicated peptide fasta, replace the "." in the ensembl ID with underscore _, and remove all content after the first space. a tab file is output to stderr.

$hIn = fopen("php://stdin", "r");
$hFastOut = fopen("php://stdout", "w");
$hTabOut = fopen("php://stderr", "w");

$arrGeneList = array();

$sSeqName = "";
$sSeq = "";
do {
        $sLn = fgets($hIn);
        if ($sLn === false) {
                fnCheck($sSeqName , $sSeq);
                break;
        }

        $sLn = trim($sLn);

        if ($sLn == "") continue;

        if ($sLn[0] == ">") {
                fnCheck($sSeqName , $sSeq);
                $sSeqName = substr($sLn, 1);
                $sSeq = "";
                continue;
        }

        $sSeq .= $sLn;             
} while(true) ;

fwrite($hTabOut , "GeneID\tLen\tTotalIsoforms\tTranscriptCount\tTranscriptID\tIncluded\tFullheader\n" );
foreach($arrGeneList as $sGeneID => &$arrTranscripts) {
	krsort($arrTranscripts , SORT_NUMERIC);
	$nCount = 0;
	$nTotalIsoforms = count($arrTranscripts);
	foreach($arrTranscripts as $nLen => &$arrSeq) {
		$nCount++;
		$sIncluded = "no";
		if ($nCount == 1) {
			fwrite($hFastOut , ">".$arrSeq['ensid']."\n".wordwrap($arrSeq['seq'], 100, "\n" , true)."\n" );
			$sIncluded = "yes";
		}

		fwrite($hTabOut , "$sGeneID\t$nLen\t$nTotalIsoforms\t$nCount\t".$arrSeq['ensid']."\t$sIncluded\t".$arrSeq['fullheader']."\n" );
	}
}

function fnCheck($sSeqName , $sSeq) {
	global $arrGeneList;
	if ($sSeqName == "" || $sSeq == "") return;

	$arrF = explode(" ", $sSeqName);
	$sEnsID = $arrF[0];
	$sEnsID = str_replace(".", "_", $sEnsID);


	$sGeneID = $sEnsID;

	if (!array_key_exists($sGeneID , $arrGeneList) ) {
		$arrGeneList[$sGeneID] = array();
	}

	$arrGeneList[$sGeneID][strlen($sSeq)] = array("ensid" => $sEnsID , "fullheader" => $sSeqName , "seq" => $sSeq);

}

?>
