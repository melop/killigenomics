<?php
require_once(dirname(__FILE__) . "/lib.php");

//get the clusters with exactly 16 species

$sIn = "out.seq.mci.I20";#$argv[1];
$sOutDir = "FullTaxaOrthologs";
$nTaxonCount = 16;
$hIn = fopen($sIn ,"r");

exec("mkdir -p $sOutDir");

$nAllTaxaIncluded = 0;
$nExactOrthologs = 0;
$nSingletons = 0;
$nClusters = 0;
$arrGoodGroups = array();

while(false !== ($sLn = fgets($hIn) ) ) {
	$sLn = trim($sLn);
	if ($sLn == '') continue;
	$arrF = explode("\t" , $sLn);

	$arrTaxa = array();
	$nNumSeq = 0;

	foreach($arrF as $sID) {
		$arrF2 = explode('|' , $sID);
		$arrTaxa[$arrF2[0]] = true;
		$nNumSeq++;
	}

	if (count($arrTaxa) ==  $nTaxonCount) {
		$nAllTaxaIncluded++;
	}

	if (count($arrTaxa) ==  $nTaxonCount && $nNumSeq==$nTaxonCount ) {
		$nExactOrthologs++;
		$arrGoodGroups[] = $arrF;
	}

	if (count($arrTaxa) == 1 && $nNumSeq==1) {
		$nSingletons++;
	} else {
		$nClusters++;
	}
}

echo("$sIn Clusters (non-singletons): $nClusters AllTaxaRepresented: $nAllTaxaIncluded ExactOrthoGroup: $nExactOrthologs  Singletons: $nSingletons\n");
//print_r($arrGoodGroups[0]);

$sTmpList = "$sOutDir/tmp.txt";
$hTmpList = fopen( $sTmpList , 'w');
foreach($arrGoodGroups as $arrGroup) {

	fwrite($hTmpList , implode(',' , $arrGroup) . "\n");

}

//execute UPhO script to get alignment:
exec("Get_fasta_from_Ref.py -o $sOutDir -q $sTmpList -R ../allproteins.fa"); 
//do alignment:
exec("cd $sOutDir; bash ../paMATRAX_ray.sh -t -c ");
?>
