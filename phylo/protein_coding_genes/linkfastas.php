<?php
//link files
$nMinTaxa = 40; //minimal number of taxa to include the cds alignment
$nMinLnCount = $nMinTaxa*2;
$sSrcFile = "../hyphyrelax/codeml_46_spp/rerun_exclude_MNM/sum_Callopanchax.txt"; //this can be any summary results from a codeml scan, this file is just going to reuse the generated fas file
$sLkDir = "lnkfasta/";

exec("mkdir -p $sLkDir");

$h = fopen($sSrcFile, 'r');

while(false !== ($sLn = fgets($h) ) ) {
	$sLn = trim($sLn);
	$arrF = explode("\t", $sLn);
	if (count($arrF) < 8 ) continue;
	$sOrthoID = $arrF[0];
	$sDir = trim($arrF[7]);
	$sFas = "$sDir/in.fas";
	if (!file_exists($sFas)) {
		echo("Warning: $sOrthoID\t$sFas not found any more.\n");
	}

	$nLnCount = intval(exec("wc -l < $sFas "));
	if ($nLnCount < $nMinLnCount) continue;

	exec("ln -s  $sFas $sLkDir/$sOrthoID.fas");

	echo("$sOrthoID\t".($nLnCount/2)."\n");
}

?>
