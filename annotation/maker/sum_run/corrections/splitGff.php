<?php
//splits the gff file by types, and also converts "match" to "mRNA", "match_part" to "exon".
$sIn = $argv[1]; //"../genome_pilon.all.gff";
$sOutPrefix = "split_";

$hIn = fopen($sIn, "r");
$arrMainModels = array('maker' => true); //this is the main model
$arrCandidateModels = array('pred_gff:genewise' => true , 'evm' => true);
$arrProtein2GenomeModels = array('protein2genome:NCBIRefseqTeleostProt' => true , 'protein2genome:TeleosteiProt' => true, 'protein2genome:CyprinodontiformesProt' => true, 'protein2genome:AtherinomorphaeProt' => true);
$arrOutFiles = array();

while(false !== ($sLn = fgets($hIn) ) ) {

	$sLn = trim($sLn);
	if ($sLn == "") continue;

	if ($sLn == "##FASTA") break;
	if ($sLn[0] == '#') continue;

	$arrF = explode("\t" , $sLn, 4);

	if (count($arrF) != 4) continue;

	$sType = $arrF[1];
	$sFeature = fnRenameFeature($arrF[2]);

	$sOutName = fnType2FileName($sType);

	if ($sOutName === false) continue;

	if (!array_key_exists($sOutName , $arrOutFiles)) {
		$arrOutFiles[$sOutName] = fopen($sOutName, 'w');
	}

	fwrite($arrOutFiles[$sOutName] , $arrF[0]."\t".fnFormatType($sType)."\t".$sFeature."\t".$arrF[3]."\n");

}

function fnType2FileName($s) {
	global $sOutPrefix, $arrMainModels, $arrCandidateModels, $arrProtein2GenomeModels;
	if (array_key_exists($s , $arrMainModels) ) {
		return $sOutPrefix."_main.gff";
	}

	if (array_key_exists($s , $arrCandidateModels) ) {
		return $sOutPrefix."_candidate.gff";
	}

	if (array_key_exists($s , $arrProtein2GenomeModels) ) {
		return $sOutPrefix."_protein2genome.gff";
	}

	return false;
}

function fnFormatType($s) {
	return preg_replace("/[^A-z0-9_]/", "_", $s);
}

function fnRenameFeature($s) {
	if ($s == "match_part") {
		return "exon";
	}

	if (strpos( $s , "match" ) !== false) {
		return "mRNA";
	} 

	return $s;
}

?>
