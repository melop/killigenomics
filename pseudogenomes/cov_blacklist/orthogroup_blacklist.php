<?php
//load the ortholog definition file,
//then load the black list for each reference genome
//flag the orthologs where at least one of the gene model is black listed.

$sOut = "orthogroup_blacklist.txt";
$sOrthologList = "../UPhO/assigngenesymbol/killi_orthologs_table.ensemblinformed.txt";
//symbolic link is /beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/export_protein_multiref_23_july_2017/orthologs.improved.txt

$hOrthologList = fopen($sOrthologList , "r");
$arrOrthologs = array();
$arrOrthologMeta = array();
$arrSpp = array();

$nLn = -1;
echo("Parsing ortholog definitions...\n");
while( ($sLn=fgets($hOrthologList))!==false ) {
	$sLn = trim($sLn);
	if ($sLn == "") {
		continue;
	}
	$nLn++;
	if ($nLn==0) {
		$arrSpp = array_slice( explode("\t" , $sLn) , 3);
	
		continue; // skip header
	}
	
	$arrFields = explode("\t", $sLn);
	$arrOrthologs[] =  array_combine( $arrSpp ,  array_slice($arrFields , 3));
	$arrOrthologMeta[] = array_slice($arrFields , 0, 3);
}

echo("Loaded ". count($arrOrthologs) ." ortholog definitions\n" );

$arrBlackList = array();
foreach($arrSpp as $sSp) {
	$hBlackList = fopen("$sSp"."_blacklist.txt" , 'r');
	$arrBlackList[$sSp] = array();
	while( false !== ($sLn = fgets($hBlackList) )) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		$arrF = explode("\t", $sLn);
		$nGeneId = $arrF[0] ;
		$nExcluded = ($arrF[count($arrF)-1] == 'TRUE');
		if ($nExcluded) {
			$arrBlackList[$sSp][$nGeneId] = true;
		}
	}

	echo("Blacklisted " . count($arrBlackList[$sSp]) . " gene models for species $sSp\n" );
}

$hOut = fopen($sOut , 'w');
$nExcCount = 0;
for($nOrth=0; $nOrth < count($arrOrthologs); $nOrth++ ) {
	$arrMeta = $arrOrthologMeta[$nOrth];
	$arrBadSp = array();
	foreach($arrOrthologs[$nOrth] as $sSp => $sGeneId) {
		if (array_key_exists($sGeneId , $arrBlackList[$sSp])) {
			$arrBadSp[] = $sSp;
		}
	}

	if (count($arrBadSp) > 0) {
		fwrite($hOut, implode("\t" , $arrMeta ) . "\t". implode(",", $arrBadSp) . "\n" );
		$nExcCount++;
	}
}

echo("Excluded $nExcCount groups out of ". count($arrOrthologs) ." ortholog groups \n");
?>
