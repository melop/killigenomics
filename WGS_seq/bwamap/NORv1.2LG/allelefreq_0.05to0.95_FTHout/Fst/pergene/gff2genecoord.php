<?php
require_once(dirname(__FILE__) . "/lib.php");
$sGFF = "/beegfs/group_dv/home/RCui/killifish_genomes/denovo/discovardenovo/NOR/LG/v1.2/NORv1.2.LG.gff";
$sFasta = "/beegfs/group_dv/home/RCui/killifish_genomes/denovo/discovardenovo/NOR/LG/v1.2/Both_Crosses.fasta";
$sOrthologList = "/beegfs/group_dv/home/RCui/killifish_genomes/annotation/UPhO_final/assigngenesymbol/killi_orthologs_table.ensemblinformed.txt";
$sSp = "NOR";
$sOut = "genecoord.txt.gz";

$hOut = popen("gzip > $sOut",'w');

$oGTF = new MappedCDSGFF();
$oGTF->LoadGFF($sGFF , $sFasta , 'CDS' ); //multiple genome mode

/*============ Load the ortholog definitions  ======================================*/
$hOrthologList = fopen($sOrthologList , "r");
$arrOrthologs = array();
$arrOrthologMeta = array();
$arrSpp = array();
$nCol = 0;

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
		foreach($arrSpp as $nIdx => $sSpname) {
			if ($sSpname == $sSp) {
				$nCol = $nIdx;
				echo("Found $sSp in header column $nCol\n");
				continue 2;
			}
		}
		die("Fail to find $sSp in the header.\n");
		continue; // skip header
	}
	
	$arrFields = explode("\t", $sLn);
	$arrOrthologs[$arrFields[0]] =  array_slice($arrFields , 3)[$nCol];
	//$arrOrthologMeta[] = array_slice($arrFields , 0, 3);
}

echo("Loaded ". count($arrOrthologs) ." ortholog definitions\n" );
#print_r($arrOrthologs);

foreach($arrOrthologs as $sGroupId => $sSpGeneId) {
	$sRNAID = $oGTF->fnSpGeneId2RNAID($sSpGeneId);
	if ($sRNAID === false ) {
		echo("warning: spgeneid $sSpGeneId not found in gff\n");
		continue;
	}

	$omRNA = $oGTF->GetmRNA('maker' , $sRNAID);

	//print_r($omRNA);

	//die();

	foreach($omRNA['CDSs'] as $nStart => $nEnd) {
		fwrite($hOut, "$sGroupId\t$sSpGeneId\t$sRNAID\t".$omRNA['strand']."\t".$omRNA['scf']."\t$nStart\t$nEnd\n");
	}
}

?>
