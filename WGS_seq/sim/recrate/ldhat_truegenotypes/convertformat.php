<?php
$sPop = array_key_exists(1, $argv)? $argv[1]:"RACWET";
$sChr = array_key_exists(2, $argv)? $argv[2]:"chr1";
$sRoot = "/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/sim/recrate/vcfs";
$sBCF = "$sRoot/simORT.out.txt.raw.txt.$sPop.vcf.gz";
$sOutDir = "ldhat/$sPop/$sChr";
$sMatchToSites = "../ldhat/ldhat/$sPop/$sChr/loc.txt"; //match to site list

$nMaxAF = 0.9;
$nMinAF = 0.1;



$hLociList = popen("tail -n1 $sMatchToSites", 'r');
$arrMatchToLoci = explode(" " , trim(fgets($hLociList)));

foreach($arrMatchToLoci as $nIdx => &$nVal) {
	$nVal = intval(floatval($nVal) * 1000);
}

$arrMatchToLoci = array_flip($arrMatchToLoci);
$ak = array_keys($arrMatchToLoci);
echo($ak[0]. "\t". $ak[1]."\n");
echo("Matching to ". count($arrMatchToLoci)." loci...\n");
exec("mkdir -p $sOutDir");
$arrIndividuals = array();
$arrPos = array();
//parse genotypes:
$hIn = popen("bcftools1.2 view $sBCF " , 'r');
$nChrLen = -1;
while( false !== ($sLn = fgets($hIn) ) ) {
	$sLn = trim($sLn);
	if ($sLn == '') continue;
	if ($nChrLen==-1 && substr($sLn, 0, 9) == '##contig=') {
		//try parsing chromosome length:
		preg_match("/ID=(\S+),length=(\d+)/", $sLn, $arrM);
		if (count($arrM) != 3 ) {
			continue;
		}
		if ($arrM[1] == $sChr ) {
			$nChrLen = $arrM[2];
		}
		continue;
	}

	if ($sLn[0] == '#') continue;
	
	if ($nChrLen == -1) {
		die("Error: Chromosome length is not found in the header of the bcf file!\n");
	}

	$arrF = explode("\t", $sLn);
	$sThisChr = $arrF[0];
	if ($sThisChr != $sChr) continue;

	$nPos = $arrF[1];

	$arrRawGenotypes = array_slice($arrF , 9);
	fnPushLocus($nPos, $arrRawGenotypes);
}

echo(count($arrPos) . " Loci\n" );
$hSites = fopen("$sOutDir/sites.txt"  , 'w');
$hLoci = fopen("$sOutDir/loc.txt"  , 'w');

fwrite($hSites , count($arrIndividuals) ." ".count($arrPos)." 2\n"); 

foreach($arrIndividuals as $nIdx => $arrSeq) {
	fwrite($hSites , ">Genotype$nIdx\n" . implode('', $arrSeq)."\n" ); 
}

fwrite($hLoci , count($arrPos) ." ".($nChrLen/1000)." L\n"); //write length in kb 
fwrite($hLoci , implode(' ', $arrPos) ."\n"); 

function fnPushLocus($nPos, $arrRawGenotypes) {
	global $arrIndividuals, $arrPos, $nMaxAF , $nMinAF, $arrMatchToLoci;



	if (!array_key_exists(intval($nPos) , $arrMatchToLoci) ) {
		return;
	}

	if (count($arrIndividuals) == 0) {
		$arrIndividuals = array_fill(0, count($arrRawGenotypes) , array() );
		echo(count($arrIndividuals) . " individuals\n" );
	} else if (count($arrIndividuals) != count($arrRawGenotypes) ) {
		die("Individual count became different!\nPos: $nPos\n");
	}

	$arrEncodedGenotype = array();
	foreach($arrRawGenotypes as $nIdx => $sGenotype) {
		$nEncodedGeno = ($sGenotype=='0/0')? 0 : (($sGenotype=='1/1')? 1: (($sGenotype=='1/0' || $sGenotype=='0/1')? 2 : '?') ); //hom: 0, 1; het: 2; missing: ?
		$arrEncodedGenotype[] = $nEncodedGeno;
	}

	$nAF = array_sum($arrEncodedGenotype) / (count($arrEncodedGenotype)*2);

	if ($nAF < $nMinAF || $nAF >$nMaxAF ) {
		echo("exclude $nPos, AF = $nAF\n");
		return; // af doesn't match
	}

	foreach($arrEncodedGenotype as $nIdx => $sGenotype) {
		$arrIndividuals[$nIdx][] = $sGenotype;
	}

	$arrPos[] = $nPos/1000;

}

?>
