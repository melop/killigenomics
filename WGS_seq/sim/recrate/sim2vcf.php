<?php
//directly convert the 1/0 genotypes into a vcf file for smc++

$sIn = "simORT.out.txt.raw.txt";
$sSubPop = "ORTWET"; //"ORTDRY"; //put "both" to output both pops.
//$sIn = "simRAC.out.txt.raw.txt";
$sOut = "bgzip > vcfs/$sIn.$sSubPop.vcf.gz";
$sAncChrFile = "chr1.NOR.replaced.fa";
$sSampleNames = "cat /beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/NORv1.2LG/mapped/*.bamlist.txt";


$arrSampleNames = fnLoadSampleName($sSampleNames);

//print_r($arrSampleNames );
//die();
list($sChrName, $sAncSeq) = fnReadAnc($sAncChrFile);

list($nVarSites, $arrVarSites, $arrDevBase, $arrInd) = fnReadSim($sIn);

//write header.
$hOut = popen($sOut , 'w');

fwrite($hOut, "##fileformat=VCFv4.2\n");
fwrite($hOut, "##reference=". realpath($sAncChrFile) ."\n");
fwrite($hOut, "##contig=<ID=$sChrName,length=".strlen($sAncSeq).">\n");
fwrite($hOut, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" );

foreach($arrInd as $sPop => &$arrI) {
	if ( $sPop != $sSubPop && $sSubPop!='both') continue;
	$arrIndNames = array_keys($arrI);
	fwrite($hOut, "\t" . implode("\t", $arrIndNames) );
}
	fwrite($hOut, "\n" );

for($nSite=0; $nSite<$nVarSites; $nSite++) {
	$nSitePos = $arrVarSites[$nSite];
	$sDevBase = $arrDevBase[$nSite];
	$sRefBase = $sAncSeq[$nSitePos-1];
	fwrite($hOut, "$sChrName\t$nSitePos\t.\t$sRefBase\t$sDevBase\t999\tPASS\t.\tGT");

	foreach($arrInd as $sPop => &$arrI) {
		if ( $sPop != $sSubPop && $sSubPop!='both') continue;
		foreach($arrI as $s => &$arrGen) {
			$arrGeno = fnStr2Genotype($arrGen[0][$nSite], $arrGen[1][$nSite]);
			fwrite($hOut, "\t". $arrGeno[0]);
		}
	}
	fwrite($hOut, "\n" );

}

function fnReadAnc($sFile) {
	$hIn = fopen($sFile, 'r');
	$nRecordCount = 0;
	$sSeq = '';
	$sName = '';
	while( false !== ($sLn = fgets($hIn)) ) {
		$sLn = trim($sLn);
		if ($sLn[0] == '>') {
			$sName = substr($sLn , 1);
			$nRecordCount++;
			if ($nRecordCount > 1) {
				die("Error: the anc.fa file should contain a single fasta record for the ancestral sequence.\n");
			}
			continue;
		}

		if (preg_replace("/[ATCGatcg]/", "", $sLn) != '') {
			die("Error: DNA sequence can only contain ATCG, no ambiguous character is allowed.\n$sLn\n");
		}
		$sSeq .= $sLn;
	}

	return array($sName , strtoupper($sSeq));
}

function fnReadSim($sIn) {
	global $arrSampleNames;
	$hIn = fopen($sIn , 'r');

	$bReadPos = false;
	$bReadDeriveBase = false;
	$bReadHap = false;

	$arrPos = array();
	$arrDevBases = array();

	$sPop = "";
	$nHapCount = 0;
	$sHap1 = "";

	$arrRet = array();

	$nPos = 0;

	while(false !== ($sLn = fgets($hIn))) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		if ($sLn == '>Positions') {
			$bReadPos = true;
			continue;
		}

		if ($sLn == '>DerivedBases') {
			$bReadDeriveBase = true;
			continue;
		}

		if ($bReadPos) {
			$arrPos = explode(',', $sLn);
			$nPos = count($arrPos);
			$bReadPos = false;
			continue;
		}

		if ($bReadDeriveBase) {
			$arrDevBases = explode(',' , $sLn);
			if ($nPos != count($arrDevBases) ) {
				die("Pos count not same . devbase $nPos vs ". count($arrDevBases)."\n");
			}
			$bReadDeriveBase = false;
			$bReadHap = true;
			continue;
		}

		if ($bReadHap && $sLn[0] == '>') {
			$sHapName = substr($sLn, 1);
			$arrF = explode('.', $sHapName);
			$sThisPop = $arrF[0];
			if ($sThisPop != $sPop) {
				$nHapCount = 0;
				$sPop = $sThisPop ;
			}

			$nHapCount++;
			continue;
		}

		if ($nHapCount % 2 == 0) { // 
			$sHap2 = $sLn;	

			if ($nPos != strlen($sHap2) ) {
				die("Pos count not same. hap2, $nPos vs. ". strlen($sHap2)." \n");
			}

			$nInd = $nHapCount / 2;

			$sSName = $arrSampleNames[$sPop][$nInd-1];

			echo("$sPop $nInd $sSName\n");

			if (!array_key_exists($sPop , $arrRet)) {
				$arrRet[$sPop] = array();
			}

			$arrRet[$sPop][$sSName] = array($sHap1, $sHap2);
			$sHap1 = $sHap2 = '';
			continue;

		} else {
			$sHap1 = $sLn;
			if ($nPos != strlen($sHap1) ) {
				die("Pos count not same. hap1, $nPos vs. ". strlen($sHap1)." \n");
			}

			continue;	
		}

	}

	return array($nPos, $arrPos,$arrDevBases, $arrRet);
}

function  fnLoadSampleName($sCmd) {
	$h = popen($sCmd , 'r');
	$arrSamples = array();
	while(false !== ($sLn = fgets($h))) {
		$sLn = trim($sLn);
		$sBam = basename($sLn);
		$sName = str_replace('.bam', '', $sBam);
		$arrF = explode('_', $sName);
		$sPop = $arrF[0];
		if (!array_key_exists($sPop, $arrSamples) ) {
			$arrSamples[$sPop] = array();
		}

		$arrSamples[$sPop][] = $sName;
	}

	return $arrSamples;
}

function fnStr2Genotype($sHap1, $sHap2) {

	$arrRet = array();
	if (strlen($sHap1) != strlen($sHap2) ) {
		die("haplotype lengths differ.\n");
	}

	for($i=0;$i<strlen($sHap1);$i++) {
		$n1 = $sHap1[$i];
		$n2 = $sHap2[$i];
		$arrRet[] = ($n1!=$n2)? "0/1" : "$n1/$n2";
	}

	return $arrRet;
}
?>
