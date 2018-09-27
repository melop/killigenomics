<?php
ini_set('memory_limit', -1);
//echo(fnComputeFST(array(10,0,2), array(20,0,0)));
//die();
$sPop1 = $argv[1];
$sPop2 = $argv[2];
//$sPop = "testpop";
$sOutFolder = $argv[3];

$sInputPrefix = "/beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/bwamap/NORv1.2LG/allelefreq_lowfreq_excl/";
//$sInputPrefix = "./";

exec("mkdir -p $sOutFolder");
$hOut = popen("gzip > $sOutFolder/fst.txt.gz", 'w');

$arrVarSites1 = array( 'PolyMorphic' => "$sInputPrefix/$sPop1/ancestral_polarized_for_allelefreq.txt", 'InvarDerived' => "$sInputPrefix/$sPop1/ancestral_polarized_for_allelefreq.invariable_der.txt",
 'InvarAncestral' => "$sInputPrefix/$sPop1/ancestral_polarized_for_allelefreq.invariable_anc.txt");
$arrVarSites2 = array( 'PolyMorphic' => "$sInputPrefix/$sPop2/ancestral_polarized_for_allelefreq.txt", 'InvarDerived' => "$sInputPrefix/$sPop2/ancestral_polarized_for_allelefreq.invariable_der.txt",
 'InvarAncestral' => "$sInputPrefix/$sPop2/ancestral_polarized_for_allelefreq.invariable_anc.txt");


exec("mkdir -p $sOutFolder");

do {

	$arrVariantSites1 = fnLoadAlleleFreq($arrVarSites1);
	$arrVariantSites2 = fnLoadAlleleFreq($arrVarSites2);

	if ($arrVariantSites1 === false || $arrVariantSites2 === false) {
		echo("All done.\n");
		break;
	}

	$arrChrs = array_intersect(array_keys($arrVariantSites1['sites']) , array_keys($arrVariantSites2['sites']));

	foreach($arrChrs as $sChr) {
		$arrPos = array_intersect(array_keys($arrVariantSites1['sites'][$sChr]) , array_keys($arrVariantSites2['sites'][$sChr]));

		foreach($arrPos as $nPos) {
			$arrCounts1 = $arrVariantSites1['sites'][$sChr][$nPos]['GTcount'];
			$arrCounts2 = $arrVariantSites2['sites'][$sChr][$nPos]['GTcount'];
			$arrFst = fnComputeFST($arrCounts1 , $arrCounts2);
			if ($arrFst !== false) {
				fwrite($hOut, "$sChr\t$nPos\t". implode("\t", $arrFst) ."\n");
			}
		}
	}

} while (true);

function fnComputeFST($arrCounts1 , $arrCounts2) {
	$nPop1P = ($arrCounts1[2] * 2 + $arrCounts1[1]) / (2* ($arrCounts1[0]  + $arrCounts1[1] + $arrCounts1[2]));

	$nPop2P = ($arrCounts2[2] * 2 + $arrCounts2[1]) / (2* ($arrCounts2[0]  + $arrCounts2[1] + $arrCounts2[2]));

	if ( ($nPop1P == $nPop2P) && ($nPop1P==0 || $nPop1P == 1) ) {
		return false; //undefined.
	}

	#echo("$nPop1P $nPop2P\n");

	$nPop1Ind = $arrCounts1[0] +  $arrCounts1[1] + $arrCounts1[2];
	$nPop2Ind = $arrCounts2[0] +  $arrCounts2[1] + $arrCounts2[2];

	$h1 = $arrCounts1[1] / $nPop1Ind;
	$h2 = $arrCounts2[1] / $nPop2Ind;

	#echo("$h1 $h2\n");

	#based on AND C. C. COCKERHAM 1984

	$r = 2; //number of pops

	$n1 = $nPop1Ind * 2;
	$n2 = $nPop2Ind * 2;

	$n_bar = ($n1 + $n2) / $r; //the average sample size

	$n_c = ( ($r * $n_bar) - (pow($n1,2) + pow($n2,2))/ ($r * $n_bar ) ) / ($r-1);

	$p_bar =  ($nPop1P * $n1 + $nPop2P * $n2) / ($r*$n_bar);  //the average sample frequency of allele A

	$s_square = ( $n1 * (pow($nPop1P - $p_bar, 2 )) + $n2 * (pow($nPop2P - $p_bar , 2 )) ) / (($r - 1) * $n_bar); // the sample variance of allele A frequencies over populations

	$h_bar = ( $n1 * $h1 + $n2 * $h2) / ($r * $n_bar); //, the average heterozygote frequency for allele A.

	$a = ($n_bar/$n_c) * ( $s_square - (1/($n_bar - 1)) * ($p_bar*(1-$p_bar) - $s_square*($r-1)/$r - 0.25 * $h_bar )  );

	$b = ($n_bar / ($n_bar - 1) ) * ($p_bar*(1-$p_bar) - $s_square*($r-1)/$r - ((2*$n_bar - 1)/(4*$n_bar)) * $h_bar );

	$c = 0.5 * $h_bar;

	//echo("nbar $n_bar nc $n_c pbar $p_bar ssq $s_square hbar $h_bar\na $a b $b c $c ". pow($nPop1P - $p_bar ,2)."\n");

	if (($a + $b + $c) <= 0 ) {
		return false;
	}

	$Fst = $a / ($a + $b + $c);

	return array($Fst, $a, $b, $c);
}

function fnLoadAlleleFreq(&$arrVarSites) {

	if (!array_key_exists("file_handlers" , $arrVarSites)) {
		 $arrVarSites["file_handlers"] = array();
		foreach($arrVarSites as $sFileType => $sVarFile) {
			if ($sFileType == 'file_handlers') continue;
			echo("Open file $sVarFile \n");
			$arrVarSites["file_handlers"][$sFileType] = fopen($sVarFile , 'r'); // open the file handler, only once.
		}
	}

	$arrVariantSites= array();
	$nVarSitesCount = 0;
	$nDeriveFixSitesCount = 0;
	$nAncFixSitesCount = 0;
	$nEffectiveChrTotal = 0;

	$arrChrMaxCoord = array();

	echo("Loading allele frequencies...\n");
	$sCurrChr = "";
	$bFileNotEnd = false;
	foreach($arrVarSites as $sFileType => $sVarFile) {

		if ($sFileType == "file_handlers") {
			continue;
		}

		$hVarFile = $arrVarSites["file_handlers"][$sFileType]; // use the opened file handler.
		$nFilePointer = ftell($arrVarSites["file_handlers"][$sFileType]);
		
		while( false !== ($sLn = fgets($hVarFile) ) ) {
			$sLn = trim($sLn);
			$arrF = explode("\t", $sLn);
			if ($arrF[0] == 'scaffold') continue;

			//echo("currchr $sCurrChr arrF[0] $arrF[0]\n");
			if ($sCurrChr == '') {
				$sCurrChr = $arrF[0];
			} else if ($sCurrChr != $arrF[0]) {
				//reset to previous line, and return results
				fseek($arrVarSites["file_handlers"][$sFileType], $nFilePointer ,SEEK_SET);
				$bFileNotEnd = true;
				break;
			}
			$nFilePointer = ftell($arrVarSites["file_handlers"][$sFileType]);

			$arrLocus = array( 0, 0, 0); //number of counts for genotypes. hom-anc, het, hom-derived
			if ($arrF[0] == 'scaffold') continue;
			//if ($arrF[0] != 'chr1') break;

			$nMAF = 0; //default is fixed for ancestral
			$nEffectChrN = 0;

			if ($sFileType == 'PolyMorphic') {
				$nMAF = $arrF[8];
				//die($nMAF);
				$nVarSitesCount++;
				//$nEffectChrN = round($arrF[5]);
				$arrLocus[0] = round($arrF[5] * 0.5 * $arrF[9]);
				$arrLocus[1] = round($arrF[5] * 0.5 * $arrF[10]);
				$arrLocus[2] = round($arrF[5] * 0.5 * $arrF[11]);
				$nEffectChrN = ($arrLocus[0] + $arrLocus[1] + $arrLocus[2]) * 2;
			} else if ($sFileType == 'InvarDerived') {
				$nMAF = 1; //derived allele freq
				$nDeriveFixSitesCount++;
				$nEffectChrN = round($arrF[5]);
				$arrLocus[0] = 0;
				$arrLocus[1] = 0;
				$arrLocus[2] = $nEffectChrN / 2;
			} else {
				$nAncFixSitesCount++;
				$nEffectChrN = round($arrF[4]);
				$nMAF = 0; //derived allele freq
				$arrLocus[0] = $nEffectChrN / 2;
				$arrLocus[1] = 0;
				$arrLocus[2] = 0;

			}

			$nEffectiveChrTotal += $nEffectChrN;
			if (!array_key_exists($arrF[0], $arrVariantSites)) {
				$arrVariantSites[$arrF[0]] = array();
			}

			$arrVariantSites[$arrF[0]][$arrF[1]] = array( 'GTcount' => $arrLocus );

		}

	}

	if (!$bFileNotEnd ) return false; //all files have reached the end.

	$nTotalSites = ($nVarSitesCount+$nDeriveFixSitesCount+$nAncFixSitesCount);
	$nAvgEffChr = round($nEffectiveChrTotal/ $nTotalSites );
	$nTheta = $nVarSitesCount/($nVarSitesCount+$nDeriveFixSitesCount+$nAncFixSitesCount)/fnHarmonicNum(1, $nAvgEffChr - 1 );
	echo("$sCurrChr: Loaded polymorphic $nVarSitesCount sites, $nDeriveFixSitesCount fixed, derived sites and $nAncFixSitesCount fixed ancestral sites.\nTotal = " . $nTotalSites ."; Theta = $nTheta; EffectChr = $nAvgEffChr\n");
	echo("sorting...\n");
	foreach($arrVariantSites as $sChr => &$arrChrVariantSites) {
		ksort($arrChrVariantSites);
		$arrPos = array_keys($arrChrVariantSites); 
		$arrChrMaxCoord[$sChr] = $arrPos[count($arrPos)-1]; 
	}

	return array( 'sites' => $arrVariantSites , 'maxcoord' => $arrChrMaxCoord );

}


function fnHarmonicNum($nLow, $nHigh) {
	$nSum = 0;
	for($i=$nLow; $i<= $nHigh; $i++) {
		$nSum += 1/$i;
	}

	return $nSum;
}

?>
