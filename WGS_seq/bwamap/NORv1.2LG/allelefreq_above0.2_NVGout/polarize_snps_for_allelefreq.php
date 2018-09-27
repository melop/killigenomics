<?php
require_once(dirname(__FILE__) . "/lib.php");

//this script filters the frequencies estimated by GFE package.
//for computing pair-wise LD.
//This script produces an output format conforming to the "l" output format of GFE
/*
Average coverage per population.
The population coverage is x * num_samples
Exclude coverage > x * num_samples * 2
  Group.1        x    NormX
1  ORTDRY 2.229171 1.000000
2  ORTWET 2.230603 1.000642
3  RACDRY 2.572075 1.153826
4  RACWET 2.251057 1.009818

*/
$sDIR = $argv[1];//"ORTDRY";
$nHWCriticalValue = 3.841; // critical chisquare value for p = 0.05 , ignore HWE, because this is haploid data
$nPolyCriticalValue = 3.841; // critical chisquare value for p = 0.05
$nMinPopCoverage = 20;
$nMaxPopCoverage = floatval($argv[2]);// 2.229171 * 60 * 2;
$nMinEffectiveChromosomes = 20;
$nMinEffectiveIndividuals = 10;
$nParts = 100;
$nMinMinorAlleleFreq = 0.2;  //min allele freq. if lower than this frequency it is ignored (won't be treated as fixed ancestral).
$nMaxMinorAlleleFreq = 0.95; //max allele freq.


$nMinMajorAlleleFreq = 0.95; // will be treated as fixed if above this frequency.
$sOutgroupFasta = "/beegfs/group_dv/home/RCui/killifish_genomes/denovo/discovardenovo/NOR/LG/v1.2/NVG.chr.fa";

$oOutgroup = new FastaHack();
$oOutgroup->SetDb($sOutgroupFasta);

//$sGFEInput = "$sDIR/My_in_GFE.txt";

$sOut = "$sDIR/ancestral_polarized_for_allelefreq.txt";
$sOutInvariableAnc = "$sDIR/ancestral_polarized_for_allelefreq.invariable_anc.txt";
$sOutInvariableDer = "$sDIR/ancestral_polarized_for_allelefreq.invariable_der.txt";
$hOut = fopen($sOut , 'w');
$hOutInvariableAnc = fopen($sOutInvariableAnc , 'w');
$hOutInvariableDer = fopen($sOutInvariableDer , 'w');

$arrPassLoci = array(); // records the selected locus's positions
$arrOpenFiles = array();
for($nPart=1;$nPart<=$nParts;$nPart++) {
	$sInFile = "$sDIR/out_f_GFE_".$nPart."_of_".$nParts.".txt";
	if (!file_exists($sInFile)) {
		die("$sInFile not found\n");
	}

	$arrOpenFiles[] = fopen( "$sDIR/out_f_GFE_".$nPart."_of_".$nParts.".txt", 'r');
}


$nEOFCount = 0;
$nCurrPart = 0;
$bHeaderWritten = false;
$bHeaderWrittenAnc = false;
$bHeaderWrittenDer = false;
echo("reading...\n");
while ($nEOFCount < $nParts ) {

	if ($nCurrPart >= $nParts) {
		$nCurrPart = 0;
	}

	if (feof($arrOpenFiles[$nCurrPart])) {
		$nCurrPart++;
		continue;
	}

	$sLn = fgets($arrOpenFiles[$nCurrPart]);
	if ($sLn === false ) {
		echo($nCurrPart." ended\n");
		$nEOFCount++;
		$nCurrPart++;

		continue;
	}

	$arrF = explode("\t" , trim($sLn));
	$arrOutput = array( 'scaffold' => $arrF[0] , 'site' => $arrF[1], 'ancestralbase' => 'N', 
				'derivedbase' => 'N' , 'popcov' => $arrF[5] , 
				'effective_chr' => $arrF[6], 'effective_ind' =>$arrF[7],
				'ancestral_freq' => 0, 'derived_freq' => 0, 'hom_anc' => 0, 'het' => 0 , 'hom_der' => 0, 'pol_llstat' => $arrF[20] , 'HWE_llstat' =>  $arrF[21] );
	$arrF[21] = ($arrF[21] =='NA')? 0 : $arrF[21];
	$arrF[20] = ($arrF[20] =='NA')? 0 : $arrF[20];
	$arrF[11] = ($arrF[11]=='NA')? 0 : $arrF[11];

	if ($arrF[0] == 'scaffold' || $arrF[3] == 'NA' || $arrF[5] < $nMinPopCoverage 
		|| $arrF[5] > $nMaxPopCoverage || $arrF[6] < $nMinEffectiveChromosomes || $arrF[7] < $nMinEffectiveIndividuals 
		|| $arrF[21] > $nHWCriticalValue
	 ) {
		$nCurrPart++;
		//echo("exclude: $sLn\n");
		continue;
	}

	$arrBaseFreq = array( $arrF[3] => 1 );
	if ($arrF[4] != 'NA' && $arrF[20] > $nPolyCriticalValue) { // reliably polymorphic
		$arrBaseFreq[$arrF[4]] = $arrF[11];
		$arrBaseFreq[$arrF[3]] = 1-$arrF[11];
	}

	$arrBases = array_keys($arrBaseFreq ); 
	$sOutgroupBase = $oOutgroup->GetCDSRegion($arrF[0],$arrF[1], $arrF[1] );

	if ($sOutgroupBase == 'N') {
		$nCurrPart++;
		continue; // this site contains 3 different states, exclude.
	}

	if ( (!array_key_exists($sOutgroupBase , $arrBaseFreq)) && count($arrBaseFreq) == 2 ) {
		$nCurrPart++;
		continue; // this site contains 3 different states, exclude.
	}

	$hToWrite = $hOut;
	if ( (array_key_exists($sOutgroupBase , $arrBaseFreq)) && ( count($arrBaseFreq) == 1 || $arrBaseFreq[$sOutgroupBase] > $nMinMajorAlleleFreq ) ) { //fixed ancestral sites.
		/*
		$arrOutput['ancestralbase'] = $sOutgroupBase;
		$arrOutput['derivedbase'] = 'N';
		$arrOutput['ancestral_freq'] = 1;
		$arrOutput['derived_freq'] = 0;
		$arrOutput['hom_anc'] = 1;
		$arrOutput['het'] = 0;
		$arrOutput['hom_der'] = 0;
		*/
		$arrOutput['ancestralbase'] = $sOutgroupBase;
		unset($arrOutput['derivedbase']);
		unset($arrOutput['ancestral_freq']);
		unset($arrOutput['derived_freq']);
		unset($arrOutput['hom_anc']);
		unset($arrOutput['het']);
		unset($arrOutput['hom_der']);
		unset($arrOutput['pol_llstat']);
		unset($arrOutput['HWE_llstat']);

		$hToWrite = $hOutInvariableAnc; //fixed ancestral sites
	} else if ( ((!array_key_exists($sOutgroupBase , $arrBaseFreq)) && count($arrBaseFreq) == 1) || ( array_key_exists($sOutgroupBase , $arrBaseFreq) && (1-$arrBaseFreq[$sOutgroupBase]) > $nMinMajorAlleleFreq ) ) { //divergent site with the outgroup species.
		/*$arrOutput['ancestralbase'] = $sOutgroupBase;
		$arrOutput['derivedbase'] = $arrBases[0];
		$arrOutput['ancestral_freq'] = 0;
		$arrOutput['derived_freq'] = 1;
		$arrOutput['hom_anc'] = 0;
		$arrOutput['het'] = 0;
		$arrOutput['hom_der'] = 1;
		*/

		if ( array_key_exists($sOutgroupBase , $arrBaseFreq) &&  (1-$arrBaseFreq[$sOutgroupBase]) > $nMinMajorAlleleFreq  ) {
			//exclude.
			$nCurrPart++;
			continue;
		}

		$arrOutput['ancestralbase'] = $sOutgroupBase;
		$arrOutput['derivedbase'] = $arrBases[0];
		unset($arrOutput['ancestral_freq']);
		unset($arrOutput['derived_freq']);
		unset($arrOutput['hom_anc']);
		unset($arrOutput['het']);
		unset($arrOutput['hom_der']);
		unset($arrOutput['pol_llstat']);
		unset($arrOutput['HWE_llstat']);

		$hToWrite = $hOutInvariableDer; //divergent sites.
	} else if ((array_key_exists($sOutgroupBase , $arrBaseFreq)) && count($arrBaseFreq) == 2) {
		$arrDiffBase = array_values(array_diff($arrBases , array( 0 => $sOutgroupBase) )); 
		if (!array_key_exists(0, $arrDiffBase) ) {
			print_r($arrBases);
			print_r($arrBaseFreq);
			print_r($arrDiffBase);	
			echo($sLn."\n");
			die("error occured!\n");
		}

		$sDeriveBase = $arrDiffBase[0];

		//check of derive allele freq is within the requested range of allele frequencies:

		if ($arrBaseFreq[$sDeriveBase] < $nMinMinorAlleleFreq || $arrBaseFreq[$sDeriveBase] > $nMaxMinorAlleleFreq ) {
			//exclude.
			$nCurrPart++;
			continue;
		}

		$arrOutput['ancestralbase'] = $sOutgroupBase;
		$arrOutput['derivedbase'] = $sDeriveBase;
		$arrOutput['ancestral_freq'] = $arrBaseFreq[$sOutgroupBase];
		$arrOutput['derived_freq'] = $arrBaseFreq[$sDeriveBase];

		if ($arrBaseFreq[$sOutgroupBase] > 0.5) { //if the ancestral base is the major allele
			$arrOutput['hom_anc'] = $arrF[16];
			$arrOutput['het'] = $arrF[17];
			$arrOutput['hom_der'] = $arrF[18];
		} else {
			$arrOutput['hom_anc'] = $arrF[18];
			$arrOutput['het'] = $arrF[17];
			$arrOutput['hom_der'] = $arrF[16];
		}

	}

	//scaffold	site	ref_nuc	n1	n2	pop_coverage	best_p	best_error	pol_llstat	HWE_llstat
	if (!$bHeaderWritten && ( $hToWrite == $hOut) ) {
		fwrite($hOut , implode("\t" , array_keys($arrOutput)) . "\n");
		$bHeaderWritten = true;
	}

	if (!$bHeaderWrittenAnc && ($hToWrite == $hOutInvariableAnc) ) {
		fwrite($hOutInvariableAnc , implode("\t" , array_keys($arrOutput)) . "\n");
		$bHeaderWrittenAnc = true;
	}

	if (!$bHeaderWrittenDer && ($hToWrite == $hOutInvariableDer) ) {
		fwrite($hOutInvariableDer , implode("\t" , array_keys($arrOutput)) . "\n");
		$bHeaderWrittenDer = true;
	}


	fwrite($hToWrite , implode("\t" , $arrOutput) . "\n");
	$nCurrPart++;
}

echo(count($arrPassLoci) ." loci passing\n" );

?>
