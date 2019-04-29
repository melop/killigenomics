<?php
$sIndPopMap = "samples.panel";
$sOutDir = "polarized";
$BCFTOOLS = "/beegfs/group_dv/software/source/bcftools-1.6/bcftools";
$sOutgroupSeq = "/beegfs/group_dv/home/RCui/killifish_genomes/human_alfred/gorilla/gorilla_pseudo.fa"; // this is the outgroup pseudogenome (exact coordinates with the ref).
$sRefGenome = "/beegfs/group_dv/data/human_GRCh37.p13/GRCh37.p13_chr_rename.fa"; // the genome used for mapping.

$nMinAlleleFreq = 0.045 ; //min is 5%, assuming 10 individual or 20 chromosomes, if lower than this, it's treated as ancestrally fixed
$nMinChr = 20;

$nThisPart = 0;
$nTotalPart = 1;

$sChrPrefix = ""; //the prefix in vcf file, typically "chr" 

$arrWarnings = array();

while(count($argv) > 0) {
    $arg = array_shift($argv);
    switch($arg) {
        case '-o':
            $sOutDir  = trim(array_shift($argv));
            break;
        case '-N':
            $nTotalPart  = intval(trim(array_shift($argv)));
            break;
        case '-f':
            $nThisPart  = intval(trim(array_shift($argv)));
            break;

    }
}



exec("mkdir -p $sOutDir");

$arrRefGenome = fnParseFasta($sRefGenome);
$arrOutgroupGenome = fnParseFasta($sOutgroupSeq);




$arrChrs = range(1, 22);
$arrChrs[] = "X";
$arrChrs[] = "Y";

$arrVCFs = "../";
$arrPops = fnParsePops($sIndPopMap);

$arrDivergentPop2Outgroup = array(); //key - pop, values = array with divergent sites from outgroup for each pop.


//print_r($arrPops);
$nJobCount = -1;
foreach($arrChrs as $sChr) {
	$nJobCount++;
	if ($nJobCount % $nTotalPart != $nThisPart) continue;

	$arrDivergentRefOutgroup = fnGetDivergentSites($arrRefGenome , $arrOutgroupGenome, $sChr); // get divergent sites, 0-based
	$arrDivergentPop2Outgroup = fnInitializeDivergentPop2Outgroup($arrDivergentRefOutgroup, $arrPops, $sChr) ; 

        $arrFile = glob("../ALL.chr$sChr".".*genotype*.gz");
        if (count($arrFile) != 1) {
                print_r($arrFile);
                echo("../ALL.chr$sChr".".*genotype*.gz\n");
                die("Error grabing vcf for chr $sChr\n");
        }

	$sF = $arrFile[0];

	echo("Reading $sF...\n");

	$h = popen("$BCFTOOLS view -r $sChrPrefix$sChr $sF" , 'r');

	while(false !== ($sLn = fgets($h) ) ) {
		$sLn = trim($sLn);
		$arrF = explode("\t", $sLn);
		if (count($arrF) < 9) continue;
		if ($arrF[0] == "#CHROM") {
			fnAssignPopCol($arrF, $arrPops, $sChr);
			continue;
		}

		//process the line:
		fnProcess($arrF, $arrPops );
	}

}

function fnProcess($arrF, $arrPops ) {
	global $nMinAlleleFreq, $nMinChr, $arrRefGenome, $arrOutgroupGenome, $arrDivergentPop2Outgroup, $sChrPrefix;
	$sChr = $arrF[0];
	$sChr = str_replace($sChrPrefix , '', $sChr);

	$nSite = $arrF[1];

	if (strlen($arrF[3]) != 1 || strlen($arrF[4]) != 1 ) {
		return;
	}

	if ($arrF[3] == 'N' || $arrF[4] == 'N') {
		return;
	}

	if (!array_key_exists($sChr, $arrRefGenome) ) {
		Warn("$sChr not found in reference genome, skip\n");
		return; // this chromosome is not found in the genome
	}

//treat the gorilla sequence as the ancestral 
/*
	preg_match('/AA=(\\S)/', $arrF[7], $arrM);
	if (count($arrM)!=2) {
		return;
	}

	$sAncAllele = strtoupper($arrM[1]);
*/
	//assert reference allele is correct:
	if ($arrF[3] != $arrRefGenome[$sChr][$nSite-1]) {
		die("Error: $sChr:$nSite, reference base differ from fasta file! $arrF[3] vs. ".$arrRefGenome[$sChr][$nSite-1]."\n");
	}

	$sAncAllele = $arrOutgroupGenome[$sChr][$nSite-1]; //change to 0-based

	if ($sAncAllele == 'N') {
		return; //unable to polarize.
	}

	if ($sAncAllele != $arrF[3] &&  $sAncAllele != $arrF[4] ) {
		return; //triallelic sites
	}

	$bFlipFreq = false;
	$sDevAllele = $arrF[4];
	if ($sAncAllele == $arrF[4]) { // flip the alleles
		$bFlipFreq = true;
		$sDevAllele = $arrF[3];
	}

	foreach($arrPops as $sPop => $arrPopInfo) {
		$hOut = $arrPopInfo['stdout'];
		$arrGenotypes = array_intersect_key($arrF , $arrPopInfo['samplescolskeys'] );
		list( $n1Count, $nTotalCount, $nEffectiveInd, $nHom0, $nHet, $nHom1) = fnGenotype2Freq($arrGenotypes); //allele freq for "1"s
		$nDevAF = $n1Count / $nTotalCount;
		if ($bFlipFreq) {
			$nDevAF = 1 - $nDevAF;
			$nTmp = $nHom0;
			$nHom0 = $nHom1;
			$nHom1 = $nTmp;
		}

		if ($nDevAF < $nMinAlleleFreq || $nTotalCount < $nMinChr) {
			if (($nDevAF < $nMinAlleleFreq)  && 
				array_key_exists($sChr, $arrDivergentPop2Outgroup[$sPop] ) && array_key_exists($nSite-1, $arrDivergentPop2Outgroup[$sPop][$sChr] ) ) {
				unset($arrDivergentPop2Outgroup[$sPop][$sChr][$nSite-1]); //don't count as divergent if this pop is fixed for ancestral allele.
			}

			continue;
		}

		if ($nDevAF == 1 ) { // if the derived freq is 1, then set as divergent
			if (!array_key_exists($sChr, $arrDivergentPop2Outgroup[$sPop]) )  {
				$arrDivergentPop2Outgroup[$sPop][$sChr] = array();
			}

			$arrDivergentPop2Outgroup[$sPop][$sChr][$nSite-1] = array($sAncAllele , $sDevAllele, $nTotalCount * 10, $nTotalCount, $nEffectiveInd);

			continue;
		}	

		$nAncFreq = 1 - $nDevAF;

		fwrite($hOut, "$sChr\t$nSite\t$sAncAllele\t$sDevAllele\t".($nTotalCount*10)."\t$nTotalCount\t$nEffectiveInd\t$nAncFreq\t$nDevAF\t$nHom0\t$nHet\t$nHom1\n");
	}


}

	//write divergent sites
foreach($arrDivergentPop2Outgroup as $sPop => $arrPopDiverge) {


	foreach($arrPopDiverge as $sChr => $arrChrDiverge) {
		$hOut = popen("gzip -c > $sOutDir/der_".$sPop.".chr".$sChr.".txt.gz" , 'w');
		fwrite($hOut, "scaffold\tsite\tancestralbase\tderivedbase\tpopcov\teffective_chr\teffective_ind\n");

		foreach($arrChrDiverge as $n0Site => $arrInfo) {
			$nSite = $n0Site + 1; //change 0-based to 1 based for output
			if (count($arrInfo) == 2) { // default write
				fwrite($hOut, "$sChr\t$nSite\t$arrInfo[0]\t$arrInfo[1]\t-1\t-1\t-1\n");
			} else if (count($arrInfo) == 5) {
				fwrite($hOut, "$sChr\t$nSite\t$arrInfo[0]\t$arrInfo[1]\t$arrInfo[2]\t$arrInfo[3]\t$arrInfo[4]\n");
			}
		}
	}
}



function fnGenotype2Freq($arrGenotypes) {
	$n1Count = 0;
	$nTotalCount = 0;
	$nEffectiveInd = 0;
	$nHom0 = 0;
	$nHet = 0;
	$nHom1 = 0;

	foreach($arrGenotypes as $sG) {
		$arrF = preg_split('/[\|\/\\\\]/', substr($sG,0,3) , -1);
		if (count($arrF) > 2) {
			die("Genotype format unexpected: $sG\n");
		}

		$bIndCounted = false;
		$bMissing = false;
		$nSum = 0;
		foreach($arrF as $nVal) {
			if ($nVal == '.') {
				$bMissing = true;
				continue;
			}

			if ($nVal == 1) {
				$n1Count++;
			}

			$nSum += $nVal;

			$nTotalCount++;
			if (!$bIndCounted) {
				$nEffectiveInd++;
				$bIndCounted = true;
			}
		}

		if (!$bMissing ) {
			if ($nSum == 0) {
				$nHom0++;
			}
			if ($nSum == 1) {
				$nHet++;
			}
			if ($nSum == 2) {
				$nHom1++;
			}

		}
	}

	return array($n1Count , $nTotalCount, $nEffectiveInd, $nHom0, $nHet, $nHom1);
}

function fnAssignPopCol($arrF, &$arrPops, $sChr) {
	global $sOutDir;
	foreach($arrPops as $sPop => $arrSamples) {
		foreach($arrSamples as $sIdx => $s) {
			$arrPops[$sPop]['samplescols'][$sIdx] = -1;
		}

		foreach($arrF as $nCol => $sSampleID) {
			if (array_key_exists($sSampleID, $arrSamples)) {
				$arrPops[$sPop]['samplescols'][$sSampleID] = $nCol;
			}
		}

		$arrPops[$sPop]['samplescolskeys'] = array_flip($arrPops[$sPop]['samplescols']);
		$arrPops[$sPop]['stdout'] = popen("gzip -c > $sOutDir/".$sPop.".chr".$sChr.".txt.gz" , 'w');
		fwrite($arrPops[$sPop]['stdout'] , "scaffold\tsite\tancestralbase\tderivedbase\tpopcov\teffective_chr\teffective_ind\tancestral_freq\tderived_freq\thom_anc\thet\thom_der\n");

	}
}

function fnParsePops($sIndPopMap) {
	global $sOutDir;
	$h = fopen($sIndPopMap , 'r');
	$arrRet = array();
	while( false !== ($sLn = fgets($h) ) ) {
		$sLn = trim($sLn);
		$arrF = explode("\t", $sLn);
		if (count($arrF) < 2) continue;
		if ($arrF[1] == 'pop') continue;
		if (!array_key_exists($arrF[1], $arrRet) ) {
			$arrRet[$arrF[1]] = array('samplescols' => array(), 'stdout' => false);
			//$arrRet[$arrF[1]] = array('samplescols' => array(), 'stdout' => fopen("$sOutDir/".$arrF[1].".txt" , 'w'));

		}

		$arrRet[$arrF[1]][$arrF[0]]['samplescols'] = -1;
	}

	return $arrRet;
}

function fnParseFasta($sRefGenome) {
	//============================================================================
	$hRefGenome = fopen($sRefGenome, "r");
	$arrRefGenome = array(); //save ref genome into the RAM

	$sLn = "";

	$sSeqName = "";
	$sSeq = "";
	$bRecordedLnWidth = false;
	echo("Loading genome fasta...\n");
	do {
		$sLn = fgets($hRefGenome);

		if ($sLn==false) {
		        //fnProcessRecord($sSeqName, $sSeq);
			if ($sSeqName != '') {
		        	$arrRefGenome[$sSeqName] = $sSeq;
			}

		        break;
		}

		$sLn = trim($sLn);
		if (substr($sLn , 0,1) == ">") {
			if ($sSeqName != '') {
		        	$arrRefGenome[$sSeqName] = $sSeq;
			}

			//if ($sSeqName == '1') break;

		        $sSeqName = substr($sLn , 1);
		        $sSeq = "";
		        
		} else {
			if (!$bRecordedLnWidth) {
				$nFastLnWidth = strlen($sLn);
				$bRecordedLnWidth = true;
				echo("Fasta line width: $nFastLnWidth\n");
			}

		        $sSeq .= strtoupper($sLn);
		}

	} while (true);

	return $arrRefGenome;

}

function Warn($s) {
	global $arrWarnings;
	if (!array_key_exists($s, $arrWarnings)) {
		echo("Warning: $s\n");
		$arrWarnings[$s] = true;
	}
}

function fnGetDivergentSites($arrRefGenome , $arrOutgroupGenome, $sLookChr) {
	$arrRet = array();
	echo("Calculating divergence between ref and outgroup...\n");
	$nTotalBases = 0;
	$nDivBases = 0;
	foreach($arrRefGenome as $sChr => $sRefSeq) {
		if ($sChr != $sLookChr) continue;

		if (!array_key_exists($sChr , $arrOutgroupGenome) ) {
			die("$sChr not found in outgroup genome, the outgroup genome has to match the reference genome exactly\n");
		}

		if (strlen($arrRefGenome[$sChr]) != strlen($arrOutgroupGenome[$sChr]) ) {
			die("$sChr length differs between ref and outgroup, the outgroup genome has to match the reference genome exactly\n");
		}

		$arrRet[$sChr] = array();
		
		$nLen = strlen($arrRefGenome[$sChr]);
		for($nBase=0; $nBase < $nLen; $nBase++) {
			$sOutgroupBase = $arrOutgroupGenome[$sChr][$nBase];
			$sRefBase = $arrRefGenome[$sChr][$nBase];
			if ($sOutgroupBase == 'N' || $sRefBase == 'N' ) {
				continue;
			}

			if ($sOutgroupBase != $sRefBase) {
				$arrRet[$sChr][$nBase] = array($sOutgroupBase , $sRefBase);
				$nDivBases++;
			}

			$nTotalBases++;
		}
	}

	echo("divergence: $nDivBases / $nTotalBases = ".($nDivBases / $nTotalBases)."\n");

	return $arrRet;

}

function fnInitializeDivergentPop2Outgroup($arrDivergentRefOutgroup, $arrPops, $sChr) {
	$arrRet = array();
	foreach($arrPops as $sPop => $arrSamples) {
		$arrRet[$sPop] = $arrDivergentRefOutgroup;
	}

	return $arrRet;
}
?>
