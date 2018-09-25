<?php
require_once(dirname(__FILE__) . "/lib.php");
//allow complete masking sequence with stop codon!
$sGenomeDef = "genomes.txt";
$arrGenomes = array();
$sGTFFile = "./August6_Maker_gene_models.explicit.noseq.GENESYMBOL_15-07-2014.gff3";
$sOutputPrefix = "output/ret";
$nTotalParts = 100;
$nThisPart = 0; //index starts from 0
$nMissingCutoff = 0.5; //at most allow 50% missing
$sMissingFileRules = "missingrules_relax.txt"; // taxon on each row is OR, each line is AND. if this is not empty, the $nMissingCutoff will be ignored.
$bMaskSeqWStop = true;
$bOnlyCreateFastaIndex  = false;
while(count($argv) > 0) {
    $arg = array_shift($argv);
    switch($arg) {
        case '-R':
            $sGenomeDef = trim(array_shift($argv));
            break;
        case '-G':
            $sGTFFile  = trim(array_shift($argv));
            break;
	case '-o':
	    $sOutputPrefix  = trim(array_shift($argv));
	    break;
	case '-N':
	    $nTotalParts  = trim(array_shift($argv));
	    break;
	case '-f':
	    $nThisPart  = trim(array_shift($argv));
	    break;
	case '-m':
	    $bMaskSeqWStop  = !(trim(array_shift($argv))=="no");
	    break;
	case '-I':
	    $bOnlyCreateFastaIndex  = true;
	    break;

    }
}



/*========================================================================*/

$hGenomeDef = fopen($sGenomeDef, "r");
//$hOutTable = fopen($sOutputPrefix."_counts.txt" , "w");
$sOutputPrefix = $sOutputPrefix."_par".$nThisPart."_of_".$nTotalParts;
$hOutAA = fopen($sOutputPrefix."_AA.fasta" , "w");
$hOutDNA = fopen($sOutputPrefix."_DNA.fasta" , "w");
$hOutHyphyRelax = fopen($sOutputPrefix."_RELAX.txt" , "w");

$arrRefGenomes = array();
$arrRefFastaHack = array();
$nLn=-1;

echo("Parsing genome definitions...\n");
while( ($sLn=fgets($hGenomeDef))!==false ) {
	$sLn = trim($sLn);
	if ($sLn == "") {
		continue;
	}
	$nLn++;
	if ($nLn==0) {
		continue; // skip header
	}
	
	$arrFields = explode("\t", $sLn);
	$arrRefGenomes[$arrFields[0]] = $arrFields[1];
	$arrRefFastaHack[$arrFields[0]] = new FastaHack();
	$arrRefFastaHack[$arrFields[0]]->SetDb( $arrFields[1]);
	if ($bOnlyCreateFastaIndex) {
		$arrRefFastaHack[$arrFields[0]]->ForceCreateIndex();
	}
}


echo("Parsed in ".count($arrRefGenomes)." genome definitions\n");
if ($bOnlyCreateFastaIndex) {
	die();
}

// read missing rules
$arrMissingRules = array(); //each line is an element.
if ($sMissingFileRules!="") {
	$hMsF = fopen($sMissingFileRules , "r");
	if (!$hMsF) {
		die("Cannot open Missing rule file\n");
	}

	fgets($hMsF); //discard header
	while( false !== ($sLn = fgets($hMsF))) {
		$sLn = trim($sLn);
		if ($sLn =="") continue;
		list($nRatio, $sTaxa) = explode("\t" , $sLn); 
		$arrTaxa = explode("," , $sTaxa);
		$arrMissingRules[] = array('Ratio' => $nRatio , 'Taxa' => $arrTaxa );
	}
} else {

	foreach($arrRefGenomes as $sTaxonName => $sTaxonFile) {
		$arrMissingRules[] = array('Ratio' => $nMissingCutoff , 'Taxa' => array($sTaxonName) );
	}

}

//print_r($arrMissingRules);
//fwrite($hOutTable, "gene_id\tline\texon\ttotal_codons\tmorethanonebase_change\tidentical_codons\tsyn\tnsyn\tsyn_codon1\tsyn_codon2\tsyn_codon3\tnsyn_codon1\tnsyn_codon2\tnsyn_codon3\tmulti_codon1\tmulti_codon2\tmulti_codon3\tstopcodon_in_exon\tstop_codon_in_prev_exons".PHP_EOL);
fwrite($hOutHyphyRelax, "GeneName\tCodonLen" 
."\t"."MG94xREV"
."\t"."NULL" 
."\t"."Alternative" 
."\t"."K" 
."\t"."LR" 
."\t"."P" 
."\t"."R_omega0" 
."\t"."R_omega0_prop" 
."\t"."R_omega1" 
."\t"."R_omega1_prop" 
."\t"."R_omega2" 
."\t"."R_omega2_prop" 
."\t"."T_omega0" 
."\t"."T_omega0_prop" 
."\t"."T_omega1" 
."\t"."T_omega1_prop" 
."\t"."T_omega2" 
."\t"."T_omega2_prop" 
."\t"."Comment" 
.PHP_EOL);


foreach (glob($sGTFFile) as $sFileName) {
   if (strpos($sFileName, "_blast_results") !== false) {
		continue;
   }

    if (strpos($sFileName, "sites_matchup") !== false) {
		continue;
   }
   
   if (strpos($sFileName, "exon-direction") !== false) {
		continue;
   }
   
   if (strpos($sFileName, "exon+direction") !== false) {
		continue;
   }
   
   
	$sGTF = $sFileName;//ENSXMAG00000001050";//ENSXMAG00000001400";//ENSXMAG00000000500//ENSXMAG00000000995//ENSXMAG00000001272

	if (!file_exists($sGTF)) {
		echo("Warning: $sGTF not found\n");
		continue;
	}
	else {
		echo("Loading GTF: $sGTF ...\n");
		//die();
	}

	$oGTF = new NothosGTF();//XiphoGTF();
	$oGTF->LoadGTF($sGTF, "CDS"); //use "exon" features

	/*print_r($oGTF->arrGenes);
	die();*/
	$nGeneCounter = -1;
	foreach($oGTF->arrGenes as $sGeneName => $arrCDS) {
		$nGeneCounter++;

		if ($nGeneCounter % $nTotalParts != $nThisPart) {
			continue;
		}
		echo("Gene: ".$sGeneName.PHP_EOL);

		$sOrientation = $arrCDS[0]["orientation"];
		
		$arrFullCodingSeq = array();
		$bAllSppHaveGeneSeq = true; //all species contain sequence in this region?
		
		foreach($arrCDS as $oExon ) {
			//$sExon = $oExon["exon"] ;
			//print_r($oExon);
			//die();
			$sExon = $oExon["exon"] ;
			$nExonStart = intval($oExon["start"]);
			$nExonEnd = intval($oExon["end"]);
			
			$arrsDNA = array();
			$bOnReverseStrand = false;
			//$nSqLen = -1; //debug variable.
			foreach($arrRefFastaHack as $sSpecies => $oGenome1) {
				//echo("calling...\n");
				$sDNA1 = $oGenome1->GetCDSRegion($oExon["chr"], $nExonStart , $nExonEnd );
				if ($sDNA1 === false) { //not found in fasta file
					$bAllSppHaveGeneSeq = false;
					$arrsDNA = array();
					$arrFullCodingSeq = array();
					echo("Taxon $sSpecies doesn't have coverage in genomic region ".$oExon["chr"].":".$nExonStart."..".$nExonEnd .PHP_EOL);
					break 2;
				}
				if ( ($sOrientation=="-" && $nExonStart <$nExonEnd) || ($sOrientation=="+" && $nExonStart >$nExonEnd) ) {
					$sDNA1 = fnReverseComplement($sDNA1 );
					$bOnReverseStrand = true;
				}
				//$sDNA1 = substr($sDNA1 , $oExon["frame"]); This is not needed.
				$arrsDNA[$sSpecies] = $sDNA1;
				/*if ($nSqLen == -1) {
					$nSqLen = strlen($sDNA1);
					continue;
				} 
				if ($nSqLen!=strlen($sDNA1)) {
					echo("size inconsistent. \n");
					print_r($arrsDNA);
					die();
				}*/
			}

						
			$arrCleanedDNA = $arrsDNA ; //fnExludeCodonsWithN($arrsDNA);
			
			/*
			echo("1: ".fnDNA2Peptide($sDNA1)[0].PHP_EOL);
			echo("2: ".fnDNA2Peptide($sDNA2)[0].PHP_EOL);	
			
			echo("clean 1: ".fnDNA2Peptide($arrCleanedDNA[0])[0].PHP_EOL);
			echo("clean 2: ".fnDNA2Peptide($arrCleanedDNA[1])[0].PHP_EOL);
			*/
			
			//print_r(fnCountdNdS($arrCleanedDNA));
			
			$sExonFullID = $sGeneName."_exon".$sOrientation."direction_exon_".$sExon."_";
					
			
			foreach($arrCleanedDNA as $sTaxon => $sSeq) {

				fwrite($hOutDNA , ">".$sExonFullID."_".$sTaxon.PHP_EOL.$sSeq.PHP_EOL );
				if (!array_key_exists($sTaxon,$arrFullCodingSeq )) {
					$arrFullCodingSeq[$sTaxon] = "";
				}
				
				$arrFullCodingSeq[$sTaxon] .= $sSeq;
				
				
			}
			
		}

		//continue;//debug.
		if (!$bAllSppHaveGeneSeq) { //if alignment cannot be made completely due to missing taxon
			fwrite($hOutHyphyRelax, "$sGeneName\tSome taxa miss sequence in this region.".PHP_EOL);
			continue;
		}
		



		$arrFullCodingSeq = fnExcludeLastStop($arrFullCodingSeq);
		$arrAlnRet = fnTranslateAlignment($arrFullCodingSeq , false); //exclude last stop if exists.

			
		$sExonFullID = $sGeneName."_exon".$sOrientation;		
		foreach($arrAlnRet["proteins"] as $sTaxon => $sSeq) {
				fwrite($hOutAA , ">".$sExonFullID."_".$sTaxon.PHP_EOL.$sSeq.PHP_EOL );
				fwrite($hOutDNA , ">".$sExonFullID."_".$sTaxon.PHP_EOL.$arrFullCodingSeq[$sTaxon].PHP_EOL );
		}


		//Check if all sequences are identical, if so then no need to continue with codeml to save time
		$bDiff = false;
		$sTempSeq = -1;
		foreach($arrFullCodingSeq as $sTaxon => $sSeq) {
			if ($sTempSeq == -1) {
				$sTempSeq = $sSeq;
				continue;
			}
			if ($sTempSeq != $sSeq) {
				$bDiff = true;
				break;
			}
		}

		if (!$bDiff) {
			fwrite($hOutHyphyRelax, "$sGeneName\tAll sequences identical.".PHP_EOL);
			echo("All sequences identical \n");
			continue;
		}
			
		

		//exclude all sequences with a stop codon
		$sComment = "-";
		if ( $arrAlnRet["stopcodon"]  && $bMaskSeqWStop) {
			$sComment = "TaxaMaskedDueToStopCodon:";
			foreach( $arrAlnRet["stoptaxa"] as $sTaxonWStop) {
				$arrFullCodingSeq[$sTaxonWStop] = str_repeat("N", strlen($arrFullCodingSeq[$sTaxonWStop]));
			}
			$sComment .= implode(",", $arrAlnRet["stoptaxa"]);
			$arrAlnRet["stopcodon"] = false;
		}


		//check amount of missing data, if too much, don't continue.

		$arrMissingCounts = array();
		foreach($arrFullCodingSeq  as $sTaxon => $sSeq) {
			$arrMissingCounts[$sTaxon] = 1-substr_count($sSeq , "N") / strlen($sSeq) ;
		}

		foreach($arrMissingRules as $oRule) { //check each rule
			$arrTaxaToCheck = $oRule['Taxa'];
			$bPass = false;

			foreach($arrTaxaToCheck as $sTaxonToCheck) {
				$bPass = ($bPass || ($oRule['Ratio'] <= $arrMissingCounts[$sTaxonToCheck]));
			}

			if (!$bPass) {
				fwrite($hOutHyphyRelax, "$sGeneName\tAll of these taxa have extensive missing data in this region: ".implode(",",$arrTaxaToCheck).PHP_EOL);
				continue 2;
			}
		}

		

	}
	
}

function fnGetUniqueHits($arrMapInfo) {
	$arrRet = array();
	foreach($arrMapInfo as $oRecord) {
		if (!fnHitExists($arrRet,$oRecord )) {
			array_push($arrRet , $oRecord);
		}
	}
	
	return $arrRet;
}

function fnHitExists($arrRet , $arrCompare) {
	foreach($arrRet as $oRecord) {
		if (implode("%",$oRecord) == implode("%",$arrCompare)) {
			return true;
		}
	}
	return false;
}

?>
