<?php
require_once(dirname(__FILE__) . "/lib.php");
//output SNPs at 4-fold degenerate sites
//allow complete masking sequence with stop codon!
$sGenomeDef = "genomes.txt";
$arrGenomes = array();
$sGTFFile = "./August6_Maker_gene_models.explicit.noseq.GENESYMBOL_15-07-2014.gff3";
$sPremadeFasta = "../ret_DNA.fasta"; //premade fasta alignment
$sOutputPrefix = "ret";
$nMissingCutoff = 1; //at most allow 50% missing for the gene per taxon
$nMissingTaxaCutoff = 0.6; //at most 60% of the taxa missing
$bAllowN = true; //allows N
$sMissingFileRules = "missingrules_relax.txt"; // taxon on each row is OR, each line is AND. if this is not empty, the $nMissingCutoff will be ignored.
$bMaskSeqWStop = true;
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

    }
}



/*========================================================================*/

$hGenomeDef = fopen($sGenomeDef, "r");
//$h4folddegenOutTable = fopen($sOutputPrefix."_counts.txt" , "w");
$sOutputPrefix = $sOutputPrefix;

$h12codonOutTable = fopen($sOutputPrefix."_allcodons_table.txt" , "w");

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
}

echo("Parsed in ".count($arrRefGenomes)." genome definitions\n");
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

	$sGTF = $sGTFFile;//ENSXMAG00000001050";//ENSXMAG00000001400";//ENSXMAG00000000500//ENSXMAG00000000995//ENSXMAG00000001272

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


//print_r($arrMissingRules);
//fwrite($h4folddegenOutTable, "gene_id\tline\texon\ttotal_codons\tmorethanonebase_change\tidentical_codons\tsyn\tnsyn\tsyn_codon1\tsyn_codon2\tsyn_codon3\tnsyn_codon1\tnsyn_codon2\tnsyn_codon3\tmulti_codon1\tmulti_codon2\tmulti_codon3\tstopcodon_in_exon\tstop_codon_in_prev_exons".PHP_EOL);

fwrite($h12codonOutTable, "GeneName\tScfld" 
."\t"."Pos"
."\t".implode("\t",array_keys($arrRefGenomes))
.PHP_EOL);



//Now parse through the premade fasta file
$hPremadeFasta = fopen($sPremadeFasta , "r");
$sLn = "";

$sSeqName = "";
$sSeq = "";
do {
	$sLn = fgets($hPremadeFasta);

	if ($sLn==false) {
		fnProcessRecord($sSeqName, $sSeq);
		break;
	}

	$sLn = trim($sLn);
	if (substr($sLn , 0,1) == ">") {
		fnProcessRecord($sSeqName, $sSeq);
		$sSeqName = substr($sLn , 1);
		$sSeq = "";
		
	} else {
		$sSeq .= $sLn;
	}

} while (true);
fnProcessAln($sCurrentGene,$arrFullCodingSeq); //don't forget the last one

$arrFullCodingSeq = array();
$sCurrentGene = "";
function fnProcessRecord($sSeqName, $sSeq) {
	global $arrFullCodingSeq , $sCurrentGene, $arrRefGenomes;
	if (strpos($sSeqName , "direction") !== false  ) {
		return; //don't process fragments
	}

	preg_match("/(\S+)_exon[+-]_+(\S*)/", $sSeqName, $arrParsed); //bug fixed 2016-8-2

	if (count($arrParsed)!=3 ) {
		return;
	}
	$sGeneName = $arrParsed[1];
	$sSpecies = $arrParsed[2];

	if ($sCurrentGene == $sGeneName ) {
		if (!array_key_exists($sSpecies , $arrFullCodingSeq)) {
			$arrFullCodingSeq[$sSpecies] = "";
		}

		$arrFullCodingSeq[$sSpecies] .= $sSeq;
	} else {
		if (count($arrFullCodingSeq) > 0) {
			fnProcessAln($sCurrentGene,$arrFullCodingSeq);
		}
		$arrFullCodingSeq = array();
		$sCurrentGene = $sGeneName;
		$arrFullCodingSeq[$sSpecies] = $sSeq;
	}

}


function fnProcessAln($sGeneName,$arrFullCodingSeq) {
	global   $sGTFFile, $arrMissingRules, $oGTF, $h4folddegenOutTable, $h12codonOutTable, $bAllowN , $nMissingTaxaCutoff;
   



		$sScfld = "unknown";
		$arrFullCodingSeq = fnExcludeLastStop($arrFullCodingSeq);
		$arrAlnRet = fnTranslateAlignment($arrFullCodingSeq , false); //exclude last stop if exists.

	

		//check amount of missing data, if too much, don't continue.

		$arrMissingCounts = array();
		$nSeqLen = -1;
		foreach($arrFullCodingSeq  as $sTaxon => $sSeq) {
			$arrMissingCounts[$sTaxon] = 1-substr_count($sSeq , "N") / strlen($sSeq) ;
			if ($nSeqLen == -1) {
				$nSeqLen = strlen($sSeq) ;
			}
		}

		foreach($arrMissingRules as $oRule) { //check each rule
			$arrTaxaToCheck = $oRule['Taxa'];
			$bPass = false;

			foreach($arrTaxaToCheck as $sTaxonToCheck) {
				$bPass = ($bPass || ($oRule['Ratio'] <= $arrMissingCounts[$sTaxonToCheck]));
			}

			if (!$bPass) {
				//fwrite($h4folddegenOutTable, "$sGeneName\tAll of these taxa have extensive missing data in this region: ".implode(",",$arrTaxaToCheck).PHP_EOL);
				return;
			}
		}



		if (!$arrAlnRet["stopcodon"] ) { //there is no stop codon, export SNP
			
			echo("Exporting SNPs seqlen=$nSeqLen... name=$sGeneName\n");
		
			for($nPos=0;($nPos+2)<$nSeqLen;$nPos+=3) {


				$arrFoundSNPs1 = array(); //1st and 2nd positions
				$arrFoundSNPs2 = array(); //1st and 2nd positions
				$arrFoundSNPs3 = array(); //1st and 2nd positions
				$bWriteThis = true;
				foreach($arrFullCodingSeq  as $sTaxon => $sSeq) {
					$sCodon = substr($sSeq  , $nPos, 3);
					if (strpos($sCodon , 'N') !== false && (!$bAllowN ) ) { //if anything is N, do not write this position.
						$bWriteThis = false;
						break;
					}
					$arrFoundSNPs1[] = substr($sCodon, 0, 1); //put the 3rd position
					$arrFoundSNPs2[] = substr($sCodon, 1, 1); //put the 3rd position
					$arrFoundSNPs3[] = substr($sCodon, 2, 1); //put the 3rd position

				}

				//print_r($arrFoundSNPs);
				if ( $bWriteThis  ) { 

					if (fnPercentN($arrFoundSNPs1) <= $nMissingTaxaCutoff) {
						fwrite($h12codonOutTable, "$sGeneName\t"
						.implode("\t", fnGenePosToGenomePos($sGeneName, $nPos+1))
						."\t".implode("\t" , $arrFoundSNPs1)
						.PHP_EOL);	
					}


					if (fnPercentN($arrFoundSNPs2) <= $nMissingTaxaCutoff) {
						fwrite($h12codonOutTable, "$sGeneName\t"
						.implode("\t", fnGenePosToGenomePos($sGeneName, $nPos+2))
						."\t".implode("\t" , $arrFoundSNPs2)
						.PHP_EOL);
					}	
					if (fnPercentN($arrFoundSNPs3) <= $nMissingTaxaCutoff) {
						fwrite($h12codonOutTable, "$sGeneName\t"
						.implode("\t", fnGenePosToGenomePos($sGeneName, $nPos+3))
						."\t".implode("\t" , $arrFoundSNPs3)
						.PHP_EOL);
					}		

				}







			}
			
		}
		else {
			//fwrite($h4folddegenOutTable, "$sGeneName\tStopCodonFoundIn:".implode(",",$arrAlnRet["stoptaxa"]).PHP_EOL);
		}
		


	
}

function array_count_values_of($value, $array) {
    $counts = array_count_values($array);
    return array_key_exists($value, $counts)? $counts[$value] : 0;
}

function fnPercentN($arr) {

	return array_count_values_of("N", $arr) / count($arr);
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

function fnGenePosToGenomePos($sGeneName, $nPos) { //$nPos is 1 based
	global $oGTF;

	if (!array_key_exists($sGeneName , $oGTF->arrGenes )) {
		return array("unknown" , "0");
	}

	$arrCDS = $oGTF->arrGenes[$sGeneName];

	$nAccuLen = 0;

	$sOrientation = $arrCDS[0]["orientation"];
		foreach($arrCDS as $oExon ) {
			//$sExon = $oExon["exon"] ;
			//print_r($oExon);
			//die();
			$sExon = $oExon["exon"] ;
			$nExonStart = intval($oExon["start"]);
			$nExonEnd = intval($oExon["end"]);
			$sScfld = $oExon["chr"];
			$bOnReverseStrand = false;

			

			if ( ($sOrientation=="-" && $nExonStart <$nExonEnd) || ($sOrientation=="+" && $nExonStart >$nExonEnd) ) {
				$bOnReverseStrand = true;
			}

			$nAccuLen +=  abs($nExonStart - $nExonEnd)+1;// - $oExon["frame"]; don't bother with frame

			if ($nPos <= $nAccuLen) {
				if (!$bOnReverseStrand) {
					$nTruePos = $nExonEnd - ( $nAccuLen - $nPos);
				} else {
					$nTruePos = $nExonStart + ( $nAccuLen - $nPos);
				}

				return(array($sScfld , $nTruePos));
			} 
			
		}


				print_r($arrCDS);
				echo("$sGeneName, $nPos".PHP_EOL);

	die("The above gene's coordinate cannot be mapped to genome coordinate according to the GTF file. Check upstream scripts!\n");
	return array("unknown" , "0");

	
}

?>
