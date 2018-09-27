<?php
require_once(dirname(__FILE__) . "/lib.php");
//output SNPs at synonymous and non-synonymous sites
//it is necessary to also provide a list of variant sites (biallelic), and the alternative bases at the site.
//the alternative bases must be given based on the orginal strand in the reference genome, not the sense strand of the gene!
//$sPop = "ORTDRY.excludehybrid";
//$sPop = "ORTWET";
//exclude CMD sites from the list

$sPop = "RACDRY.excludehybrid";
$sRefGenome = "/beegfs/group_dv/home/RCui/killifish_genomes/denovo/discovardenovo/NOR/LG/v1.2/Both_Crosses.fasta";
$sGTFFile = "/beegfs/group_dv/home/RCui/killifish_genomes/denovo/discovardenovo/NOR/LG/v1.2/NORv1.2.LG.gff";


while(count($argv) > 0) {
    $arg = array_shift($argv);
    switch($arg) {
        case '-p':
            $sPop = trim(array_shift($argv));
            break;
        case '-R':
            $sRefGenome = trim(array_shift($argv));
            break;
        case '-G':
            $sGTFFile  = trim(array_shift($argv));
            break;
	case '-o':
	    $sOutputPrefix  = trim(array_shift($argv));
	    break;

    }
}

#Edit the allelefreq_xxx folder here!
$sVariantSites = array("../bwamap/NORv1.2LG/allelefreq_above0.2_NVGout/$sPop/ancestral_polarized_for_allelefreq.invariable_der.txt" , "../bwamap/NORv1.2LG/allelefreq_above0.2_NVGout/$sPop/ancestral_polarized_for_allelefreq.txt");
$sOutputPrefix = "ret_exclCMD_pop_$sPop";


$arrVariantSites = array();
echo("Loading variant sites\n");
$nVarSitesCount = 0;
foreach($sVariantSites as $sVarFile) {
	$hVarFile = fopen($sVarFile , 'r');
	while( false !== ($sLn = fgets($hVarFile) ) ) {
		$sLn = trim($sLn);
		$arrF = explode("\t", $sLn, 5);
		if (count($arrF) < 4) continue;
		if ($arrF[0] == 'scaffold') continue;

		if (!array_key_exists($arrF[0], $arrVariantSites)) {
			$arrVariantSites[$arrF[0]] = array();
		}

		$arrVariantSites[$arrF[0]][$arrF[1]] = array($arrF[2], $arrF[3]);
		$nVarSitesCount++;
	}

	fclose($hVarFile);
}

echo("Loaded $nVarSitesCount variant sites in ". count($arrVariantSites) . " scaffolds\n"  );

/*========================================================================*/

$hRefGenome = fopen($sRefGenome, "r");
$arrRefGenome = array(); //save ref genome into the RAM
//$hSynOutTable = fopen($sOutputPrefix."_counts.txt" , "w");
$sOutputPrefix = $sOutputPrefix;
$hSynOutTable = fopen($sOutputPrefix."_syn_table.txt" , "w");
$hNonSynOutTable = fopen($sOutputPrefix."_nonsyn_table.txt" , "w");



	$sGTF = $sGTFFile;//ENSXMAG00000001050";//ENSXMAG00000001400";//ENSXMAG00000000500//ENSXMAG00000000995//ENSXMAG00000001272

	if (!file_exists($sGTF)) {
		echo("Warning: $sGTF not found\n");
		continue;
	}
	else {
		echo("Loading GTF: $sGTF ...\n");
		//die();
	}

	$oGTF = new MappedCDSGFF();//XiphoGTF();
	$oGTF->LoadGFF($sGTF, '', 'CDS'); //use "exon" features




$sLn = "";

$sSeqName = "";
$sSeq = "";
echo("Loading genome fasta...\n");
do {
	$sLn = fgets($hRefGenome);

	if ($sLn==false) {
		//fnProcessRecord($sSeqName, $sSeq);
		$arrRefGenome[$sSeqName] = $sSeq;
		break;
	}

	$sLn = trim($sLn);
	if (substr($sLn , 0,1) == ">") {
		$arrRefGenome[$sSeqName] = $sSeq;
		$sSeqName = substr($sLn , 1);
		$sSeq = "";
		
	} else {
		$sSeq .= $sLn;
	}

} while (true);




//now go over the GTF file:
$arrGenes = $oGTF->GetGenes('maker');

foreach($arrGenes as $oGene) {
	$arrRNAIDs = array_keys($oGene['mRNAs']);
	if (count($arrRNAIDs) ==0 ) {
		continue;
	}
	$arrRNA = $oGene['mRNAs'][$arrRNAIDs[0]]; // only use the first mRNA
	$bOnReverse = ($arrRNA["strand"]=="-");
	$sScfld = $oGene["scf"];
	$arrExclusionList = fnBuildExclusionList($arrRNA,  array('3primeuncertain', '5primeuncertain', 'hugeinsertionlist' )); // exclude the AA positions potentially caused by mis-annotation.

	if (count($arrExclusionList) > 0) {
		echo("Exclude following AA positions due to unreliable gene annotation $arrRNAIDs[0] :  ");
			foreach($arrExclusionList as $key => $value) {
			  echo("$key-$value ");
			}
		echo("\n");
	}

	//sort:
	if ($bOnReverse) {
		krsort($arrRNA['CDSs']);
	} else {
		ksort($arrRNA['CDSs']);
	}

	$nPrevStart = -1;
	$nPrevEnd = -1;
	$sPrevLeftOverBases = ""; //this is the left over bases at the end of an exon ,to be joined with the next exon.
	$nAccuAALen = 0; //this is the accumulative AA position on the protein.
	foreach($arrRNA['CDSs'] as $nExonStart => $nExonEnd ) { 


		//now go over each exon, note that if orientation is + , then exon starts should monotonically increase, otherwise decrease
		if ($nPrevStart != -1) {
			if ( ($bOnReverse && $nPrevStart <$nExonStart) ||  (!$bOnReverse && $nPrevStart > $nExonStart)) {
				print_r($arrRNA);
				die("Error, in the above gene the exon order is wrong!");
			}
		}



		if (!array_key_exists($sScfld , $arrRefGenome) ) {
			die("CONTIG $sScfld not found in genome fasta.\n");
		}
		$sSeqExon = substr($arrRefGenome[$sScfld],$nExonStart-1, $nExonEnd -  $nExonStart + 1);
		if ($bOnReverse ) {
			$sSeqExon = fnReverseComplement($sSeqExon);
		}

		$nBpAppendedFromPrev = strlen($sPrevLeftOverBases); // appended this many bases at the begining.
		$sSeqExon = $sPrevLeftOverBases . $sSeqExon;
		$nLeftOverBases = strlen($sSeqExon) % 3;
		$sLeftOverBases = substr($sSeqExon , strlen($sSeqExon) - $nLeftOverBases); #bug fixed on 21.7.2016
		$nInFrameBases = strlen($sSeqExon) - $nLeftOverBases;
		for($i=0;$i<$nInFrameBases;$i+=3) {
			$nAccuAALen++;

			if (fnIsExcluded($nAccuAALen , $arrExclusionList)) {
				//echo($nAccuAALen." AA excluded. $arrRNAIDs[0]\n");
				//print_r($arrExclusionList);
				//die();
				continue;
			}

			$sCodon = substr( $sSeqExon, $i, 3);

			//loop over the 3 codon bases:
			$nCodonDiff = 0;
			$sSynCache = '';
			$sNonSynCache = '';
			$sAltCodon = 'NNN';

			for ($nCodonOffset=0;$nCodonOffset < 3; $nCodonOffset++) {
				$nPos = fnGetPos($i + $nCodonOffset, strlen($sPrevLeftOverBases), $bOnReverse, $nPrevStart, $nPrevEnd , $nExonStart, $nExonEnd ); //returns 1-based coordinates
				//check if this position is mutated in the variant list

				$sAltCodon[$nCodonOffset] = $sCodon[$nCodonOffset];
				if ( (!array_key_exists($sScfld , $arrVariantSites )) || (!array_key_exists($nPos , $arrVariantSites[$sScfld] )) ) continue;

				//ok variant exists:
				$sClassified = fnClassifyVar($sCodon , $nCodonOffset, $bOnReverse, $arrVariantSites[$sScfld][$nPos] ); 
				if ($sClassified === false) {
					continue; //cannot classify
				}
				$nCodonDiff++;
				foreach($arrVariantSites[$sScfld][$nPos] as $nBase => $sBase) { 
					if ($sBase != $sAltCodon[$nCodonOffset] ) {
						$sAltCodon[$nCodonOffset] =  $sBase;
						break;
					}
				}

				if ($sClassified == 'syn') {
					 $sSynCache = "$sScfld\t".( $nPos-1 ) . "\t".$nPos.PHP_EOL ; //bed format, start pos is 0 based, end pos is 1 based
				}

				if ($sClassified == 'non-syn') {
					 $sNonSynCache = "$sScfld\t".( $nPos-1 ) . "\t".$nPos.PHP_EOL ; //bed format, start pos is 0 based, end pos is 1 based
				}
			}

			if ($nCodonDiff > 1) {
				echo("CMD:\t$arrRNAIDs[0]\t$nAccuAALen\t$sCodon\t$sAltCodon\t$nCodonDiff\t$sScfld\t$nPos\n");
			} else if ($nCodonDiff == 1) {
				if ($sSynCache!='') {
					fwrite($hSynOutTable, $sSynCache ); //bed format, start pos is 0 based, end pos is 1 based
				} 

				if ($sNonSynCache!='') {
					fwrite($hNonSynOutTable, $sNonSynCache ); //bed format, start pos is 0 based, end pos is 1 based
				} 
			}



		}

		$sPrevLeftOverBases = $sLeftOverBases;
		$nPrevStart = $nExonStart;
		$nPrevEnd = $nExonEnd;
	}
} 

function fnGetPos($i , $nPrevLeftoverBases , $bOnReverse, $nPrevStart, $nPrevEnd , $nExonStart, $nExonEnd ) {
	//check if this is the bases from the previous exon:

	if ( $i  < $nPrevLeftoverBases ) { // return index based on previous exon
		if ($bOnReverse) {
			return $nPrevStart+ ($nPrevLeftoverBases -1) - $i;
		} else {
			return $nPrevEnd - ($nPrevLeftoverBases-1) + $i;
		}
	} else {
		$i = $i - $nPrevLeftoverBases ;
		if ($bOnReverse) {
			return $nExonEnd - $i;
		} else {
			return $nExonStart + $i;
		}
	}
}

function fnBuildExclusionList($arrRNA,  $arrExcludeTypes ) {
	$nAALen = $arrRNA['CDSLen'] / 3;
	$arrRet = array();
	//print_r($arrExcludeTypes);
	//print_r($arrRNA);
	foreach($arrExcludeTypes as $sType) {
		if ( (!array_key_exists($sType , $arrRNA['annot'])) || trim($arrRNA['annot'][$sType])=='' ) {
			continue;
		}

		if ($sType == '5primeuncertain') {
			if ($arrRNA['annot'][$sType] > 0) {
				$arrRet[1] = $arrRNA['annot'][$sType];
			}
			continue;
		}
		if ($sType == '3primeuncertain') {
			if ($arrRNA['annot'][$sType] > 0) {
				$arrRet[$nAALen-$arrRNA['annot'][$sType]+1] = $nAALen;// bug fix, 2018-07-9. this bug may cause masking of whole gene.
			}
			continue;
		}

		$arrList = explode(',', $arrRNA['annot'][$sType]);
		if (count($arrList) > 0 ) {
			foreach($arrList as $sRange) {
				
				 $arrSplitRange = explode('-' , $sRange);
				if (count( $arrSplitRange)!=2) continue;
				list($nExcStart, $nExcEnd)  = $arrSplitRange;
				$arrRet[$nExcStart] = $nExcEnd;
			}
		}
		
	}
    ksort($arrRet);	
    return $arrRet;
}

function fnIsExcluded($nAccuAALen , &$arrExclusionList) {
	if (count($arrExclusionList) == 0 ) return false;

	foreach($arrExclusionList as $nExcStart => $nExcEnd) {
		if ($nExcStart <= $nAccuAALen && $nAccuAALen <= $nExcEnd) {
			return true;
		}
	}

	return false;
	
}

function fnClassifyVar($sCodon , $nCodonOffset, $bOnReverse, $arrVarBases ) {
	$sPrevAA = '';
	$bSyn = false;
	foreach($arrVarBases as $nBase => $sBase) {
		if ($bOnReverse) {
			$sBase = fnReverseComplement($sBase); //if the gene is on the reverse strand, rev the base
		}

		$sChangedCodon = substr_replace($sCodon , $sBase , $nCodonOffset);
		if (strpos($sChangedCodon , 'N') !== false || strpos($sChangedCodon , '-') !== false) {
			return false; // cannot classify due to missing data.
		}
		list($sAA) = fnDNA2Peptide($sChangedCodon, false, false);
		if ($sAA == '*') {
			return false; // cannot classify stops.
		}

		if ($sPrevAA == '') {
			$sPrevAA = $sAA;
			continue;
		}

		$bSyn = ($sPrevAA == $sAA);
	}

	if ($bSyn) {
		return 'syn';
	} else {
		return 'non-syn';
	}
}
?>
