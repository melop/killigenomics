<?php
ini_set('memory_limit', -1);

require_once(dirname(__FILE__) . "/lib.php");
//allow complete masking sequence with stop codon!
$sRefDef = "refs.txt"; // where are the reference genomes used for mapping and where are their associated GFFs?
//$sOrthologList = "orthologs.txt"; // the gene names from each reference genome that are identified as orthologs. Alignments will only be produced for genes in this list.
$sOrthologList = "ortholog.improved.txt";
$sOutgroupList = "outgroups.txt"; //use these sequences directly, instead of exporting from fasta/gff pairs
$sGenomeDef = "genomes.txt";
$arrGenomes = array();
$sOutputPrefix = "output/ret";
$nTotalParts = 1;
$nThisPart = 0; //index starts from 0
$nMissingCutoff = 0.5; //at most allow 50% missing
$sMissingFileRules = "missingrules_relax.txt"; // taxon on each row is OR, each line is AND. if this is not empty, the $nMissingCutoff will be ignored.
$bMaskSeqWStop = true;
$bOnlyCreateFastaIndex  = false;
$nCovBlackListPCutoff = 0.01;

exec("mkdir -p output");

while(count($argv) > 0) {
    $arg = array_shift($argv);
    switch($arg) {
    case '-r':
        $sRefDef = trim(array_shift($argv));
        break;
    case '-R':
        $sGenomeDef = trim(array_shift($argv));
        break;
	case '-O':
	    $sOrthologList  = trim(array_shift($argv));
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
	case '-C':
	    $nCovBlackListPCutoff = floatval(trim(array_shift($argv)));
	    break;

    }
}

/*============ load outgroup fasta =============*/
$hOutGroupList = fopen($sOutgroupList, 'r');
$arrOutgroupProt = array();

while(false !== ($sLn = fgets($hOutGroupList) ) ) {
    $sLn = trim($sLn);
    if ($sLn == '') continue;
    if ($sLn[0] == '#') continue;
    list($sTaxon, $sFile) = explode("\t" , $sLn);
    echo("loading CDS file for $sTaxon...\n");
    $arrFasta = fnParseProtFasta($sFile);
    $arrOutgroupProt[$sTaxon] = $arrFasta;
}

//print_r($arrOutgroupProt);
//die();
/*============ Load the ref genomes ======================================*/
$hRefDef = fopen($sRefDef , "r");
$arrRefGTFs = array();

$nLn = -1;
echo("Parsing reference genome definitions...\n");
while( ($sLn=fgets($hRefDef))!==false ) {
	$sLn = trim($sLn);
	if ($sLn == "") {
		continue;
	}
	$nLn++;
	if ($nLn==0) {
		continue; // skip header
	}
	
	$arrFields = explode("\t", $sLn);
	$arrRefGTFs[$arrFields[0]] = $arrFields[2];
}




/*========================================================================*/

$hGenomeDef = fopen($sGenomeDef, "r");
//$hOutTable = fopen($sOutputPrefix."_counts.txt" , "w");
$sOutputPrefix = $sOutputPrefix."_par".$nThisPart."_of_".$nTotalParts;
$hOutAA = fopen($sOutputPrefix."_AA.fasta" , "w");
$hCleanOutAA = fopen($sOutputPrefix."_cleanAA.fasta" , "w");
$hOutDNA = fopen($sOutputPrefix."_DNA.fasta" , "w");
$hCleanOutDNA = fopen($sOutputPrefix."_cleanDNA.fasta" , "w");


$arrRefGenomes = array();
$arrRefFastaHack = array();
$nLn=-1;
$arrGTFs = array();
$arrPseudoGenomes = array();
$arrSp2Ref = array();
$arrSp2Genome = array();
$arrSpCovBlackList = array();

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
	
	list($sRef, $sSp, $sPseudoGenomeFa, $sCovBlkList ) = explode("\t", $sLn);
	if (!array_key_exists($sRef, $arrRefGTFs)) {
		die("Reference $sRef not found in reference definition file $sRefDef\n");
	}

	if (!array_key_exists( $sRef , $arrPseudoGenomes) ) {
		$arrPseudoGenomes[$sRef] = array();
	}

	$arrPseudoGenomes[$sRef][] = $sPseudoGenomeFa;
	$arrSp2Ref[$sSp] =  $sRef;
	$arrSp2Genome[$sSp] =  $sPseudoGenomeFa;
	$arrSpCovBlackList[$sSp] = fnParseBlackList($sCovBlkList);
}

//load gff
foreach($arrPseudoGenomes as $sRef => $arrPseudoGenomesForRef) {
	$arrGTFs[$sRef] = new MappedCDSGFF();
	$arrGTFs[$sRef]->LoadGFF($arrRefGTFs[$sRef] , $arrPseudoGenomesForRef , 'CDS' ); //multiple genome mode
}

echo("Parsed in ".count($arrGTFs)." genome definitions\n");

//print_r($arrRefGenomes);
//print_r($arrRefFastaHack);
//die();
if ($bOnlyCreateFastaIndex) {
	die();
}

/*============ Load the ortholog definitions  ======================================*/

$arrOrthologSpp = array();
$arrOrthologs = array();
$arrOrthologMeta = array();


$nLn = -1;
echo("Parsing ortholog definitions $sOrthologList ...\n");
$hOrthologList = fopen($sOrthologList , "r");

while( ($sLn=fgets($hOrthologList))!==false ) {
	$sLn = trim($sLn);
	if ($sLn == "") {
		continue;
	}
	$nLn++;
	if ($nLn==0) {
		$arrOrthologSpp = array_flip(array_slice( explode("\t" , $sLn) , 3));
		
		continue; // skip header
	}
	
	$arrFields = explode("\t", $sLn);
	$arrOrthologs[] =  array_slice($arrFields , 3);
	$arrOrthologMeta[] = array_slice($arrFields , 0, 3);
}

echo("Loaded ". count($arrOrthologs) ." ortholog definitions\n" );
//die();

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

//now iterate through the ortholog definition
$arrRefGenomeSpp = array_keys($arrRefGTFs);
for($nOrth=0;$nOrth<count($arrOrthologs);$nOrth++ ) {

	if ($nOrth % $nTotalParts != $nThisPart) {
		continue;
	}

	$arrOrthMeta = $arrOrthologMeta[$nOrth];
	$arrSpGeneIDs = array(); 
    
    foreach($arrRefGenomeSpp as $sRefGenomeSp) {
        if (!array_key_exists($sRefGenomeSp , $arrOrthologSpp) ) {
            echo("Warning: reference species $sRefGenomeSp not found in ortholog definition. This whole gene will be skipped\n");
	    continue 2;
        }
        
        $arrSpGeneIDs[$sRefGenomeSp] = $arrOrthologs[$nOrth][$arrOrthologSpp[$sRefGenomeSp]];
    }
    
	echo("Gene: ");
	// now check if these gene names exist in these GFF files:
	foreach($arrSpGeneIDs as $sRef => $sIDString ) {
	    list($s1, $sSpGeneID, $s3) = explode('|', $sIDString);
		if (  $arrGTFs[$sRef]->fnSpGeneId2RNAID($sSpGeneID) === false) {
			echo("The gene name $sSpGeneID is not found in the GFF file for reference genome $sRef. ID string that caused error: $sIDString. Gene will be skipped \n");
			continue 2;
		}
		echo("  $sRef:$sSpGeneID ");
	}
	echo("\n");

	$arrFullCodingSeq = array(); //array_combine(array_keys($arrSp2Ref) , array_fill(0, count($arrSp2Ref) , "" )  );
	$bAllSppHaveGeneSeq = true;
	foreach($arrSp2Ref as $sSpecies => $sRefForSp) {
		//print_r(array_keys($arrGTFs));
		//print_r($arrSp2Ref);

		if ( array_key_exists( $arrSpGeneIDs[$sRefForSp], $arrSpCovBlackList[$sSpecies])) {
			echo("Notice: $sSpecies coverage too high for gene ".$arrSpGeneIDs[$sRefForSp]." ".implode(" ", $arrOrthMeta)." , exclude\n");
			continue;
		}

		$oGTF = $arrGTFs[$sRefForSp];
		list($sTmp, $sSpGeneID , $sTmp2) =  explode('|', $arrSpGeneIDs[$sRefForSp]);
		$sRNAID = $oGTF->fnSpGeneId2RNAID($sSpGeneID);
		echo("$sSpecies / $sRefForSp / $sRNAID / " . $sSpGeneID ."\n");
		$omRNA = $oGTF->GetmRNA('maker' , $sRNAID);
		$oExtracted = $oGTF->ExtractmRNASequence('maker' , $sRNAID, $arrSp2Genome[$sSpecies]);
		$sTrimAA = $oExtracted['AA'];
		$sTrimDNA = $oExtracted['mRNA'];

		//mask huge gaps:
		if (array_key_exists('hugeinsertionlist', $omRNA['annot'])) {
			$arrHugeInsList = explode(",", $omRNA['annot']['hugeinsertionlist']);
			foreach($arrHugeInsList as $sHugeGap) {
				$sHugeGap = trim($sHugeGap);
				if ($sHugeGap == "") continue;
				list($nGapStart, $nGapEnd) = explode('-', $sHugeGap);
				$nGapLen = abs($nGapEnd-$nGapStart) + 1;
				$nDNAGapLen = $nGapLen * 3;
				$sTrimAA = substr_replace($sTrimAA , str_repeat('X', $nGapLen) , ($nGapStart-1) ,  $nGapLen);
				$sTrimDNA = substr_replace($sTrimDNA , str_repeat('N', $nDNAGapLen) , ($nGapStart-1)*3 ,  $nDNAGapLen);
			}
		}

		$n5PrimeUncertain = 0;
		$n3PrimeUncertain = 0;
		if (array_key_exists('5primeuncertain', $omRNA['annot'])) {
			$n5PrimeUncertain = intval($omRNA['annot']['5primeuncertain']);
		}

		if (array_key_exists('3primeuncertain', $omRNA['annot'])) {
			$n3PrimeUncertain = intval($omRNA['annot']['3primeuncertain']);
		}

		//trim off 3'
		$sTrimAA = substr($sTrimAA , 0 , strlen($sTrimAA) - $n3PrimeUncertain );
		$sTrimDNA = substr($sTrimDNA , 0 , strlen($sTrimDNA) - ($n3PrimeUncertain*3) );
		//trim off 5'
		$sTrimAA = substr($sTrimAA , $n5PrimeUncertain );
		$sTrimDNA = substr($sTrimDNA ,  $n5PrimeUncertain * 3 );

		$arrFullCodingSeq[$sSpecies] = $sTrimDNA;
	}


	/*if (!$bAllSppHaveGeneSeq) { //if alignment cannot be made completely due to missing taxon
		fwrite($hOutHyphyRelax, "$sSpGeneID\tSome taxa miss sequence in this region.".PHP_EOL);
		continue;
	}*/
		

	$arrFullCodingSeq = fnExcludeLastStop($arrFullCodingSeq);
	$arrAlnRet = fnTranslateAlignment($arrFullCodingSeq , false); //exclude last stop if exists.

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

	$arrNonNRatio = array();
	foreach($arrFullCodingSeq  as $sTaxon => $sSeq) {
		$arrNonNRatio[$sTaxon] = 1-substr_count($sSeq , "N") / strlen($sSeq) ;
	}


	//Added - now append the outgroup sequences
	foreach( array_keys($arrOutgroupProt) as $sOutgroupSp ) {
        //find it in the ortholog table
		if (!array_key_exists($sOutgroupSp, $arrOrthologSpp) ) {
			die("Outgroup species $sOutgroupSp is not defined in the ortholog table\n");
		}
		$sIDString = $arrOrthologs[$nOrth][$arrOrthologSpp[$sOutgroupSp]];
		list($s1, $s2, $sGeneEnsNCBIID) = explode('|', $sIDString);
		if (substr($sGeneEnsNCBIID, 0, 2) == 'EN' ) {
			list($sGeneEnsNCBIID, $s1) = explode('_', $sGeneEnsNCBIID);
		}
        
		if (!array_key_exists($sGeneEnsNCBIID ,$arrOutgroupProt[$sOutgroupSp] ) ) {
            		echo("Warning: gene $sGeneEnsNCBIID not found in cds file for $sOutgroupSp, skipped\n");
			continue;
		}
        
		$sSeq2Add = $arrOutgroupProt[$sOutgroupSp][$sGeneEnsNCBIID];
		$sLastCodon = substr($sSeq2Add , -3);
		if ($sLastCodon == "TGA" || $sLastCodon == "TAG" || $sLastCodon == "TAA" || $sLastCodon == "TAR" ) {
			$sSeq2Add = substr($sSeq2Add, 0, strlen($sSeq2Add) -3 );
		}

	        $arrFullCodingSeq[$sOutgroupSp] = $sSeq2Add ;  //add the outgroup spp.
	}



	$arrFullCodingSeq2 = array();
	foreach($arrFullCodingSeq as $sTaxon => $sSeq) {

		if (array_key_exists($sTaxon , $arrNonNRatio) && $arrNonNRatio[$sTaxon] < $nMissingCutoff) continue;
		$sRefForSp = array_key_exists($sTaxon,  $arrSp2Ref)? $arrSp2Ref[$sTaxon] : $sTaxon; //this is the reference for the species
		$sSpGeneId = ($sRefForSp ==  $sTaxon)? "000000" : (explode('|', $arrSpGeneIDs[$sRefForSp]))[1];

		$sSeqID = "Ortho:$nOrth;Sp:$sTaxon;MappedToRef:$sRefForSp;SpGeneId:$sSpGeneId;GroupId:".$arrOrthMeta[0];//.";GeneSymbol:".$arrOrthMeta[1];

		$arrFullCodingSeq2[$sSeqID] = $sSeq;

	}

	$arrTrX = fnTranslatorX($arrFullCodingSeq2);
	if (false !== $arrTrX) {

		fwrite($hOutDNA , $arrTrX[0]);
		fwrite($hCleanOutDNA , $arrTrX[1]);
		fwrite($hOutAA , $arrTrX[2]);
		fwrite($hCleanOutAA , $arrTrX[3]);
	} else {
		echo("TranslatorX failed on gene \n");
	}

}

function fnParseBlackList($sCovBlkList) {
	global $nCovBlackListPCutoff;

	if (!file_exists($sCovBlkList)) {
		echo("Warning: $sCovBlkList not found, no black list will be applied.\n");
		return array();
	}

	$hCovBlkList = fopen($sCovBlkList , 'r');
	$arrRet = array();
	while(false !== ($sLn = fgets($hCovBlkList) ) ) {
		$sLn = trim($sLn);
		$arrF = explode("\t" , $sLn);
		if (count($arrF) <9 ) continue;

		if ($arrF[8] < $nCovBlackListPCutoff) {
			$arrRet[$arrF[0]] = true;
		}
	}

	return $arrRet;
}

function fnParseProtFasta($sFile) {
    
    $hIn = fopen($sFile, 'r');
    $sSeq = "";
    $sSeqName = "";
    $arrRet = array();
    do {
        $sLn = fgets($hIn);
        if (false === $sLn || $sLn[0] == ">") {
            if ($sSeq != '') {
                $arrRet[$sSeqName] = $sSeq;
            }
            if (false === $sLn) break;
            $sSeqName = substr(trim($sLn), 1);
            $sSeq = "";
            continue;
        }
        
        $sSeq .= trim($sLn);
        
    } while(true);
    
    return $arrRet;
}
?>
