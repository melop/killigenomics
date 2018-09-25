<?php
//ini_set('memory_limit', '20000M');
var_dump(ini_get("memory_limit"));
//die();
	require_once(dirname(__FILE__) . "/config.php");
	require_once(dirname(__FILE__) . "/lib.php");
	require_once(dirname(__FILE__) . "/codons.php");
/*
2015-4-1
changed for yumi:
allow variant-only bcf format where unknown sites are treated as identical as the reference.
*/
//parse commandline

$nMaxGapSize = 1000;
$nShortestAln = 10;
$nMinDepth = 10;
$nSparseCutoff = 0.5; // if less than 50% of taxa contains character then that column is defined as "sparse".
$sOutFolder = "";
$sConfigFile = "";
$sStartFromContig = "";
$bPenalizeCoverage = false;
$nMinQual = 20;
$nAlternativeReadRatio = 1;
$sCharRequirementFile = "";
$bNAsDifference = true;
$arrCharRequirements = array();
$sTranscriptListFile = "";
$arrTranscriptList = array();
$bDivergenceFilter = true;
$bIncludeRefNPos = false;
$bCodeIndelOrMultiNonRefAlleleAsN = false;
$bTreatNAsRefBase = false;

while(count($argv) > 0) {
    $arg = array_shift($argv);
    switch($arg) {
        case '-w':
            $nMaxGapSize   = intval(array_shift($argv));
            break;
        case '-c':
            $sConfigFile = array_shift($argv);
            break;
        case '-o':
            $sOutFolder  = array_shift($argv);
            break;
	case '-s':
	    $nSparseCutoff = floatval(array_shift($argv));
	    break;
	case '-q':
	    $nMinQual = floatval(array_shift($argv));
	    break;
	case '-r':
	    $nAlternativeReadRatio = floatval(array_shift($argv));
	    break;
	case '-l':
	    $nShortestAln = intval(array_shift($argv));
	    break;
	case '-u':
	    $sStartFromContig = array_shift($argv);
	    break;
	case '-d':
	    $nMinDepth=intval(array_shift($argv));
	    break;
	case '-x':
	    $bPenalizeCoverage = (trim(array_shift($argv))=="yes");
	    break;
	case '-R':
	    $sCharRequirementFile = trim(array_shift($argv));
	    break;
	case '-T':
	    $sTranscriptListFile = trim(array_shift($argv));
	    break;
	case '-N':
	    $bNAsDifference = (trim(array_shift($argv))=="yes");
	    break;
	case '-D':
            $bDivergenceFilter = (trim(array_shift($argv))=="yes");
	    break;
	case '-n':
            $bIncludeRefNPos = (trim(array_shift($argv))=="yes");
	    break;
	case '-b':
            $bTreatNAsRefBase = (trim(array_shift($argv))=="yes");
	    break;
	case '-I':
            $bCodeIndelOrMultiNonRefAlleleAsN = (trim(array_shift($argv))=="yes");
	    break;

    }
}

//echo all the parameters
            echo("\tMaxGapSize: $nMaxGapSize".PHP_EOL);
            echo("\tConfigFile: $sConfigFile".PHP_EOL);
            echo("\tOutput folder: $sOutFolder".PHP_EOL);
            echo("\tSparseCutoff: $nSparseCutoff".PHP_EOL);
            echo("\tMin Quality: $nMinQual".PHP_EOL);
            echo("\tShortest alignment to output: $nShortestAln bp".PHP_EOL);
            echo("\tStart from contig: $sStartFromContig".PHP_EOL);
            echo("\tMinDepth: $nMinDepth".PHP_EOL);
            echo("\tPenalize Coverage: ".fnB2S($bPenalizeCoverage).PHP_EOL);
            echo("\tCharRequirementFile:  $sCharRequirementFile".PHP_EOL);
            echo("\tTranscriptListFile:  $sTranscriptListFile ".PHP_EOL);
            echo("\tCount 'N's as difference for divergent filter: ".fnB2S($bNAsDifference).PHP_EOL);
            echo("\tEnable divergent filter: ".fnB2S($bDivergenceFilter).PHP_EOL);
            echo("\tInclude positions where reference is N: ".fnB2S($bIncludeRefNPos).PHP_EOL);
	    echo("\tCode indel or more than 1 non-ref allele as N (otherwise position deleted): ".fnB2S($bCodeIndelOrMultiNonRefAlleleAsN).PHP_EOL);
	    echo("\tRatio of reads supporting the alternative allele (-r), \n\tfor example, 1 means 0 reads supporting ref,\n\t 1 supporting alt, or 0 to 2, 0 to 3 etc.\n\t if doesn't pass this, masked as N:".fnB2S($nAlternativeReadRatio).PHP_EOL);

	if ($nMaxGapSize <= 100 ) {
		echo("Missing window size. Specify by -w \n");
		Usage();
	}

	if ($sOutFolder == "" ) {
		echo("Missing output folder. Specify by -o \n");
		Usage();
	}

	if ($sConfigFile == "" ) {
		echo("Missing configuration file. Specify by -c \n");
		Usage();
	}

	if ($nSparseCutoff > 1 || $nSparseCutoff < 0) {
		echo("Sparse site cutoff should be [0,1] \n");
		Usage();
	}

	if ( $nShortestAln < 0) {
		echo("Shortest alignment to output must be >=0 \n");
		Usage();
	}

	//parse config file

	$hConfig = fopen($sConfigFile, "r");

	if (!$hConfig ) {
		echo("Could not open configuration file $sConfigFile \n");
		Usage();
	}

	$sRefGenome = "";
	$sRefTaxon = "";

	$nLnCount = 0;
	$arrBCFs = array();
	$arrCoverages = array();

	while ($sLn = fgets($hConfig)) {

		$sLn = trim($sLn);
		if ($sLn =="") continue;
		if (!$bPenalizeCoverage)
		{
    			list ($sTaxon, $sFile) = explode("\t",$sLn);
    			//... do something with the values
		}
		else {
			list ($sTaxon, $sFile, $nCoverage) = explode("\t",$sLn);
		}

		if ($nLnCount == 0) {
			$sRefTaxon = trim($sTaxon);
			$sRefGenome = trim($sFile);
			
		}
		else {
			$arrBCFs[trim($sTaxon)] = trim($sFile);
			if ($bPenalizeCoverage) {
				$arrCoverages[trim($sTaxon)] = $nCoverage;
			}
		}

		$nLnCount++;
	}
	fclose($hConfig);

	//parse character requirement file:

	if ($sCharRequirementFile !=="") {

		$hCharRqFile = fopen($sCharRequirementFile , "r");
		while ($sLn = fgets($hCharRqFile)) {
			$arrLn = explode("\t" , trim($sLn));
			array_push($arrCharRequirements , $arrLn); 
		}
		fclose($hCharRqFile);

	}


	if ($sTranscriptListFile!=="") {

		$hTListFile = fopen($sTranscriptListFile , "r");
		while ($sLn = fgets($hTListFile)) {
			$arrTranscriptList[trim($sLn)] = 1;
		}
		fclose($hTListFile);

	}

	/*echo(count($arrTranscriptList));
	die();*/


	//check config:

	if ($sRefTaxon == "" || $sRefGenome == "") {
		echo("Reference taxon and genome not present\n");
		Usage();
	}

	if (count($arrBCFs) < 1) {
		echo("No other taxon present in the configuration file other than the reference.\n");
		print_r($arrBCFs);
		Usage();
	}

	if (count($arrBCFs) < 3) {
		echo("Warning: You need at least 4 taxa (including the reference) to do AU test.\n");
		
	}

	$arrBCFUtils = array();

	foreach($arrBCFs as $sTaxon => $sPath) {

		$arrBCFUtils[$sTaxon] = new Bcftools();
		$arrBCFUtils[$sTaxon]->SetBcf($sPath); //Load the Bcf

	}
	
		
	Main();


?>

<?php

function Usage() {
	die("
		getAlignment.php 
		-w <Max gap size> 
		-c <Configuration File> 
		-s <Sparse cutoff 0,1> 
		-l <Shortest alignment to output> 
		-o <Output folder> 
		-u <Start from contig>
		-x yes|no <Penalize highly covered samples? must specify relative coverage (e.g. 4 ) compared to the lowest covered sample at the 3rd column in config file>
		-N yes|no Count N characters as differences for the variant filter. default: yes
		-R Character requirement filter file
		-T Use only these transcripts listed in this file
		-D yes|no toggle Divergence filter.
		-n yes|no whether N positions in the reference seq will be included in the alignment
		-b yes|no treat N as ref base. If you have a variant-only bcf file, then turning this on will turn all uncovered positions to ref base instead of N.
		The configuration file is a tab delimited file. Column 1 is taxon name(no space), Column two is fasta or bcf file. 
                The first line of the config file must contain the reference taxon and its fasta file. This script will treat N's as identical base to the reference.



");
}

function Main() {
	global $nMaxGapSize, $sRefGenome, $sRefTaxon, $sOutFolder , $arrBCFs,$arrCoverages, $bPenalizeCoverage , $arrBCFUtils, $sStartFromContig, $nShortestAln, $arrTranscriptList;

/*
	$oFastaHack = new FastaHack();
	$oFastaHack->SetDb("$sRefGenome");
*/
	$hRefFile = fopen($sRefGenome , "r");

	mkdir($sOutFolder);
	
	

	if (!$hRefFile ) {
		die("Cannot open Reference genome for reading: $sRefGenome ");
	}

	//$hProtein = fopen($sProteinOut  , "w");
	//$hDNAOut = fopen($sDNAOut , "w");
	

	$nRefIndex = 0;
	$sContigName = "";
	$sRefSeq = "";//ref seq for current contig
	$bStartProcess = ($sStartFromContig == "" )? true:false;

	while (($sLn = fgets($hRefFile )) !== false) { //parse through the genome ref file
		if (substr($sLn, 0, 1) == '>') { // title line
			if ($sRefSeq!="") {
				$sRefSeq = strtoupper($sRefSeq);
				if ($bStartProcess) {
					if( count($arrTranscriptList)>0 && !(array_key_exists($sContigName, $arrTranscriptList)) ) {
						echo("$sContigName not in the provided list, skipping...\n");
					}
					else{
						fnProcessContig($sContigName, $sRefSeq  );
					}
				}
				else if($sStartFromContig == $sContigName) {

					$bStartProcess = true;
					if( count($arrTranscriptList)>0 && !(array_key_exists($sContigName, $arrTranscriptList)) ) {
						echo("$sContigName not in the provided list, skipping...\n");
					}
					else{
						fnProcessContig($sContigName, $sRefSeq  );
					}
				}
				else {
					echo("Skipping contig $sContigName\n");
				}
			}
			$sContigName = trim(substr($sLn, 1));
			$sRefSeq = "";
			continue;
		}

		$sRefSeq .= trim($sLn);
	}

	

	if( count($arrTranscriptList)>0 && !(array_key_exists($sContigName, $arrTranscriptList)) ) {
		echo("$sContigName not in the provided list, skipping...\n");
	}
	else{
		fnProcessContig($sContigName, $sRefSeq  );
	}// don't forget the last contig, haha
	
	//fclose($hProtein);
	//fclose($hDNAOut);
	//echo($oFastaHack->GetCDSRegion("Scfld0", 1, 1000));

}

function fnProcessContig($sContigName, &$sRefSeq) {
	
	global $nMaxGapSize, $sRefGenome, $sRefTaxon, $sOutFolder, $arrBCFs,$arrCoverages, $bPenalizeCoverage, $arrBCFUtils, $nSparseCutoff, $nMinDepth, $nMinQual , $nShortestAln, $bIncludeRefNPos, $bCodeIndelOrMultiNonRefAlleleAsN, $bTreatNAsRefBase , $nAlternativeReadRatio;
	$nLenContig = strlen($sRefSeq);
	$nStart = 0;
	$nTaxa = count($arrBCFs) + 1;

	if ($nLenContig<$nShortestAln) {
		echo("$sContigName too short \n");
		return;
	}

	progressBar(1, $nLenContig, $sContigName);

	$nEmptyCount = 0; //when empty count exceeds the max gap limit, out put the alignment
	$nGoodColumnCount=0;
	$bBlockFound = false; //already found block?
	$nBlockStart = 0;
	$nBlockEnd = 0;
	$arrAlignment = array($sRefTaxon => ""); //multispecies alignment

	$nSparseCutoffTaxa = $nTaxa * $nSparseCutoff  ;
	$bStopContig = false;

	for ($i=0;$i< $nLenContig;$i++) {


			$arrCurrBase = array(); //key is taxon

		
			$sRefBase = substr($sRefSeq, $i, 1);


			if ( (!$bIncludeRefNPos) && $sRefBase == 'N') {
				$nEmptyCount ++;
				if ($nEmptyCount > $nMaxGapSize) {
					$nBlockEnd = $i;
					$nEmptyCount = 0;
					$bBlockFound = false;
					
					if ($nGoodColumnCount > 0) { 
						fnWriteOut($sContigName , $nBlockStart, $nBlockEnd, $arrAlignment, $nLenContig);
					}
					$nGoodColumnCount = 0;
				}
				continue;
			}

			$arrCurrBase[$sRefTaxon] = $sRefBase;
			$nNCount = 0;

			foreach($arrBCFUtils as $sTaxon => $oBCFTool) {
				//Now append the other taxa, given that they are only SNPs but not indels

				$arrPos1 = $oBCFTool->GetVariantsAt($sContigName , $i + 1); //remember that BCF is 1 based

				if ($arrPos1 === -1 ) { // error occurred
					//$bStopContig = true;
					//break;
					$arrPos1  = false; //treat error as N.
					
				}

				if ($arrPos1 === false) {
				
					$arrCurrBase[$sTaxon] = ($bTreatNAsRefBase)? $sRefBase:'N';//$sRefBase; no reads covering that region, say N
					$nNCount++;
					continue;
				}


				else if (strlen($arrPos1[0]) > 1 || strlen($arrPos1[1]) > 1 ) {
					if ($bCodeIndelOrMultiNonRefAlleleAsN) {
						$arrCurrBase[$sTaxon] = 'N';//$sRefBase; no reads covering that region, say N
						$nNCount++;
						continue;
					}
					else {
						$arrCurrBase = false; // do not use indels for phylogenetics
						break;
					}					
							
				}



				if (substr($arrPos1[0],0,1) !== $sRefBase ) {
					die("Reference base in fasta disagree with reference base given in bcf file, using the wrong ref? $sContigName : $i $sRefBase ".$arrPos1[0].PHP_EOL);
				}
				$sAltBase = ($arrPos1[1] == '.')? $sRefBase : $arrPos1[1]; // if the alternative base is identical with ref.
				
				$nRequiredMinDepth = $nMinDepth;
				if ($bPenalizeCoverage) {
					$nRequiredMinDepth *= $arrCoverages[$sTaxon];
				}

				$bExtraFilterCriteria = true;
				if ($sAltBase != $sRefBase) {
					$arrDepthForAlleles = explode(",", $arrPos1["DPR"]);
					if (count($arrDepthForAlleles) ==2) {
						$bExtraFilterCriteria = (($arrDepthForAlleles[1] / ($arrDepthForAlleles[0] + $arrDepthForAlleles[1] ))>=$nAlternativeReadRatio);
					} 
					else {
						$bExtraFilterCriteria = false;
					}
					
				}

				
				$arrCurrBase[$sTaxon] =  (($arrPos1[2] >= $nMinQual && $arrPos1["DP"]>= $nRequiredMinDepth && $bExtraFilterCriteria)? $sAltBase : 'N');

				//print_r($arrPos1);
				//print_r($bExtraFilterCriteria);
				//print_r($arrCurrBase[$sTaxon]);
				//die();

				if ( $arrCurrBase[$sTaxon] =='N') {
					$nNCount++;
				}
			}

			if ($bStopContig) {
				break;
			}




			//echo( $nNCount . "\t");
			//echo( ($nTaxa - 1)  . "\n") ;
			if ( $nNCount  > $nSparseCutoffTaxa) {				
				$nEmptyCount ++;
				//ob_end_clean();
				//echo ($nEmptyCount . " ");
				continue;
			}
			else if ( !$bBlockFound ) { // if the next block is not found, this is the start of the block
				$arrAlignment = array($sRefTaxon => ""); //multispecies alignment
				foreach($arrBCFs as $sTaxon => $sPath) {
					$arrAlignment[$sTaxon] = ""; //initialize elements for each taxon
				}
				$bBlockFound = true;
				$nEmptyCount = 0;
				$nBlockStart = $i;
				$nBlockEnd = 0;
				$nGoodColumnCount = 0;
			}
			else {
				$nGoodColumnCount++;
			}

			if ($arrCurrBase !== false) { // if this base is not deleted.
				foreach($arrCurrBase as $sTaxon => $sBase) {
					$arrAlignment[$sTaxon] .= $sBase;
				}
			}

			if ( ($nEmptyCount > $nMaxGapSize) || ($i == $nLenContig-1)) { // trigger write
				$nBlockEnd = $i;
				$nEmptyCount = 0;
				$bBlockFound = false;
				if ($nGoodColumnCount > 0) {
					fnWriteOut($sContigName , $nBlockStart, $nBlockEnd, $arrAlignment, $nLenContig);
				}
				$nGoodColumnCount = 0;
			}
			


	}

	progressBar($nLenContig, $nLenContig, $sContigName);

}

function fnWriteOut($sContigName , $nBlockStart, $nBlockEnd, $arrAlignment, $nLenContig) {
	global $nMaxGapSize, $sRefGenome, $sRefTaxon, $sOutFolder, $arrBCFs, $arrBCFUtils, $nShortestAln , $arrCharRequirements, $bNAsDifference , $bDivergenceFilter;

		$sContigName = preg_replace("/[^a-zA-Z0-9s.]/", "_", $sContigName);
		progressBar($nBlockEnd + 1, $nLenContig, $sContigName);
		//print_r($arrAlignment);
		//die();
		$arrAlignment2 = $bDivergenceFilter? fnExludeDivergentAlignment( $arrAlignment , $bNAsDifference  ) : $arrAlignment; // exclude ambiguous or divergent regions, treat "N"s as identity
		$nAlnlen = strlen($arrAlignment2[$sRefTaxon]);

		//print_r($arrAlignment2);
		//die();

		if ( $nAlnlen < $nShortestAln) {
			return;
		}

		//check effective character according to rule file:

		foreach($arrCharRequirements as $arrRule) {

			$nORTaxa = count($arrRule); // first item is the num of valid characters

			if ( $nORTaxa <= 1) continue;
			$nAtLeastChar = $arrRule[0];

			$bSatisfyFilter = false;
			for($nTaxonRule=1;$nTaxonRule<$nORTaxa;$nTaxonRule++) {
				if ( ($nAlnlen - substr_count($arrAlignment2[$arrRule[$nTaxonRule]] , 'N')) >= $nAtLeastChar) {
					$bSatisfyFilter  = $bSatisfyFilter | true;
				}
				else {
					$bSatisfyFilter  = $bSatisfyFilter | false;
				}
			}

			if (!$bSatisfyFilter) {
				return; // discard
			}
			
		}


		//$sContigName = 
		$sOutFile = "$sOutFolder/$sContigName"."_".($nBlockStart+1)."_".($nBlockEnd+1).".phy";
		$fOutFile = fopen($sOutFile , "w");

		fwrite($fOutFile , " ". count($arrAlignment2). " ". $nAlnlen. PHP_EOL);

		foreach($arrAlignment2 as $sTaxon => $sSeq) {
		
			fwrite($fOutFile , $sTaxon."           ");
			fwrite($fOutFile , $sSeq . PHP_EOL );
		

		}
		unset($arrAlignment);
		unset($arrAlignment2);
		

		fclose($fOutFile );


}

function progressBar($current=0, $total=100, $label="", $size=50) { 

    //Don't have to call $current=0 
    //Bar status is stored between calls 
    static $bars; 
    $new_bar = FALSE;
    if(!isset($bars[$label])) { 
        $new_bar = TRUE; 
        fputs(STDOUT,"$label Progress:\n"); 
    } 
    else if($current == $bars[$label]) return 0; 

    $perc = round(($current/$total)*100,2);        //Percentage round off for a more clean, consistent look 
    for($i=strlen($perc); $i<=4; $i++) $perc = ' '.$perc;    // percent indicator must be four characters, if shorter, add some spaces 

    $total_size = $size + $i + 3; 
    // if it's not first go, remove the previous bar 
    if(!$new_bar) { 
        for($place = $total_size; $place > 0; $place--) echo "\x08";    // echo a backspace (hex:08) to remove the previous character 
    } 
     
    $bars[$label]=$current; //saves bar status for next call 
    // output the progess bar as it should be 
    for($place = 0; $place <= $size; $place++) { 
        if($place <= ($current / $total * $size)) echo '[42m [0m';    // output green spaces if we're finished through this point 
        else echo '[47m [0m';                    // or grey spaces if not 
    } 

    // end a bar with a percent indicator 
    echo " $perc%"; 

    if($current == $total) { 
        echo "\n";        // if it's the end, add a new line 
        //unset($bars[$label]); 
    } 
}  

function fnB2S($b) {
	return $b? "yes":"no";
}

?>
