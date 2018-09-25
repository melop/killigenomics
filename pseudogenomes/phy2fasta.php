<?php
include("lib.php");
//Joins the output of getTranscriptAlignmentALlowN.php (phylip) format to a single fasta format, where each scaffold is one record

$sPhylipPrefix = "testaln/";
$nFrom = 0;
$nTo = -1;
$sIndLabel = "YA_muscle1";
$sOutFile = "YA_muscle1.pseudogenome2.fasta";
$sRefGenome = "";



while(count($argv) > 0) {
    $arg = array_shift($argv);
    switch($arg) {
        case '-p':
            $sPhylipPrefix = trim(array_shift($argv));
            break;
	case '-f':
	    $nFrom  = trim(array_shift($argv));
	    break;
	case '-t':
	    $nTo  = trim(array_shift($argv));
	    break;
	case '-L':
	    $sIndLabel = trim(array_shift($argv));
	    break;
	case '-o':
	    $sOutFile = trim(array_shift($argv));
	    break;
	case '-R':
	    $sRefGenome = trim(array_shift($argv));
	    break;

    }
}



$hOut = fopen($sOutFile , "w");
$arrScfldNames = fnGetFastaNames($sRefGenome );
$oRefFasta = new FastaHack();
$oRefFasta->SetDb($sRefGenome);

echo("Found ".count($arrScfldNames)." scaffold names in reference\n");
if ($nTo == -1) {
	$nTo = count($arrScfldNames) - 1;
}



$arrAllFiles = glob($sPhylipPrefix."*.phy");
$arrParsedFiles = array();
//print_r($arrAllFiles);

foreach($arrAllFiles as $sFile) {
	preg_match('/(\S+)_[0-9]+_[0-9]+.phy/', basename($sFile), $arrMatches);
	if (count($arrMatches)!=2) {
		die("Parse error: $sFile \n");
	}

	if (array_key_exists($arrMatches[1], $arrParsedFiles) ) {
		die("Two files found for the same contig, please set window size to a larger number! \n".$arrMatches[1]."\n" );
	}

	//echo($arrMatches[1].PHP_EOL);
	$arrParsedFiles[$arrMatches[1]] = $sFile;
}


for ($nContig=$nFrom; $nContig<=$nTo; $nContig++) {
	
	/*$arrFiles = glob($sPhylipPrefix.$nContig."_*.phy");
	if (count($arrFiles)>1) {
		die("Error: more than 1 phylip file found for contig $nContig. Make sure when you run getTranscriptAlignmentALlowN.php you specify a window size larger than any contig.\n");
	}
	*/

	$sScfldName = $arrScfldNames[$nContig];
	fwrite($hOut , ">".$sScfldName .PHP_EOL);


	if (!array_key_exists($sScfldName  , $arrParsedFiles) ) {
		echo("Warning: contig $nContig ($sScfldName) doesn't exist, will put N instead.\n");
		$sNs = str_repeat("N", strlen($oRefFasta->GetContig($sScfldName))); //get length of the contig
		fwrite($hOut , wordwrap($sNs, 75, PHP_EOL, true).PHP_EOL);
	} else {

		$sPhyFile = $arrParsedFiles[$sScfldName] ;//$arrFiles[0];
		echo($sPhyFile." ...".PHP_EOL);
		$hPhyFile = fopen($sPhyFile , "r");
		while(false!==($sLn = fgets($hPhyFile))) {
			$sLn = trim($sLn);
			if ($sLn == "") continue;
			if (strpos($sLn, $sIndLabel." ")===0) {
				//echo("a\n");
				fwrite($hOut , wordwrap(trim(explode(" ", $sLn, 2)[1]), 75, PHP_EOL, true).PHP_EOL);
			}
		}

	}

	
}

function fnGetFastaNames($sFasta ) {
	$hFas = fopen($sFasta, "r");
	$arrNames = array();
	while(false !== ($sLn = fgets($hFas))) {
		$sLn = trim($sLn);
		if (substr($sLn , 0,1) == ">" ) {
			$arrNames[] = substr($sLn , 1) ;
		}
	}

	return $arrNames;
}
?>
