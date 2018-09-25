<?php
include("lib.php");
//Joins the output of getTranscriptAlignmentALlowN.php (phylip) format to a single fasta format, where each scaffold is one record

$sPhylipPrefix = "testaln/GapFilledScaffold_";
$sScfldNamePattern = "/([^\\\/]_[0-9]+)_\S+.phy/";
$sScaffoldPrefix = "GapFilledScaffold_";
$nFrom = 1;
$nTo = 46729;
$sIndLabel = "YA_muscle1";
$sOutFile = "YA_muscle1.pseudogenome2.fasta";
$sRefGenome = "/beegfs/group_dv/home/RCui/killifish_genomes/map_bwamem/samtools1.2_relaxed/notfurScflds.fa";



while(count($argv) > 0) {
    $arg = array_shift($argv);
    switch($arg) {
        case '-p':
            $sPhylipPrefix = trim(array_shift($argv));
            break;
        case '-S':
            $sScaffoldPrefix  = trim(array_shift($argv));
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

    }
}



$hOut = fopen($sOutFile , "w");
$oRefFasta = new FastaHack();
$oRefFasta->SetDb($sRefGenome);

//echo($oRefFasta->GetContig($sScaffoldPrefix."1"));
//die();

$arrAllFiles = glob($sPhylipPrefix."*.phy");
$arrParsedFiles = array();
//print_r($arrAllFiles);

foreach($arrAllFiles as $sFile) {
	preg_match('/([^\/_]+_[0-9]+)_\S+.phy/', $sFile, $arrMatches);
	if (count($arrMatches)!=2) {
		die("Parse error: $sFile \n");
	}

	if (array_key_exists($arrMatches[1], $arrParsedFiles) ) {
		die("Two files found for the same contig, please set window size to a larger number! \n".$arrMatches[1]."\n" );
	}

	$arrParsedFiles[$arrMatches[1]] = $sFile;
}


for ($nContig=$nFrom; $nContig<=$nTo; $nContig++) {
	
	/*$arrFiles = glob($sPhylipPrefix.$nContig."_*.phy");
	if (count($arrFiles)>1) {
		die("Error: more than 1 phylip file found for contig $nContig. Make sure when you run getTranscriptAlignmentALlowN.php you specify a window size larger than any contig.\n");
	}
	*/


	fwrite($hOut , ">".$sScaffoldPrefix.$nContig.PHP_EOL);


	if (!array_key_exists($sScaffoldPrefix.$nContig , $arrParsedFiles) ) {
		echo("Warning: contig $nContig doesn't exist, will put N instead.\n");
		$sNs = str_repeat("N", strlen($oRefFasta->GetContig($sScaffoldPrefix.$nContig))); //get length of the contig
		fwrite($hOut , wordwrap($sNs, 75, PHP_EOL, true).PHP_EOL);
	} else {

		$sPhyFile = $arrParsedFiles[$sScaffoldPrefix.$nContig] ;//$arrFiles[0];
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
?>
