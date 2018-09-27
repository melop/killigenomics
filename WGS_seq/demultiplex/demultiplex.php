<?php
//demultiplex the following structured barcodes:
//i7 in a separate file
//i5 is embedded within read 1

$sSampleSheet = "samplesheet.txt";
$sRunDef = "rundefs.txt";

$nTrimBaseBeforeBarcodeI5 = 0; //trim off this many bases before taking the i5 barcode.
$nTrimBaseAfterBarcodeI5 = 1; //trim off the fixed "T" from read 1 after barcode.
$nAllowedMismatch = 1; //allow 1 base mismatch.
$nI5BarcodeLen = 0;
$nI7BarcodeLen = 0;
$sOutFolder = "demultiplexed/";

//read in sample sheet
$hSampleSheet = fopen($sSampleSheet , "r");
$arrBarcodes = array(); //first dimension is i5 barcode, second is i7
$arrI5List = array();
$arrI7List = array();

exec("mkdir -p $sOutFolder");

while(false!== ($sLn = fgets($hSampleSheet) )) {
	$sLn = trim($sLn);
	if ($sLn == "") continue;
	if ($sLn[0] == '#') continue;

	$arrF = explode("\t" , $sLn);

	if ($arrF[6] == "#N/A") continue; // 

	$sSampleID = trim($arrF[0]);
	$sI5Name = trim($arrF[6]);
	$sI7Name = trim($arrF[7]);
	$sI5seq = strtoupper(trim($arrF[8]));
	$sI7seq = strtoupper(trim($arrF[9]));
	$sSpecies = trim($arrF[12]);
	$sPop = trim($arrF[13]);

	if ($nI5BarcodeLen == 0) $nI5BarcodeLen=strlen($sI5seq);
	if ($nI7BarcodeLen == 0) $nI7BarcodeLen=strlen($sI7seq);

	if ($nI5BarcodeLen != strlen($sI5seq)) {
		die("Error, variable i5 barcode length not supported!\n");
	}

	if ($nI7BarcodeLen != strlen($sI7seq)) {
		die("Error, variable i7 barcode length not supported!\n");
	}

	if (!array_key_exists($sI5seq , $arrBarcodes) ) {
		$arrBarcodes[$sI5seq] = array();
		$arrI5List[] = $sI5seq;
	}

	if (array_key_exists($sI7seq , $arrBarcodes[$sI5seq]) ) {
		die("Error: barcode combination $sI5seq + $sI7seq already occupied by a previous sample!\n");
	}

	if (!in_array($sI7seq , $arrI7List )) {
		$arrI7List[] = $sI7seq;
	}

	$arrBarcodes[$sI5seq][$sI7seq] = array('i5name' => $sI5Name , 'i7name' => $sI7Name, 'species' => $sSpecies, 'pop' => $sPop , 'id' => $sSampleID, 'hRead1' => false, 'hRead2' => false, 'npairs' => 0, 'nbases1' => 0, 'nbases2' => 0, 'nbasetotal' =>0);
}

//print_r($arrBarcodes);
//print_r($arrI5List);
//print_r($arrI7List);

$arrMinDiff1 = fnCheckBarcodeDistance($arrI5List);

//print_r($arrMinDiff1);
 
$arrMinDiff2 = fnCheckBarcodeDistance($arrI7List); 

//print_r($arrMinDiff2);

if ($nAllowedMismatch >= $arrMinDiff1['mindiff'] || $nAllowedMismatch >= $arrMinDiff2['mindiff'] ) {
	die("The current allowed mismatches is too large for the set of barcodes. Set this to at most ".(min($arrMinDiff1['mindiff'] , $arrMinDiff2['mindiff'] )-1).PHP_EOL );
}

//perform checks to make sure that the number of mismatched bases in the barcodes are higher than the set number of allowed mismatches

//NOW read in the run definitions. 
$hRunDef = fopen( $sRunDef, "r");
$arrRuns = array();
while(false!== ($sLn = fgets($hRunDef) )) {
	$sLn = trim($sLn);
	if ($sLn == "") continue;
	if ($sLn[0] == '#') continue;

	$arrRuns[] = explode("\t" , $sLn);
}

//print_r($arrRuns);
//echo(fnCorrectBarcode( 'TAAGGCGT' , $arrI7List));

//now go over all runs

//make the unassigned output
$hAmb1 = fopen("$sOutFolder/ambiguous_R1.fq", "w");
$hAmb2 = fopen("$sOutFolder/ambiguous_R2.fq", "w");
$hAmb3 = fopen("$sOutFolder/ambiguous_R3.fq", "w");

foreach($arrRuns as $arrRun) {
	echo("Demultiplexing RUN: ".$arrRun[0].PHP_EOL);
	$sRead1 = trim($arrRun[1]);
	$sRead2 = trim($arrRun[2]);
	$sIdxRead = trim($arrRun[3]);

	// now use zcat -f to open these files as streams. this also transparently deals with gz format

	echo("Opening fastq file for reading. $sRead1 $sRead2 $sIdxRead\n");
	$pR1 = popen("zcat -f $sRead1", "r");
	$pR2 = popen("zcat -f $sRead2", "r");
	$pR3 = popen("zcat -f $sIdxRead", "r");

	/*if ( (!is_resource($pR1)) || (!is_resource($pR2)) || (!is_resource($pR3))) {
		die("Fail to open fastq file for reading. $sRead1 $sRead2 $sIdxRead\n");
	}*/

	//now read
	echo("Reading...\n");
	$nLnInRecord = -1; // first line 0 , sequence line 1, third 2, quality line 3
	$arrCurrRecord = array();
	while(true) {

		$sLn1 = fgets($pR1);
		$sLn2 = fgets($pR2);
		$sLn3 = fgets($pR3);

		if ($sLn1===false || $sLn2===false || $sLn3===false) {
			break;
		}
		$sLn1 = trim($sLn1);
		$sLn2 = trim($sLn2);
		$sLn3 = trim($sLn3);

		//echo("Hello");

		if ($sLn1 == '' && $sLn2 == '' && $sLn3 == '') continue;

		$nLnInRecord++;

		if ($nLnInRecord == 0) {
			if ( $sLn1[0] != "@" || $sLn2[0] != "@" || $sLn3[0] != "@") {
				die("Error: expecting fastq header line, instead found:\n$sLn1\n$sLn2\n$sLn3\n");
			}	
			//check if read name matches
			$sReadName1 = preg_split("/[\s\\\\\/\.]/", $sLn1);
			$sReadName2 = preg_split("/[\s\\\\\/\.]/", $sLn2);
			$sReadName3 = preg_split("/[\s\\\\\/\.]/", $sLn3);
			if ($sReadName1[0] != $sReadName2[0] || $sReadName2[0] != $sReadName3[0]) {
				die("Error: fastq files not synchronized: \n$sLn1\n$sLn2\n$sLn3\n");
			}
			$arrCurrRecord = array();
			$arrCurrRecord['readname1'] = $sReadName1[0] ;
			$arrCurrRecord['readname2'] = $sReadName2[0] ;
			$arrCurrRecord['readname3'] = $sReadName3[0] ;
			continue;
		}

		if ($nLnInRecord == 1) {
			$arrCurrRecord['seq1'] = $sLn1;
			$arrCurrRecord['seq2'] = $sLn2;
			$arrCurrRecord['seq3'] = $sLn3;
			continue;
		}

		if ($nLnInRecord == 2) {
			if ( $sLn1[0] != "+" || $sLn2[0] != "+" || $sLn3[0] != "+") {
				die("Error: expecting fastq quality line header, instead found:\n$sLn1\n$sLn2\n$sLn3\n");
			}	

			continue;
		}

		if ($nLnInRecord == 3) {
			$arrCurrRecord['qual1'] = $sLn1;
			$arrCurrRecord['qual2'] = $sLn2;
			$arrCurrRecord['qual3'] = $sLn3;
			fnProcessRecord($arrCurrRecord);
			$nLnInRecord = -1;
			continue;
		}

	}

}

//finally output report:

$hReport = fopen("$sOutFolder/reports.txt" , "w");
fwrite($hReport, "SampleID\tSpecies\tPop\ti5\ti7\ti5seq\ti7seq\tReadPairs\tRead1Bases\tRead2Bases\tTotalBases\n");
foreach($arrBarcodes as $sI5 => $arrSamples) {
	foreach($arrSamples as $sI7 => $oSample) {
		fwrite($hReport, $oSample['id']."\t".$oSample['species']."\t".$oSample['pop']."\t".$oSample['i5name']."\t".$oSample['i7name']."\t$sI5\t$sI7\t".$oSample['npairs']."\t".$oSample['nbases1']."\t".$oSample['nbases2']."\t".($oSample['nbases1'] + $oSample['nbases2'])."\n");
	}
}

function fnProcessRecord(&$arrCurrRecord) {
	global $nTrimBaseBeforeBarcodeI5 , $nTrimBaseAfterBarcodeI5, $nAllowedMismatch, $nI5BarcodeLen , $nI7BarcodeLen , $arrBarcodes, $arrI5List , $arrI7List , $sOutFolder;

	$sI5Barcode = strtoupper( substr($arrCurrRecord['seq1'], $nTrimBaseBeforeBarcodeI5, $nI5BarcodeLen));
	$sI7Barcode = strtoupper( $arrCurrRecord['seq3'] );
	if (strlen($sI7Barcode) != $nI7BarcodeLen) {
		die("Error: The length of the I7 barcode in read 3 differs from the sample sheet definition!\n");
	}

	if (!array_key_exists($sI5Barcode , $arrBarcodes) ) { // try correcting
		$sI5Barcode = fnCorrectBarcode($sI5Barcode , $arrI5List);
		if ($sI5Barcode === false) { //correction failed
			fnWriteNotFound($arrCurrRecord);
			return;
		}
	}

	if (!array_key_exists($sI7Barcode , $arrBarcodes[$sI5Barcode]) ) { // try correcting
		$sI7Barcode = fnCorrectBarcode($sI7Barcode , $arrI7List);
		if ($sI7Barcode === false) { //correction failed
			fnWriteNotFound($arrCurrRecord);
			return;
		}
	}

	if (!array_key_exists($sI7Barcode , $arrBarcodes[$sI5Barcode]) ) {
		fnWriteNotFound($arrCurrRecord);
		return;
	}

	$oSample = $arrBarcodes[$sI5Barcode][$sI7Barcode];

	//'i5name' => $sI5Name , 'i7name' => $sI7Name, 'species' => $sSpecies, 'pop' => $sPop , 'id' => $sSampleID, 'hRead1' => false, 'hRead2' => false

	if (false === $arrBarcodes[$sI5Barcode][$sI7Barcode]['hRead1'] ) {

		$sR1File = "$sOutFolder/".$oSample['species'].'_'.$oSample['pop'].'_'.$oSample['id'].'_'.$oSample['i5name'].'_'.$oSample['i7name'].'_R1.fq';
		$sR2File = "$sOutFolder/".$oSample['species'].'_'.$oSample['pop'].'_'.$oSample['id'].'_'.$oSample['i5name'].'_'.$oSample['i7name'].'_R2.fq';
		$arrBarcodes[$sI5Barcode][$sI7Barcode]['hRead1'] = fopen($sR1File , "w");
		$arrBarcodes[$sI5Barcode][$sI7Barcode]['hRead2'] = fopen($sR2File , "w");
	}

	$sTrimmedSeq1 = substr($arrCurrRecord['seq1'] , $nTrimBaseBeforeBarcodeI5 + $nI5BarcodeLen + $nTrimBaseAfterBarcodeI5);
	$sTrimmedQual1 = substr($arrCurrRecord['qual1'] , $nTrimBaseBeforeBarcodeI5 + $nI5BarcodeLen + $nTrimBaseAfterBarcodeI5);

	fwrite($arrBarcodes[$sI5Barcode][$sI7Barcode]['hRead1'], $arrCurrRecord['readname1']."\n".$sTrimmedSeq1."\n+\n".$sTrimmedQual1."\n");
	fwrite($arrBarcodes[$sI5Barcode][$sI7Barcode]['hRead2'], $arrCurrRecord['readname2']."\n".$arrCurrRecord['seq2']."\n+\n".$arrCurrRecord['qual2']."\n");

	//WRITE Down stats //'npairs' => 0, 'nbases1' => 0, 'nbases2' => 0, 'nbasetotal' =>0
	$arrBarcodes[$sI5Barcode][$sI7Barcode]['npairs'] += 1;
	$arrBarcodes[$sI5Barcode][$sI7Barcode]['nbases1'] += strlen($sTrimmedSeq1) ;
	$arrBarcodes[$sI5Barcode][$sI7Barcode]['nbases2'] += strlen($arrCurrRecord['seq2']) ;
	//$arrBarcodes[$sI5Barcode][$sI7Barcode]['nbasetotal'] += $arrBarcodes[$sI5Barcode][$sI7Barcode]['nbases1'] + $arrBarcodes[$sI5Barcode][$sI7Barcode]['nbases2'];
}

function fnWriteNotFound($arrCurrRecord) {
	global $hAmb1, $hAmb2, $hAmb3 ;

	fwrite($hAmb1 , $arrCurrRecord['readname1']."\n".$arrCurrRecord['seq1']."\n+\n".$arrCurrRecord['qual1']."\n");
	fwrite($hAmb2 , $arrCurrRecord['readname2']."\n".$arrCurrRecord['seq2']."\n+\n".$arrCurrRecord['qual2']."\n");
	fwrite($hAmb3 , $arrCurrRecord['readname3']."\n".$arrCurrRecord['seq3']."\n+\n".$arrCurrRecord['qual3']."\n");
	
}

function fnCheckBarcodeDistance($arrBarcodes) {
	$nMinDiff = 100;
	$sPair = "";
	for($i=0;$i<count($arrBarcodes)-1;$i++) {
		$sB1 = $arrBarcodes[$i];
		for($j=$i+1;$j<count($arrBarcodes);$j++) {
			$sB2 = $arrBarcodes[$j];
			$nDiff = strlen($sB1) - substr_count( $sB1 ^ $sB2 , "\0") ;// use XOR on the 2 strings. identical characters becomes \0
			if ($nDiff < $nMinDiff) {
				$nMinDiff = $nDiff;
				$sPair = "$sB1  -  $sB2";
			}
		}
	}

	return array('mindiff' => $nMinDiff , 'pair' => $sPair);
}

function fnCorrectBarcode($sBarcode, &$arrBarcodes) {
	global $nAllowedMismatch;
	$sBarcode = strtoupper($sBarcode);
	foreach($arrBarcodes as $sExpBarcode) {	
		$nDiff = strlen($sExpBarcode) - substr_count( $sExpBarcode ^ $sBarcode , "\0") ;	
		if ($nDiff <= $nAllowedMismatch) {
			return $sExpBarcode;
		}
	}

	return false;
}


?>
