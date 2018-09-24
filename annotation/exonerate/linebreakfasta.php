<?php
//change line length in fasta file.
$nLineWidth = 100;
$sIn = "plp.fa";
$sOut = "plp_reformat.fa";

$hInFasta = fopen($sIn, "r");
$hOutFasta = fopen($sOut, "w");
$sSeqName = "";
$sSeq = "";
$nIndex = 0;


do {

	$sLn = fgets($hInFasta);
	//echo($sLn);
	if ( $sLn===false || substr($sLn , 0, 1) == '>' ) {
		//process current record.
		//echo("Header line!\n");
		if ($sSeq != "") {
				fnProcess($nIndex , $sSeqName , $sSeq); 

		}
		if ($sLn===false) {
			break;
		}

		$nIndex++;
		$sSeqName = trim($sLn);
		$sSeq = "";

	} else {
		$sLn = trim($sLn);
		$sSeq .= $sLn ;
	}

} while (true);

function fnProcess($nIndex , $sSeqName , $sSeq) {
	global $hOutFasta, $nLineWidth;
	fwrite( $hOutFasta  , $sSeqName.PHP_EOL . wordwrap ( $sSeq, $nLineWidth, "\n" , true ).PHP_EOL);
}

?>
