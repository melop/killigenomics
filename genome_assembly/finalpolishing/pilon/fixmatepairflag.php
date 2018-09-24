<?php
// bwa marks all mate pair reads as improperly paired. This is problematic for pilon
// This program sets the properly pair flag for mate pairs, and improper pair flag for innies
$hIn = fopen('php://stdin', 'r');
$nPECutoff = 500; //if separation less than this, treat as improper.

const SAM_PROPER_PAIR =  2;
const SAM_READ_UNMAPPED = 4;
const SAM_MATE_UNMAPPED = 8;
const SAM_FIRST_IN_PAIR= 64;
const SAM_SECOND_IN_PAIR= 128;
const SAM_READ_REV = 16;
const SAM_MATE_REV = 32;
const SAM_PCR_DUP = 1024;

$sPrevReadScf = "";
$sPrevMateScf = "";
$nPrevReadPos = -1;
$nPrevMatePos = -1;


while(false !== ($sLn = fgets($hIn ) )) {
	$sLn = trim($sLn);
	if ($sLn == "") continue;
	if ($sLn[0] == '@') {
		echo($sLn."\n");
		continue;
	}

	$arrM = explode("\t", $sLn);

	$F = intval($arrM[1]);
	$nLen = $arrM[8];

	$sReadScf = $arrM[2];
	$sMateScf = ($arrM[6] == '=')? $sReadScf : $arrM[6];
	$nReadPos = $arrM[3];
	$nMatePos = $arrM[7];

	

	if ( ( ($F & SAM_READ_UNMAPPED ) > 0 ) || ( ($F & SAM_MATE_UNMAPPED) > 0) ) {
		//echo("unmapped flag: $F\n");
		echo($sLn."\n");
		continue;
	}

	//echo("$sReadScf $nReadPos $sMateScf $nMatePos \n");
	//check if is PCR dup
	$bPCRDup = ( ($sReadScf == $sPrevReadScf) && ($sMateScf  == $sPrevMateScf) && ($nReadPos == $nPrevReadPos) && ($nMatePos == $nPrevMatePos) ) ;

	if ($bPCRDup) {
		$F = $F | SAM_PCR_DUP; //mark as pcr duplicate.
	}

	$sPrevReadScf = $sReadScf;
	$sPrevMateScf = $sMateScf;
	$nPrevReadPos = $nReadPos;
	$nPrevMatePos = $nMatePos;


	//both are mapped:
	$bFirstInPair = (($F & SAM_FIRST_IN_PAIR ) > 0 );
	$bSecondInPair = (($F & SAM_SECOND_IN_PAIR ) > 0);
	$bReadRev = (($F & SAM_READ_REV ) > 0);
	$bMateRev = (($F & SAM_MATE_REV ) > 0);

	//echo("first in pair: $bFirstInPair \t sec: $bSecondInPair \t rerev: $bReadRev \t materev: $bMateRev\t");

	$bOuttie = ($bFirstInPair && $bReadRev && ( $nLen > 0 ) ) || ($bSecondInPair && $bMateRev && ($nLen < 0) ) || ($bFirstInPair && $bMateRev && ( $nLen < 0 ) ) || ($bSecondInPair && $bReadRev && ($nLen > 0) );

	$bNewFlag = 0;

	if ($bOuttie && (abs($nLen) > $nPECutoff) ) {
		$bNewFlag = $F | SAM_PROPER_PAIR; 
		//echo("Oldflag : $F Newflag: $bNewFlag \t");
	} else {
		//echo("not pass: \t");
		$bNewFlag = $F & ~SAM_PROPER_PAIR;  // delete the proper pair flag for PE reads and pairs with very short inserts
	}

	$arrM[1] = $bNewFlag;

	echo(implode("\t" , $arrM)."\n" );
}

?>
