<?php
//Sometimes the assembled MT DNA is presented in the fasta in such a way that it breaks a gene.
//This script adjust the fasta using MITOS annotations, such that the sequence starts from the 12s RNA and ends after cytb.

$sMitosOutDir = "out";
$sOutDir = "mt_pos_adjusted";

exec("mkdir -p $sOutDir");

$sFasta = $argv[1];

//$sFasta = "mt/PLP.fasta";

if (!file_exists($sFasta) ) {
	die("Fasta file cannot be found!\n");
}


$sFileName = basename($sFasta);

$sOutSpDir = "$sMitosOutDir/$sFileName";

$sBedFile = "$sOutSpDir/ret.bed";
echo($sBedFile . "\n");
//parse fasta
$sSeq = file_get_contents($sFasta );
$arrLines = explode("\n", $sSeq);
$sSeq = trim(implode("", array_slice($arrLines, 1) ));
echo("Len: " . strlen($sSeq) );

//parse bed
$arrGenes = array();
$hBed = fopen($sBedFile , 'r');
$nCutPoint = 0;
$bLookStrand = "+";
$nPrevEnd = 0;
$arrCandidateBreaks = array( );
$arrCandidateBreaks[] = 0;

while(false !== ($sLn = fgets($hBed))) {
	$sLn = trim($sLn);
	if ($sLn == '') continue;
	$arrF = explode("\t" , $sLn);
	$arrGenes[$arrF[1] ] = array('end' => $arrF[2]-1, 'gene' => $arrF[3]);

	$nCandidateBreak = (($arrF[1] - $nPrevEnd)-1)/2;

	if ($nCandidateBreak >= 6) {
		$nCandidateBreak = $arrF[1] - round($nCandidateBreak);
		$arrCandidateBreaks[] = $nCandidateBreak;
	} 

	$nPrevEnd = $arrF[2]-1;


	if (strpos( $arrF[3] , "rrnS") !== false ) {
		if (strpos( $arrF[3] , "_") !== false) { //cob is split
			$sSplitPart = (explode('_' , $arrF[3]))[1];
			if ($sSplitPart == 'a') {
				$bLookStrand = $arrF[5];
				$nCutPoint = ($bLookStrand=='+')? $arrF[1]:($arrF[2]-1);
			}
		} else {

				$bLookStrand = $arrF[5];
				$nCutPoint = ($bLookStrand=='+')? $arrF[1]:($arrF[2]-1);
		}
	}
} 

$arrCandidateBreaks[] =  strlen($sSeq)-1;

ksort($arrGenes);
print_r($arrGenes);
echo("$nCutPoint $bLookStrand \n");
print_r($arrCandidateBreaks);
$nSelectedBreakpoint = 0;

$nAddDirection = ($bLookStrand == '+')? -1:1;
$nStartLook = ($bLookStrand == '+')? count($arrCandidateBreaks)-1 :0;
for($i=$nStartLook; $i>=0 && $i<=(count($arrCandidateBreaks)-1); $i += $nAddDirection) {
	$nDiff = $nCutPoint - $arrCandidateBreaks[$i] ;
	if ($bLookStrand == '+' && $nDiff >= 0) {
		$nSelectedBreakpoint = $arrCandidateBreaks[$i];
		break;
	}

	if ($bLookStrand == '-' && $nDiff <= 0) {
		$nSelectedBreakpoint = $arrCandidateBreaks[$i];
		break;
	}
}

echo("break and rearrange at: $nSelectedBreakpoint\n");
if ($nSelectedBreakpoint == 0 || ($nSelectedBreakpoint == strlen($sSeq)-1) ) {
	echo("No change is needed.\n");
	exec("mkdir -p out_adjusted/");
	exec("ln -sf `realpath $sOutSpDir` out_adjusted/");
	die();
}

$sNewSeq = substr($sSeq , $nSelectedBreakpoint);
$sNewSeq .= substr($sSeq , 0, $nSelectedBreakpoint);

if (strlen($sNewSeq) != strlen($sSeq) ) {
	die("Something is wrong\n");
}

$sAdjustedFasta = "$sOutDir/$sFileName";
$hAdjF = fopen($sAdjustedFasta , 'w');
fwrite($hAdjF , ">mt\n$sNewSeq\n");

?>
