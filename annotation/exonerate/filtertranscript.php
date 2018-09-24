<?php
//filter transcript: type=complete, length >= 500bp

$sIn = "Trinity.cdhitout.fasta"; #output from cdshit
$sInGFF = "Trinity.fasta.transdecoder.gff3"; #transcoder annotations

$sOut = "Trinity.cdhitout.filtered.fa";
$sOutGFF = "Trinity.fasta.transdecoder.filtered.gff3";
$sOutAnnot = "Trinity.fasta.transdecoder.filtered.annot";
$nMinLen = 500;

$hIn = fopen($sIn, "r");
$hInGFF = fopen($sInGFF , "r");

$hOut = fopen($sOut, "w");
$hOutGFF = fopen($sOutGFF, "w");
$hOutAnnot = fopen($sOutAnnot , "w"); //this is for exonerate

$bOutput = false;

$arrKept = array();
$nTotalLen = 0;
while(false !== ($sLn = fgets($hIn)) ) {
	if (substr($sLn, 0, 1) == '>' ) {
		$bOutput = false;
		if ( strpos($sLn , "type:complete") === false) {
			continue;
		}

		preg_match("/len:(\d+)/", $sLn, $arrMatches);

		if (count( $arrMatches ) !=2) {
			echo("Warning, length not found in header : $sLn");
		}

		if ( ($arrMatches[1] * 3) >= $nMinLen) {
			list($sTransName, $sTmp) = explode("|", $sLn, 2);
			$arrKept[substr($sTransName,1)] = true;
			fwrite($hOut , $sTransName.PHP_EOL); 
			$bOutput = true;
			$nTotalLen += $arrMatches[1] * 3;
		} 
	} else if ($bOutput) {
		fwrite($hOut , $sLn); 
	}
}

echo("Kept ". count($arrKept). " transcripts, average length:" .  $nTotalLen / count($arrKept) .PHP_EOL);
//print_r($arrKept);

$sCurrTransName = "";
$nmRNAStart = 0;
while(false !== ($sLn = fgets($hInGFF)) ) {

	$sLn = trim($sLn);
 	if ($sLn == "") continue;

	list($sTransName, $sTmp) = explode("\t", $sLn, 2);

	if (array_key_exists($sTransName , $arrKept) ) {
		fwrite($hOutGFF, $sLn.PHP_EOL); 

		$arrFields = explode("\t", $sLn);

		if ($arrFields[2] == "mRNA") {
			$nmRNAStart = ($arrFields[6]=="+")? $arrFields[3] : $arrFields[4] ;
			$sCurrTransName = $arrFields[0];
		}

		if ($arrFields[2] == "CDS" && $arrFields[0] == $sCurrTransName) {
			$nLen = abs($arrFields[4] - $arrFields[3])+1;
			$nStart = ($arrFields[6]=="-")? $nmRNAStart - $arrFields[4] + 1: $arrFields[3] - $nmRNAStart + 1 ;
			fwrite($hOutAnnot , $arrFields[0]."\t+\t".$nStart."\t".$nLen.PHP_EOL);
		}
	}
}

?>
