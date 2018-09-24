<?php
$sOldContig = "contigs.fa"; //"contigs.reduced.fa";//"contigs.fa";//"contigs.reduced.fa";

$sBlackList = "contaminated.list.txt";

$sOut = $sOldContig.".decontaminated.fa";
$sOutCont = $sOldContig.".contaminants.fa";

$hBlackList = fopen($sBlackList ,"r");
$hOldContig = fopen($sOldContig ,"r");
$hOut = fopen($sOut ,"w");
$hOutCont = fopen($sOutCont ,"w");

$arrBlackList = array();

while(false !== ($sLn = fgets($hBlackList)) ) {
	$arrF = explode("\t" , trim($sLn) );
	$arrBlackList[$arrF[0]] = $arrF[3];
}

//print_r($arrBlackList);
echo(count($arrBlackList)." blacklist items loaded\n");
$bOut = true;

while( false !== ($sLn = fgets($hOldContig) )) {
	if ($sLn[0] == ">") {
		$sContigName = substr(trim($sLn), 1);
		if (array_key_exists($sContigName , $arrBlackList) ) {
			echo("$sContigName excluded...\n");
			$bOut = false;
			fwrite($hOutCont , $sLn);
		} else {
			fwrite($hOut , $sLn);
			$bOut = true;
		}

		continue;
	}

	if ($bOut) {
			fwrite($hOut , $sLn);
	} else {
			fwrite($hOutCont , $sLn);
	}
}


?>
