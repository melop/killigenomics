<?php
$sDir = "Relax_Pachypanchax";
$sOut = "Pachypanchax.omega.txt";
$hIn = popen("cat $sDir/ret_part* | grep -v \"BUSTED error\"", 'r');
$hOut = fopen($sOut,'w');

fwrite($hOut, "GeneName\tParameterization\tSuccess\tGeneSymbol\tAALen\ttmpfolder\tSuccess\tLRT\tP\tConstraint_LK\tUnconstraint_LK\tMG94_LK\tMG94_double_LK\tMG94_triple_LK\tGeneFullName\tTest_omega3\tTest_omega3_prop\tBg_omega3\tBg_omega3_prop\n");

while(false !== ($sLn = fgets($hIn) ) ) {
	$sLn = trim($sLn);
	if ($sLn == '') continue;
	$arrF = explode("\t", $sLn);
	if (count($arrF) < 6  ) continue;
	if ($arrF[2] != 'Success') continue;
	$sBUSTEDDIR = $arrF[5];
	$sRet = "$sBUSTEDDIR/in.nex.BUSTED.json";
	if (!file_exists($sRet)) {
		echo("$sRet no long exists\n");
		continue;
	}	

	$arrJSONret = json_decode(file_get_contents($sRet) , true);
	if (is_null($arrJSONret)) {
		    echo 'Last error: ', json_last_error_msg(), PHP_EOL, PHP_EOL;

		echo("$sRet cannot be loaded as json\n");
		continue;
	}

	$nFgOmega3 = $arrJSONret["fits"]["Unconstrained model"]["Rate Distributions"]["Test"]['2']['omega'];
	$nFgOmega3Prop = $arrJSONret["fits"]["Unconstrained model"]["Rate Distributions"]["Test"]['2']['proportion'];

	$nBgOmega3 = $arrJSONret["fits"]["Unconstrained model"]["Rate Distributions"]["Background"]['2']['omega'];
	$nBgOmega3Prop = $arrJSONret["fits"]["Unconstrained model"]["Rate Distributions"]["Background"]['2']['proportion'];
	fwrite($hOut, implode("\t", $arrF)."\t$nFgOmega3\t$nFgOmega3Prop\t$nBgOmega3\t$nBgOmega3Prop\n");
	

	
}



?>
