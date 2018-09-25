<?php
$sFolder = "genfamily_sim";
$sErrModelFolder = realpath("errmodels/");
$arrFiles = glob("$sFolder/rnd_*.tab");
$sOutDIR = "compute_sims";
$sNullTemplate = "cafe1ratelambdamu.sh.template";
$sAlternativeTemplate = "cafe2ratelambdamu.sh.template";


exec("mkdir -p $sOutDIR/logs/");

$nTotalParts = 1;
$nThisPart = 0;

while(count($argv) > 0) {
    $arg = array_shift($argv);
    switch($arg) {
        case '-N':
            $nTotalParts  = trim(array_shift($argv));
            break;
        case '-f':
            $nThisPart  = trim(array_shift($argv));
            break;

    }
}



for($nFile=0;$nFile<count($arrFiles);$nFile++) {


        if ($nFile % $nTotalParts != $nThisPart) {
                continue;
        }


	$sInFile = $arrFiles[$nFile];
	$sNullScript = "$sOutDIR/rep$nFile.null.sh";
	$sAltScript = "$sOutDIR/rep$nFile.alt.sh";
	$sInputTable = realpath($sInFile);
	$sNullLogFile = realpath("$sOutDIR/logs/")."/rep$nFile.null.log";
	$sAltLogFile = realpath("$sOutDIR/logs/")."/rep$nFile.alt.log";

	$sNullFlag = realpath("$sOutDIR/logs/")."/rep$nFile.null.done";
	$sAltFlag = realpath("$sOutDIR/logs/")."/rep$nFile.alt.done";





	$sErrModelFolder = realpath($sErrModelFolder);
	$arrNullVars = array('%COUNT_TABLE%' => $sInputTable , '%LOG_FILE%' => $sNullLogFile, '%ERR_MODEL_FOLDER%' => $sErrModelFolder);
	$arrAltVars = array('%COUNT_TABLE%' => $sInputTable , '%LOG_FILE%' => $sAltLogFile, '%ERR_MODEL_FOLDER%' => $sErrModelFolder);

	fnTemplateWrite($sNullTemplate , $sNullScript , $arrNullVars );
	fnTemplateWrite($sAlternativeTemplate , $sAltScript , $arrAltVars );

	if (!file_exists($sNullFlag)) {
		exec("if [ -e $sNullLogFile ]; then rm $sNullLogFile; fi;");
		exec("$sNullScript && touch $sNullFlag");
	} else {
		echo("$sNullLogFile is done. skip\n");
	}

	if (!file_exists($sAltFlag)) {
		exec("if [ -e $sAltLogFile ]; then rm $sAltLogFile; fi;");
		exec("$sAltScript && touch $sAltFlag");
	} else {
		echo("$sAltLogFile is done. skip\n");
	}




	//die();
}

function fnTemplateWrite($sTemplateFile, $sOutputFile, $arrVars) {
	$sRet = file_get_contents($sTemplateFile);
	$sRet =	str_replace( array_keys($arrVars), array_values($arrVars), $sRet );
	$hOut = fopen( $sOutputFile , "w");
	fwrite($hOut , $sRet);
	fclose($hOut);
	exec("chmod 755 $sOutputFile");
	
}

?>
