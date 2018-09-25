<?php
$sObsNull = "reports/1rate.log.txt";
$sObsAlt = "reports/2rate.log.txt";

$arrFolder = array( "compute_sims/logs");
$sOut = "sim_ln_diff.txt";
$arrNull = array();
foreach($arrFolder as $sFolder) {
	$arrNull = array_merge( $arrNull  , glob("$sFolder/rep*.null.log"));
}

//print_r($arrNull);

$hOut = fopen($sOut, 'w');

list( $nObsNullLk, $sObsNullEst) = fnGetLk($sObsNull);
list( $nObsAltLk , $sObsAltEst)  = fnGetLk($sObsAlt);

fwrite($hOut , "num\tnull_model\taltermodel\tLogLikelihoodRatio\tnull_est\talt_est\n");
fwrite($hOut , "Obs\t$nObsNullLk\t$nObsAltLk\t". ($nObsAltLk - $nObsNullLk)."\t$sObsNullEst\t$sObsAltEst\n");

$nSimCount = 0;
$nBiggerThanObs = 0;
foreach($arrNull as $sNull) {
	//echo("Reading $sNull\n");
	$sAlt = str_replace(".null.log", ".alt.log", $sNull);
	if (!file_exists($sAlt)) {
		echo("File $sAlt not found\n");
		continue;
	}


	list($nNullLk, $sNullEst )= fnGetLk($sNull);
	list($nAltLk , $sAltEst) = fnGetLk($sAlt);

	if ($nNullLk === false || $nAltLk === false) {
		if ($nNullLk === false ) {
			echo("cannot parse: $sNull \n");
		}

		if ($nAltLk === false) {
			echo("cannot parse: $sAlt \n");
		}
		continue;
	}
	$nSimCount++;

	if ( ($nAltLk - $nNullLk) > ($nObsAltLk - $nObsNullLk ) ) {
		$nBiggerThanObs++;
	}

	fwrite($hOut , "$nSimCount\t$nNullLk\t$nAltLk\t". ($nAltLk - $nNullLk)."\t$sNullEst\t$sAltEst\n");
}

echo("\n$nBiggerThanObs / $nSimCount p = " . ($nBiggerThanObs / $nSimCount) );
echo("\nwritten to $sOut\n");
function fnGetLk($sLog) {
	$hLog = fopen($sLog,'r');

	$arrScoreIndex = array();
	$nLk = 0;
	$sEst = "";
	while(false!==($sLn = fgets($hLog))) {

		preg_match("/Lambda : (\S+) Mu : (\S+) & Score: (\S+)/",  $sLn, $arrM);
		if (count($arrM) == 4) { //result line
			$sKey = "lambda: ". trim($arrM[1]) . "_mu: " .trim($arrM[2]) ;
			$arrScoreIndex[$sKey] = $arrM[3];
			continue;
		}

		preg_match("/Lambda : (\S+) & Score: \d+Mu : (\S+) & Score: \d+/", $sLn, $arrM);
		if (count($arrM) == 3) { //result line
			$sKey = "lambda: ". trim($arrM[1]) . "_mu: " .trim($arrM[2]) ;
			if (!array_key_exists($sKey , $arrScoreIndex)) {
				echo("Error in $sLog: $sKey not found\n");
				return array(false, false);
			}

			$nLk = $arrScoreIndex[$sKey];
			if ($nLk == '-inf') {
				echo("Error in $sLog: model did not converge.\n");
				return array(false, false);
			}

			$sEst = $sKey;
		}
	}

	if ($nLk == 0) {
		echo("Log file not finished: $sLog\n");
		return array(false, false);
	}

	return array( $nLk , $sEst);
}

?>
