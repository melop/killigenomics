<?php
$sURL="http://mitos.bioinf.uni-leipzig.de/result.py";
$sDownloadURL="http://mitos.bioinf.uni-leipzig.de/";
$sOutDir = "out";

if (array_key_exists(2, $argv)) {
        $sOutDir = $argv[2];
}

$sFasta = $argv[1];

if (!file_exists($sFasta) ) {
	die("Fasta file cannot be found!\n");
}


$sFileName = basename($sFasta);

$sOutSpDir = "$sOutDir/$sFileName";

$sJobIdFile = "$sOutSpDir/jobid.txt";
$hJobIdF = fopen($sJobIdFile , 'r');
$sJobId = fgets($hJobIdF);
fclose($hJobIdF);

if ($sJobId === false) {
	die("Error: cannot get job id for job $sFasta");
}

$sJobId = trim($sJobId);

$sDoneFlag = "$sOutSpDir/ret.gff";
while (true) {
	if (file_exists($sDoneFlag)) {
		echo("Finished. $sFasta $sJobId finished.");
		break;
	}

	$nRet = fnTryDownload($sOutSpDir , $sJobId);

	if ($nRet === 0) {
		sleep(30);
		continue;
	}

	if ($nRet === -1) {
		die("Cannot parse output for $sFasta $sJobId  \n");
	}
}

function fnTryDownload($sOutSpDir , $sJobId) {
	global $sURL, $sDownloadURL;
	$sRequest = "$sURL?hash=$sJobId";
	$sResponse = file_get_contents($sRequest );

	if ($sResponse === false) {
		echo("Server not responding, wait 1 min...\n");
		sleep(60);
		return 0;
	}

	//echo($sResponse);
	if (strpos($sResponse , "notified when your results are ready") !== false) {
		return 0;
	}

	preg_match_all('/(download\.py\?\S+)\"/', $sResponse , $arrM);

	if (count($arrM)!=2) {
		return -1;
	}

	foreach($arrM[1] as $sGetURL) {
		preg_match("/type=([^&]+)/", $sGetURL, $arrM2);
		$sType = $arrM2[1];
		echo("\"$sDownloadURL/$sGetURL\"\n");
		exec("cd $sOutSpDir; wget -O ret.$sType \"$sDownloadURL/$sGetURL\"");
	}

	preg_match_all('/(picture\.py\?\S+big)\"/', $sResponse , $arrM);

	foreach($arrM[1] as $sGetURL) {
		exec("cd $sOutSpDir; wget -O figure.png \"$sDownloadURL/$sGetURL\"");
	}

	exec("touch $sOutSpDir/mitos.done");

}

?>
