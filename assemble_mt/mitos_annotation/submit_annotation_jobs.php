<?php
//preg_match('|result.py\?hash\=(\S+)\\\\"|', 'result.py?hash=2234ayesdsWWH\"', $arrM);
//echo($arrM[1]);
//die();
$sURL="http://mitos.bioinf.uni-leipzig.de/upload.py";
$sOutDir = "out";
if (array_key_exists(2, $argv)) {
	$sOutDir = $argv[2];
}

exec("mkdir -p $sOutDir");

$sFasta = $argv[1];

if (!file_exists($sFasta) ) {
	die("Fasta file cannot be found!\n");
}

$sFileName = basename($sFasta);

$sOutSpDir = "$sOutDir/$sFileName";
exec("mkdir -p $sOutSpDir");

$sRewriteFasta = "$sOutSpDir/in.fas";
$hRewriteFasta = fopen($sRewriteFasta , 'w');
$sSeq = file_get_contents($sFasta );
$arrLines = explode("\n", $sSeq);
$sHeader = $arrLines[0];
$sSeq = trim(implode("", array_slice($arrLines, 1) ));
$sSeq = preg_replace("/[^ATCGatcg]/", 'N', $sSeq);

fwrite($hRewriteFasta , "$sHeader\n$sSeq\n");

if (file_exists($sOutSpDir."/jobid.txt")) { //check if previous submission is there
	$hJobIdTest = fopen($sOutSpDir."/jobid.txt" , 'r');
	$sLn1 = fgets($hJobIdTest);
	fclose($hJobIdTest);
	if ($sLn1 !== false && strpos($sLn1 , 'DOCTYPE HTML PUBLIC')===false ) {
		die("$sOutSpDir already submitted\n");
	}
}



$sOut = "$sOutSpDir/jobid.txt";
$hOut = fopen($sOut , 'w');

if (!function_exists('curl_file_create')) {
    echo("Fallback on old curl function\n");
    function curl_file_create($filename, $mimetype = '', $postname = '') {
        return "@$filename;filename="
            . ($postname ?: basename($filename))
            . ($mimetype ? ";type=$mimetype" : '');
    }
}


# http://php.net/manual/en/curlfile.construct.php


// Create a CURLFile object / oop method 
$oPostfile = curl_file_create($sRewriteFasta,'text/plain',$sFileName); // uncomment and use if the upper procedural method is not working.

// Assign POST data
$oPostData = array('name' => "Rongfeng Cui", 
'email' => "rcui@age.mpg.de",
'jobid' => $sFileName, 
 'code' => '2',
'proceed' => '',
'advanced' => '',
'prot' => true,
'trna' => true,
'rrna' => true,
'multi' => true,
'evalue' => '2',
'cutoff' => '50',
'maxovl' => '20',
'clipfac' => '10',
'fragovl' => '20',
'fragfac' => '10',
'ststrange' => '6',
'finovl' => '35',
'myFile' => $oPostfile);

$curl = curl_init();
curl_setopt($curl, CURLOPT_URL, $sURL);
curl_setopt($curl, CURLOPT_USERAGENT,'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Ubuntu Chromium/53.0.2785.143 Chrome/53.0.2785.143 Safari/537.36');
curl_setopt($curl, CURLOPT_HTTPHEADER,array('User-Agent: Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Ubuntu Chromium/53.0.2785.143 Chrome/53.0.2785.143 Safari/537.36','Referer: http://mitos.bioinf.uni-leipzig.de/index.py','Content-Type: multipart/form-data'));
curl_setopt($curl, CURLOPT_SSL_VERIFYPEER, false); // stop verifying certificate
curl_setopt($curl, CURLOPT_RETURNTRANSFER, true); 
curl_setopt($curl, CURLOPT_POST, true); // enable posting
curl_setopt($curl, CURLOPT_POSTFIELDS, $oPostData); // post images 
curl_setopt($curl, CURLOPT_FOLLOWLOCATION, true); // if any redirection after upload
$r = curl_exec($curl); 
curl_close($curl);

preg_match('|result.py\?hash\=(\S+)\\\\"|', $r, $arrM);
if (count($arrM) != 2 ) {
	echo("Cannot parse response from server, please check output directly\n");
	fwrite($hOut, $r);
} else {
	echo("Job hash code is: $arrM[1]\n");
	fwrite($hOut, $arrM[1]."\n");
	fwrite($hOut, $r."\n");
}


?>
