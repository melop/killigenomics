<?php
$sFileList = "datalist.txt";
$sOutBash = "_download.sh";

$sOutDir = "data";

exec("mkdir -p $sOutDir");

$hOut = fopen($sOutBash , 'w');
$hIn = fopen($sFileList, 'r');

while( false !== ($sLn = fgets($hIn)) ) {
	$sLn = trim($sLn);

	if ($sLn == '') continue;
	if ($sLn[0] == '#') continue;

	$arrF = explode("\t", $sLn );
	list($sFishId , $nTissueId) = explode('_', trim($arrF[0]));
	$sStrain = trim($arrF[2]);
	$nAge = trim($arrF[3]);
	$sTissue = trim($arrF[4]);
	$sSRRID = trim($arrF[5]);

	$sOutFastqName = "$sOutDir/$sStrain.$sTissue.$nAge.wk.$sFishId.fq.gz";
	$sFlag = "$sOutDir/$sOutFastqName.done";
	$sCmd = "if [ -e $sFlag ]; then\n echo $sFlag finished; \nelse\n fastq-dump -Z --skip-technical  --readids --read-filter pass --dumpbase --clip --defline-qual '+' --defline-seq '@\$sn' $sSRRID | gzip -c > $sOutFastqName && touch $sFlag &\n fi \n";

	fwrite($hOut , $sCmd);
}
	fwrite($hOut , "\nwait\n");
fclose($hOut);
//exec($sOutBash);

?>
