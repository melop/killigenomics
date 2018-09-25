<?php
$sRAxMLCmd = "/software/source/standard-RAxML-8.2.9/raxmlHPC-PTHREADS-SSE3 -T 40 -f a -p 767 -x 23333 -#100 -m PROTGAMMAJTTX";
$nTotalParts = 1;
$nThisPart = 0;
$sFileList = "cleaned_fa_list.txt";
$sDIR = "ClusteRs";

while(count($argv) > 0) {
    $arg = array_shift($argv);
    switch($arg) {
        case '-N':
            $nTotalParts  = trim(array_shift($argv));
            break;
        case '-f':
            $nThisPart  = trim(array_shift($argv));
            break;
        case '-L':
            $sFileList  = trim(array_shift($argv));
            break;
        case '-D':
            $sDIR  = trim(array_shift($argv));
            break;

    }
}

$nJobCounter = -1;
exec("mkdir -p $sDIR");

$hFileList = fopen($sFileList , "r");

while( false!== ($sLn = fgets($hFileList) ) ) {
	$sLn = trim($sLn);

	if ($sLn == "") continue;

        $nJobCounter++;

        if ($nJobCounter % $nTotalParts != $nThisPart) {
                continue;
        }

	$sOutStem = str_replace( ".fa" , ".out", $sLn);
	$sRAxMLTree = $sDIR . "/RAxML_bipartitions.".$sOutStem;

	if (file_exists($sRAxMLTree) ) continue;
	echo($sRAxMLTree." calculating ...\n");
	//otherwise, delete all raxmloutput of this group

	exec("rm $sDIR/RAxML_*.$sOutStem");
	exec("cd $sDIR; $sRAxMLCmd -s ../UPhO_Seqs/$sLn -n $sOutStem");

}

?>
