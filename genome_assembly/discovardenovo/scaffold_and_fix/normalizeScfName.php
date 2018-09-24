<?php
#this scripts changes the scaffold names to a simple, numbered format, and output a corresponding table
$sIn = "";
$sOut = "";
$sPrefix = "";

while(count($argv) > 0) {
    $arg = array_shift($argv);
    switch($arg) {
	case '-I':
	    $sIn  = trim(array_shift($argv));
	    break;

	case '-O':
	    $sOut  = trim(array_shift($argv));
	    break;

	case '-P':
	    $sPrefix  = trim(array_shift($argv));
	    break;

    }
}

$sLog = "$sIn.old.new.scfnames.txt";

$hIn = fopen($sIn  , "r");
$hOut = fopen($sOut  , "w");
$hLog = fopen($sLog , "w");

fwrite($hLog , "#Oldname\tNewname"); 
$nScfld=0;
while( false !== ($sLn = fgets($hIn) ) ) {
	$sLn = trim($sLn);
	if ($sLn == "") continue;

	if ($sLn[0] == ">") {
		$nScfld++;
		$sNewName = $sPrefix.$nScfld;
		$sName = substr($sLn, 1);
		fwrite($hOut , ">$sNewName\n");
		fwrite($hLog , "$sName\t$sNewName\n"); 
		continue;
	}
	fwrite($hOut , "$sLn\n");
}
?>
