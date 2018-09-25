<?php
$sRefDef = "../hyphyrelax/export_protein_multiref_23_july_2017/refs.txt";
$sBamRoot = "../out/bams/";
$sSH = "batchrun.sh";

$hSH = fopen($sSH, 'w');

/*============ Load the ref genomes ======================================*/
$hRefDef = fopen($sRefDef , "r");
$arrRefGTFs = array();

$nLn = -1;
echo("Parsing reference genome definitions...\n");
while( ($sLn=fgets($hRefDef))!==false ) {
        $sLn = trim($sLn);
        if ($sLn == "") {
                continue;
        }
        $nLn++;
        if ($nLn==0) {
                continue; // skip header
        }
        
        $arrFields = explode("\t", $sLn);
        $arrRefGTFs[$arrFields[0]] = array( $arrFields[1], $arrFields[2]);
}


foreach($arrRefGTFs as $sRefSp => $arrInfo) {
	$arrSppFolders = glob($sBamRoot."/$sRefSp/*");

	foreach($arrSppFolders as $sFolder) {
		$sFolder = realpath($sFolder);
		$sSp = basename($sFolder);
		$sCmd = "hhvm get_exon_cov.php -G " . $arrRefGTFs[$sRefSp][1] . " -g " . $arrRefGTFs[$sRefSp][0] . " -o ". $sRefSp . " -S " . $sSp . " -b " . $sFolder . "/sorted.bam &\n";
		fwrite( $hSH, $sCmd);
	}
}

fwrite($hSH, "wait\n");
exec("bash $sSH");
?>
