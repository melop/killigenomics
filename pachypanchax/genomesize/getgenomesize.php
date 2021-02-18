<?php
$nCores = 40;
$nCountingCore = 1;
$COUNTQUAL="/beegfs/group_dv/software/source/countquality/dist/Release/GNU-Linux/countquality";
$nMappingCore = $nCores - $nCountingCore;
$sLibDef = "libdef.txt";
//$sLibDef = "libdef_FFT.txt";
$sRefDef = "refs.txt";

$sMapDIR = "mapped";
$sCountDIR = "basecount";
$SAMTOOLS = "samtools1.2";

$hRefDef = fopen($sRefDef , "r");

$arrRefs = array();
while(false!==($sLn = fgets($hRefDef) )) {
	$sLn = trim($sLn);
	if ($sLn == "") continue;
	if ($sLn[0] == '#') continue;
	list($sRefname, $sRefseq, $sBuscoBed) = explode("\t", $sLn);
	fnAssertFile($sRefseq);
	fnAssertFile($sBuscoBed);
	$sFastaOrder = "$sRefseq.scforder.txt";
	$sSortedBuscoBed = "$sBuscoBed.sorted.bed";
	$sSortedBuscoBedFlag = "$sBuscoBed.sorted.bed.ok";

	if (!file_exists($sFastaOrder) ) {
		fnShell("$SAMTOOLS faidx $sRefseq; cut -f1,2 $sRefseq.fai > $sFastaOrder");
	}
	
	if (!file_exists($sSortedBuscoBedFlag) || !file_exists($sSortedBuscoBed)) {
		fnShell("hhvm sortbed.php $sFastaOrder $sBuscoBed $sSortedBuscoBed && touch $sBuscoBed.sorted.bed.ok");
	}
 	$arrRefs[$sRefname] = array('ref' => $sRefseq , 'buscobed' => $sSortedBuscoBed, 'fastaorder' => $sFastaOrder);
}

//print_r($arrRefs);

$hLibDef = fopen($sLibDef , "r");
$arrLibs = array();
while(false!==($sLn = fgets($hLibDef) )) {
	$sLn = trim($sLn);
	if ($sLn == "") continue;
	if ($sLn[0] == '#') continue;
	list($sSpeciesCode, $sSample , $sRef , $sRead1, $sRead2) = explode("\t", $sLn);
	fnAssertFile($sRead1);
	fnAssertFile($sRead2);
	if (array_key_exists($sSpeciesCode , $arrLibs)) {
		die("$sSpeciesCode is defined twice in the library description file!\n");
	}
	if (!array_key_exists($sRef , $arrRefs)) {
		die("$sRef is undefined in the reference definition file!\n");
	}
	if ($sRead1 == $sRead2) {
		die("$sRead1 and $sRead2 are the same\n");
	}

 	$arrLibs[$sSpeciesCode] = array('refname' => $sRef , 'r1' => $sRead1, 'r2' => $sRead2 , 'reffasta' => $arrRefs[$sRef]['ref'], 'refbuscobed' =>  $arrRefs[$sRef]['buscobed'] , 'reffastaorder' => $arrRefs[$sRef]['fastaorder']);
}

//print_r($arrLibs);

foreach($arrLibs as $sSpeciesCode => &$arrSample) {
	$sMapWD = "$sMapDIR/ref_". $arrSample['refname']."/";
	$sCountWD = "$sCountDIR/$sSpeciesCode";
	$sRefFasta = $arrSample['refname'];
	$sBUSCOBed = $arrSample['refbuscobed'];
	$sFastaOrder = $arrSample['reffastaorder'];
	fnShell("mkdir -p $sMapWD ; mkdir -p $sCountWD");
	$sCMD = "( if [ ! -e $sMapWD/$sSpeciesCode.ok ]; then bash makebams.sh $sMapWD \"".$arrSample['reffasta']."\" $nMappingCore \"".$arrSample['r1']."\" \"".$arrSample['r2']."\" $sSpeciesCode && touch $sMapWD/$sSpeciesCode.ok; fi; if [ ! -e $sMapWD/$sSpeciesCode.cov.ok ]; then bedtools coverage -a \"$sBUSCOBed\" -b $sMapWD/$sSpeciesCode.bam  -sorted -g $sFastaOrder  -d > $sMapWD/$sSpeciesCode.cov.txt && touch $sMapWD/$sSpeciesCode.cov.ok; fi ) &\n";
	$sCMD .= "( cd $sCountWD ; if [ ! -e done.ok ];then $COUNTQUAL 30  \"".$arrSample['r1']."\" \"".$arrSample['r2']."\" counts.txt 33 && touch done.ok ; fi ) & \n";
	$sCMD .= "wait;\n";
	fnShell($sCMD);
}

function fnShell($s) {
	echo($s."\n");
	exec($s);
}

function fnAssertFile($sF) {
	if (!file_exists($sF)) {
		die("Error: $sF does not exist!\n");
	}
}
?>
