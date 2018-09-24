<?php
//use blast to accelerate exonerate search.


$arrOpts = getopt("i:a:n:t:r:w:o:");

$sInFasta = $arrOpts["i"];
$sInAnnotation = $arrOpts["a"];
$nThisPart = $arrOpts["n"];
$nTotalPart = $arrOpts["t"];
$sRefGenome = $arrOpts["r"];
$sWorkingDIR = $arrOpts["w"]; 
$bOverWrite = ($arrOpts["o"] == "yes"); #over write the results?
$sExonerate = "/beegfs/group_dv/software/source/maker3.00/exe/exonerate/bin/exonerate"; # full path to the exonerate command

exec("mkdir -p \"$sWorkingDIR\"");

$hInFasta = fopen($sInFasta, "r");
$sSeqName = "";
$sSeq = "";
$nIndex = -1; //0based

echo("Doing part $nThisPart of $nTotalPart...\n");
do {

	$sLn = fgets($hInFasta);
	//echo($sLn);
	if ( $sLn===false || substr($sLn , 0, 1) == '>' ) {
		//process current record.
		//echo("Header line!\n");
		if ($sSeq != "") {
			//echo("currindex $nIndex\n");
			if ( ($nIndex % $nTotalPart + 1) == $nThisPart) {
				//echo("Processing $sSeqName\n");
				fnProcess($nIndex+1 , $sSeqName , $sSeq); 
			}

		}
		if ($sLn===false) {
			break;
		}

		$nIndex++;
		$sSeqName = trim($sLn);
		$sSeq = "";

	} else {
		$sLn = trim($sLn);
		$sSeq .= $sLn ;
	}

} while (true);

function fnProcess($nIndex , $sSeqName , $sSeq) {
	global $bOverWrite, $sRefGenome , $sWorkingDIR, $sInAnnotation, $sExonerate;
	echo("Processing $sSeqName\n");
	$sSeqDIR = $sWorkingDIR."/".$nIndex."/";
	$sDoneFile = $sSeqDIR."done.0";
	$sQueryFile = $sSeqDIR."query.fa";
	$sBlastOut = $sSeqDIR."blast_ret.txt";
	$sCoordOut = $sSeqDIR."coordinates.txt";
	$sStrandOut = $sSeqDIR."strands.txt";
	$sOutGFF = $sSeqDIR."out.gff";
	$sFinalGFF = $sSeqDIR."picked.gff";

	//trim polyA tail

	$sSeq = fnTrimPolyATail($sSeq);

	if (!is_dir($sSeqDIR)) {
		exec("mkdir \"$sSeqDIR\"");
	}

	if (file_exists($sDoneFile)) {
		if ($bOverWrite) {
			echo("Warning, overwrite previous results in $sSeqDIR\n");
			exec("rm -r $sSeqDIR/*");
		} else {
			echo("Previous results unchanged, doing nothing $sSeqDIR\n");
			return;
		}

	}

	//write query fasta
	$hQueryFile = fopen($sQueryFile , "w");
	fwrite($hQueryFile , "$sSeqName\n".$sSeq);
	fclose($hQueryFile);

	//now run tblastx 
	fnExec("tblastx  -db $sRefGenome -outfmt 7 -query $sQueryFile > $sBlastOut");
	fnExec("php getbesthit.php -i $sBlastOut -o $sCoordOut -c $sStrandOut -x 300000 -s 300000 -f $sQueryFile");

	//check if any hits found:
	if ( (!file_exists($sCoordOut) ) || (filesize($sCoordOut)==0) ) {
		echo("Warning: no blast hit found for seq $sSeqName. \n");
		fclose(fopen($sDoneFile, "w"));
		return;
	}


	fnExec("fastahack $sRefGenome -c < $sCoordOut > $sSeqDIR/fasthackseq.txt");
 	fnExec("paste <(sed -e 's/^/\>/' $sCoordOut) $sSeqDIR/fasthackseq.txt | tr \\\"\\t\\\" \\\"\n\\\" > $sSeqDIR/fasthackseq.fas");
 	fnExec("mkdir $sSeqDIR/splitted; split -l 2 $sSeqDIR/fasthackseq.fas $sSeqDIR/splitted/seq;");
 	//fnExec("for f in $sSeqDIR/splitted/seq*;  do mv \\\"\$f\\\" \\\"\$f.fa\\\"; done");
	//rename
	$arrSplitFas = glob("$sSeqDIR/splitted/seq*");
	fclose(fopen($sOutGFF , "w"));
	foreach($arrSplitFas as $sF) {
		rename($sF , $sF.".fa");
		$sTargetFa = $sF.".fa";
		//now perform exonerate :

		fnExec("$sExonerate --model cdna2genome --query $sQueryFile --target $sTargetFa --querytype dna --targettype dna --showvulgar yes --showalignment yes --showtargetgff yes --bestn 1   --softmasktarget yes --dpmemory 5000 --fsmmemory 5000 --annotation $sInAnnotation >> $sOutGFF ");
	}

	//now run exonerate on every ref.

	fnExec("php pickbestmodel2.php -f $sOutGFF -o $sFinalGFF -k 1");
	fclose(fopen($sDoneFile, "w"));
}

function fnExec($sCmd) {
	echo("executing: $sCmd\n");
	shell_exec("/bin/bash -c \"$sCmd\"");
}

function fnTrimPolyATail($s) {
	return preg_replace("/([aA]{2})([aA]{5,})$/", "$1", $s);
}

?>
