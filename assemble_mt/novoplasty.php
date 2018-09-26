<?php
$sMainList = "list.txt"; // list of input data
$sMainList = "list2.txt";

$sOutDIR = "out";
$nPerJobThread = 2;
$nPerJobMemGB = 50; 
$sNovoPlasty = "/software/source/NOVOPlasty/NOVOPlasty2.6.3.pl"; //novoplasty script
$SCRATCHDIR=getcwd()."/tmp/";

exec("mkdir -p $sOutDIR");
exec("mkdir -p $SCRATCHDIR");

$sNovoPlastyDIR = "$sOutDIR/novoplasty";

exec("mkdir -p $sNovoPlastyDIR");

list($arrRefs , $arrSpecies) = fnLoadList($sMainList); 

//print_r($arrRefs);
//print_r($arrSpecies);

$arrJobIDs = array();
//spawn mapping jobs:
foreach($arrSpecies as $sRef => $arrSppInRef) {
	$sRefGenome = $arrRefs[$sRef];

	foreach( $arrSppInRef as $sSp => $arrReads) {
		$sSpNovoPDIR = "$sNovoPlastyDIR/$sRef/$sSp";

		exec("mkdir -p $sSpNovoPDIR");
		$sSpNovoPDIR = realpath($sSpNovoPDIR);
		$sNovoPBatch = "$sSpNovoPDIR/nvp.sbatch";
		$sNovoPConfig = "$sSpNovoPDIR/config.txt";
		$bNovoPDoneFlag = "$sSpNovoPDIR/assembly.done";

		$sCmd = "";
		$sR1File = "";
		$sR2File = "";
		if (count($arrReads['r1']) == 1 ) {
			$sR1File = $arrReads['r1'][0];
			$sR2File = $arrReads['r2'][0];
		} else {
			$sR1File = "$sSpNovoPDIR/concat_R1.fq";
			$sR2File = "$sSpNovoPDIR/concat_R2.fq";
			$sCmd .= "( if [ -e $sR1File ]; then echo $sR1File exists; else zcat -f " . implode(" " , $arrReads['r1']) . " > $sR1File fi; ) & \n" ;
			$sCmd .= "( if [ -e $sR2File ]; then echo $sR2File exists; else zcat -f " . implode(" " , $arrReads['r2']) . " > $sR2File fi; ) & \nwait;\n" ;
		}

		//write config.
		exec("( zcat -f $sR1File | head -n 5000 | tail -n 1 | wc -c ) 2>/dev/null" , $arrRet);
		if (count($arrRet) == 0 || intval($arrRet[0]) == 0 ) {
			die("Fail to get read length for file $sR1File \n");
		}

		$nReadLen = intval($arrRet[0]) - 1;
		print_r($arrRet);
		$nReadLen = ($nReadLen <= 105)? 100 : (($nReadLen <= 155)? 150 : (($nReadLen <= 255)? 250 :  $nReadLen));

		$nInsert = ($nReadLen == 100)? 300 : ( ($nReadLen == 150)? 350 : 400);

		echo("$sSp $nReadLen $nInsert\n");
		$hConfig = fopen($sNovoPConfig , 'w');
		fwrite($hConfig , 
"Project name        = $sSp
Insert size          = $nInsert
Insert size aut      = yes
Read Length          = $nReadLen
Type                 = mito
Genome Range         = 12000-22000
K-mer                = 39
Insert Range         = 1.6
Insert Range strict  = 1.2
Single/Paired        = PE
Max memory           = $nPerJobMemGB
Extended log         = 0
Save assembled reads = no
Combined reads       = 
Forward reads        = $sR1File
Reverse reads        = $sR2File
Seed Input           = $sRefGenome
Reference            = 
Chloroplast sequence = 

" );

		fclose($hConfig );

		$sCmd .= "if [ -e $bNovoPDoneFlag ];then echo $bNovoPDoneFlag was finished. skip.; else $sNovoPlasty -c $sNovoPConfig && touch $bNovoPDoneFlag; fi;";
		if (!file_exists($bNovoPDoneFlag)) {
			fnSubmitSlurm($sCmd , $nPerJobThread, $nPerJobMemGB.'G', $sNovoPBatch, $sSpNovoPDIR, $arrJobIDs); //append job id to job ids
		} else {
			echo("$bNovoPDoneFlag is done, skip.\n");
		}

	}
}

echo("Mapping jobs submitted, waiting for return...\n");
fnWaitForJobs($arrJobIDs);
echo("All Mapping jobs returned...\n");


function fnLoadList($sMainList) {

	$hIn = fopen($sMainList , 'r');
	fgets($hIn); //discard header line
	$arrRetRef = array();
	$arrRetSp = array();
	while( false !== ( $sLn = fgets($hIn) )) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		$arrF = explode("\t" , $sLn);
		if (count($arrF) != 7) continue;
		if ($arrF[1] == $arrF[2]) { // ref def.
			$arrRetRef[$arrF[1]] = $arrF[5];
			if (!file_exists($arrF[5].".sa") ) {
				die("Reference genome not yet indexed by bwa, please do that! : " . $arrF[5] ."\n");
			} 
		} else {
			if (!array_key_exists($arrF[1] , $arrRetSp )) {
				$arrRetSp[$arrF[1]] = array();
			}

			if (array_key_exists($arrF[2], $arrRetSp[$arrF[1]]) ) {
				die("Error: species " . $arrF[2] . " already defined!\n");
			}

			$arrRetSp[$arrF[1]][$arrF[2]] = array('r1' => explode(',',$arrF[5]) , 'r2' => explode(',' , $arrF[6]) );
		}
	}

	return array($arrRetRef , $arrRetSp);
}

function fnSubmitSlurm($sCmd , $nPerJobThread, $sMem, $sBatch, $sWD = false, &$arrJobIDs) {
	if ($sWD === false) {
		$sWD = getcwd();
	}
	$sHeader = "#!/bin/bash\n#SBATCH -n 1\n#SBATCH -c $nPerJobThread\n#SBATCH --mem=$sMem\n#SBATCH -p blade\ncd " . $sWD ."\n";
	$hBatch = fopen($sBatch , 'w');
	fwrite($hBatch , $sHeader . $sCmd ."\n");
	//return;
	exec("cd $sWD; sbatch $sBatch" , $arrOut);
	if (count($arrOut) < 1) {
		die("Failed to submit slurm job.\n$sCmd\n");
	}

	preg_match("/\d+/" , $arrOut[0], $arrM);
	if (count($arrM) != 1) {
		die("Cannot understand sbatch output\n");
	}

	$sJobID = $arrM[0];
	echo("submitted: $sJobID\n");
	$arrJobIDs[$sJobID] = true;
}

function fnWaitForJobs(&$arrJobIDs) {
	while(count($arrJobIDs) > 0 ) {
		foreach($arrJobIDs as $sJobID => $bDum) {
			exec("squeue -j $sJobID 2>/dev/null | wc -l " , $arrOut);
			if (trim($arrOut[0]) == '0') { // if job is done
				unset($arrJobIDs[$sJobID]);
			}
		}

		sleep(30);
	}
}
?>
