<?php
require_once(dirname(__FILE__) . "/lib.php");
$sMakerGFF = "maker.finalannot.gff";
$sWorkDir = "./improve_frag/";
$sGenome = "query.fa";
$sOutGFF = "maker.finalannot.improved.gff";

while(count($argv) > 0) {
    $arg = array_shift($argv);
    switch($arg) {
        case '-R':
            $sGenome = trim(array_shift($argv));
            break;
	case '-M':
	    $sMakerGFF  = trim(array_shift($argv));
	    break;
	case '-o':
	    $sWorkDir  = trim(array_shift($argv));
	    break;

    }
}

if (!file_exists($sWorkDir)) {
	die("Cannot find $sWorkDir\n");
}

if (!file_exists($sMakerGFF)) {
	die("Cannot find $sWorkDir\n");
}

$arrScfs = fnGetScfs($sMakerGFF) ;
$arrRNAToDelete = array();
$arrGeneToDelete = array();
$oMakerGFF = new MappedCDSGFF();
$oMakerGFF->LoadGFF($sMakerGFF , $sGenome);

$hOutGFF = fopen($sOutGFF, 'w');

$arrNewAnnotations = array();
$nFixedGenes = 0;
foreach( $arrScfs as $sScf => $bDummy ) {
	$sWD = "$sWorkDir/$sScf";
	$sFinishedFlag = "$sWD/scf.done";
	$sAnnotFile = "$sWD/annotations.notes";
	$sFixedGFF = "$sWD/fixed.gff";
	if (!file_exists($sFinishedFlag)  ) {
		echo("$sScf skipped...\n");
		continue;
	}

	if (!file_exists($sAnnotFile)  ) {
		echo("$sScf skipped...\n");
		continue;
	}
	
	$hAnnotFile = fopen($sAnnotFile , 'r');
	while(false !== ($sLn = fgets($hAnnotFile) ) )  {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		$arrFxx = explode("\t" , $sLn);
		$sRNAID = $arrFxx[0];
		$sKeep = $arrFxx[1];
		$sNewAnnot = (count($arrFxx) == 3)? $arrFxx[2]:'.';


		if ($sKeep == 'delete') {
			$arrRNAToDelete[$sRNAID] = true;
			$oRNA = $oMakerGFF->GetmRNA('maker', $sRNAID);
			$sDNAID = $oRNA['annot']['Parent'];
			$arrGeneToDelete[$sDNAID] = true;
		} 
	}

	if (file_exists($sFixedGFF)) {
		$sFixedGFFContent = file_get_contents($sFixedGFF);
		$nFixedGenes += substr_count($sFixedGFFContent , "\tgene\t");
		fwrite($hOutGFF , $sFixedGFFContent); 
	}

}



$hMakerGFF = fopen($sMakerGFF , "r");


while( false !== ($sLn = fgets($hMakerGFF) )) {

	$sLn = trim($sLn);
	if ($sLn == "") continue;
	if ($sLn[0] == '#') {
		fwrite($hOutGFF , $sLn."\n");
	}

	$arrF = explode("\t" , $sLn);

	if ($arrF[2] == "gene") {
		preg_match("/ID=([^;]+)/", $arrF[8] , $arrM);
		if (count($arrM) == 2) {
			$sID = trim($arrM[1]);
			if (array_key_exists($sID , $arrGeneToDelete)) {
				echo("remove gene $sID\n");
				continue;
			}
		}
	} else if ($arrF[2] == "mRNA") {
		preg_match("/Parent=([^;]+)/", $arrF[8] , $arrM);
		if (count($arrM) == 2) {
			$sParentID = trim($arrM[1]);
			if (array_key_exists($sParentID , $arrGeneToDelete)) {
				
				echo("remove mRNA of $sParentID\n");
				continue;
			}
		}

	} else {
		preg_match("/Parent=([^;]+)/", $arrF[8] , $arrM);
		if (count($arrM) == 2) {
			$sID = trim($arrM[1]);
			if (array_key_exists($sID , $arrRNAToDelete)) {
				echo("remove part $sID\n");
				continue;
			}
		}
	}

	if ($arrF[2] == "mRNA") {
		preg_match("/ID=([^;]+)/", $arrF[8] , $arrM);
		if (count($arrM) == 2) {
			$sID = trim($arrM[1]);
			if (array_key_exists($sID , $arrNewAnnotations)) {
				
				echo("Update annotation for mRNA $sID ...\n");
				$arrF[8] = trim($arrNewAnnotations[$sID]);
				fwrite($hOutGFF , implode("\t" , $arrF) ."\n");
				continue;
			}
		}

	} 

	fwrite($hOutGFF , $sLn."\n");

}

echo(count($arrGeneToDelete) . " fragmented gene models replaced with $nFixedGenes fixed gene models\n" );
echo(count($arrNewAnnotations) . " quality annotations updated\n" );

function fnGetScfs($sGFF) {
	$hIn = popen("cut -f1,1 $sGFF | sort | uniq" ,'r');
	$arr = array();
	while(false !== ($sLn = fgets($hIn))) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		$arr[$sLn] = true;
	}

	return $arr;
}


?>
