<?php
require_once(dirname(__FILE__) . "/lib.php");
require_once(dirname(__FILE__) . "/scores.lib.php");
$sMakerGFF = "maker.genesymbol.supple.gff"; //this gff must have ensemblorthologs annotation. the completeness will be assessed by alignment the gene model to the ensembl ortholog
$sGenome = "query.fa";
$sProteinFasta = "/beegfs/group_dv/home/RCui/killifish_genomes/annotation/UPhO_final/allproteins.fa";
$sWorkDir = "./finalannot/";
$nMinBlastPIdentity = 50; //filter HSP if lower than this.
$nMaxBlastPIdentityDiff = 20; //If identity of an HSP block is 20% smaller than other blocks in the same alignment, then exclude this hsp
$nTreatAsBigGap = 2; // if larger than this, treat as "big gap".
$nTreatAsHugeGap = 10; // if larger than this, treat as "huge gap".
$sDoOnly = ""; //only do this md5 mRNA name
$sForceRedo = false;

$nThisPart = 0;
$nTotalParts = 1;


while(count($argv) > 0) {
    $arg = array_shift($argv);
    switch($arg) {
        case '-R':
            $sGenome = trim(array_shift($argv));
            break;
	case '-M':
	    $sMakerGFF  = trim(array_shift($argv));
	    break;
	case '-P':
	    $sProteinGFF  = trim(array_shift($argv));
	    break;
	case '-C':
	    $sCandidateGFF  = trim(array_shift($argv));
	    break;
	case '-D':
	    $sDoOnly  = trim(array_shift($argv));
	    break;
	case '-p':
	    $sProteinFasta  = trim(array_shift($argv));
	    break;
	case '-o':
	    $sWorkDir  = trim(array_shift($argv));
	    break;
	case '-N':
	    $nTotalParts  = intval(trim(array_shift($argv)));
	    break;
	case '-f':
	    $nThisPart  = intval(trim(array_shift($argv)));
	    break;
	case '-m':
	    $nMinBlastPIdentity  = floatval(trim(array_shift($argv)));
	    break;
	case '-d':
	    $nMaxBlastPIdentityDiff  = floatval(trim(array_shift($argv)));
	    break;
	case '-g':
	    $nTreatAsBigGap  = intval(trim(array_shift($argv)));
	    break;
	case '-G':
	    $nTreatAsHugeGap  = intval(trim(array_shift($argv)));
	    break;
	case '-F':
	    $sForceRedo = true;

    }
}

exec("mkdir -p $sWorkDir");

$oProteinFasta = new FastaHack();
$oProteinFasta->SetDb($sProteinFasta);

$oMakerGFF = new MappedCDSGFF();
$oMakerGFF->LoadGFF($sMakerGFF , $sGenome);




$arrGenes = $oMakerGFF->GetGenes('maker'); // get all maker genes

$nGeneCount = -1;

foreach($arrGenes as $sGeneID => &$oGene) {

	$nGeneCount++;
	if ( $nGeneCount % $nTotalParts != $nThisPart) continue;

	echo("Looking at $sGeneID ...\n");

	$sGeneDIR = "$sWorkDir/".$oGene['scf']."/".md5($sGeneID);

	foreach($oGene['mRNAs'] as $smRNAID => &$omRNA) {

		$smRNAIDMD5= md5($smRNAID);

		if ($sDoOnly != '') {
			if ($sDoOnly != $smRNAIDMD5) {
				continue;
			}
		}

		//check if annotation already exists
		$oAnnot = $omRNA['annot'];
		if (array_key_exists( 'score',$oAnnot)) {
			echo("model already scored. Skip\n");
			continue;
		} else if ( (!array_key_exists( 'ensemblorthologs',$oAnnot)) || $oAnnot['ensemblorthologs'] == '') {
			echo("model not scored. it doesn't have ensembl ortholog, cannot score.\n");
			continue;
		}

		echo("scoring...\n");

		$sRNADIR = $sGeneDIR."/".$smRNAIDMD5;
		$sScoreNotes = "$sRNADIR/annotation.notes";

		if (!$sForceRedo) {
			if ( file_exists($sScoreNotes) ) { //if there are already results
				if (filesize($sScoreNotes) > 0) {
					echo("computed, skip.\n");
					continue;
				}
			}
		}

		exec("mkdir -p $sRNADIR");

		$hScoreNotes = fopen($sScoreNotes , 'w');

		echo(" -- mRNA  $smRNAID in dir $sRNADIR ...\n");

		$arrExtracted = $oMakerGFF->ExtractmRNASequence('maker', $smRNAID);

		// now score this main model:
		$sMainModelProt = $arrExtracted['AA'];
		$arrEnsemblIDs = explode(',' , $oAnnot['ensemblorthologs']);

		$nHighScore = -100000;
		$oHighScore = false;
		foreach($arrEnsemblIDs as $sEvidEnsemblID) {
			$sEnsemProt = fnGetEnsemblProt( $sEvidEnsemblID );
			$oScore = fnScoreExonerateModel($sMainModelProt , $sEnsemProt, $sRNADIR , true); // ignore internal stops 
			if ($oScore['score'] > $nHighScore ) {
				$oHighScore = $oScore;
				$nHighScore = $oScore['score'];
			}
		}

		if ($oHighScore !== false ) {
			fwrite( $hScoreNotes , "$smRNAID\tannot\t".fnArr2Annot($oHighScore)."\n");
		} else {
			fwrite( $hScoreNotes , "$smRNAID\tannotfailed\tNA\n");
		}
		
		
	}
}

function fnGetEnsemblProt( $sEvidEnsemblID ) {
	global $sProteinFasta, $oProteinFasta;
	exec("grep \"$sEvidEnsemblID\" $sProteinFasta ", $arrOut);
	if (count( $arrOut ) < 1 ) {
		return false;
	}

	$sID = trim(substr($arrOut[0] , 1));
	return $oProteinFasta->GetContig($sID);	

}



?>
