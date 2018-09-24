<?php
/*
This is run after orthology identification.
It scans through the gff file, and look at gene models that are incomplete (<90%).
If an ensembl ortholog has been identified, then try running the global genewise model
to see if a more complete model can be obtained.
*/
require_once(dirname(__FILE__) . "/lib.php");
require_once(dirname(__FILE__) . "/scores.lib.php");
$sGBlastPath = "/beegfs/group_dv/software/source/genblast/genblast_v139";
putenv("GBLAST_PATH=$sGBlastPath");



$sProtFile = "/beegfs/group_dv/home/RCui/killifish_genomes/annotation/UPhO_final/allproteins.fa";
$sMakerModels = "maker.finalannot.gff";
$sGenome = "query.fa";
$sMaskedGenome = "query.masked.fa";//"../../../../repeatmasker/PLP_v1.0_pilon/scf.fa.masked";
$sOutDir = "./improve_frag";



$HMMER = "/beegfs/group_dv/software/source/hmmer-2.3.2/src/hmmbuild";

$nThisPart = 0;
$nTotalParts = 1;
$bForceRedo = false;

while(count($argv) > 0) {
    $arg = array_shift($argv);
    switch($arg) {
	case '-N':
	    $nTotalParts  = intval(trim(array_shift($argv)));
	    break;
	case '-f':
	    $nThisPart  = intval(trim(array_shift($argv)));
	    break;
	case '-F':
	    $bForceRedo  = true;
	    break;
	}
}

$arrScfs = fnGetScfs($sMakerModels) ;
$oProteinFasta = new FastaHack();
$oProteinFasta->SetDb($sProtFile);

$oMakerGFF = new MappedCDSGFF();
$oMakerGFF->LoadGFF($sMakerModels , $sGenome);

$arrEnsemblProteinMap = fnGetEnsemblProteinNames($sProtFile);

$nScfCount = -1;
//$arrScfs = array("NORv12scf102" => true);

foreach( $arrScfs as $sScf => $bDummy ) {


	$nScfCount++;
	if ( $nScfCount % $nTotalParts != $nThisPart) continue;

	$sWD = "$sOutDir/$sScf";
	$sFinishedFlag = "$sWD/scf.done";
	if (file_exists($sFinishedFlag) && (!$bForceRedo ) ) {
		echo("$sScf already done. skip\n");
		continue;
	}

	echo("Doing $sScf ...\n");


	$sScfGFF = "$sWD/scf.gff";
	$sScfmRNA = "$sWD/mRNA.gff";
	$sScfSortedmRNA = "$sWD/mRNA.sorted.gff";
	$sScfMaskedFa = "$sWD/scf_masked.fa";
	$sScfFa = "$sWD/scf.fa";

	$sRegionGenewiseGFF = "$sWD/genewise_region.gff";


	fnExec("mkdir -p $sWD");

	fnExec('grep -P "'.$sScf.'\t" '. $sMakerModels. ' > '.$sScfGFF);
	fnExec('grep -P "\tmRNA\t" '. $sScfGFF .' > '.$sScfmRNA);

	fnExec("sort -nk4 $sScfmRNA > $sScfSortedmRNA");

	fnExec('echo ">'.$sScf.'" > '.$sScfMaskedFa.'; fastahack -r "'. $sScf .'" "'.$sMaskedGenome.'" >> '.$sScfMaskedFa);
	fnExec('echo ">'.$sScf.'" > '.$sScfFa.'; fastahack -r "'. $sScf .'" "'.$sGenome.'" >> '.$sScfFa);

	fnExec( "sed -n 2p < $sScfFa | wc -m" , $arrRet);
	$nScfLen = intval($arrRet[0]); 

	echo("ScfLen = $nScfLen \n");


	//now, go over the sorted mRNA list (already deduplicated to remove isoforms) to identify fragment models
	//these fragment models are identified by : augustus abinitio model with no scores (quality=unkown)
	//or for models with scores, but completeness is < 50%, and the neighboring fragments have a increasing "missing 5 prime" and decreasing "missing 3 prime" if on the forward strand, and vice versa on reverse. 
	$hScfSortedmRNA = fopen($sScfSortedmRNA , 'r');

	$arrmRNAonScf = array();

	while(false !== ($sLn = fgets($hScfSortedmRNA) )) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		$arrF = explode("\t" , $sLn);
		if (count($arrF) != 9) continue;


		$arrmRNAonScf[] = fnExtractRNAInfo( $arrF);
	}

	$arrmRNAonScf[] = array('qualityunknown' => false, 'id' => 'dummy', 'geneid' =>  'dummy' , 'genewise' => false, 'augustus_evid' => false , 'augustus_abinit' => false, 'evm' => false, 'exonerate' => false , 'completeness' => 1, '5primemiss' => 0, '3primemiss' => 0, 'start' => $nScfLen, 'end' => $nScfLen, 'strand' => '?' );

	$nPrevEnd = 1;

	$nNumGenes = count($arrmRNAonScf);


	$arrProblemRegions = array();
	for($i=0; $i<($nNumGenes-1); $i++) {
		$oCurrGene = $arrmRNAonScf[$i];


		$bReliable = ( $oCurrGene['completeness'] >=0.90 ); 

		

		if ( !$bReliable  ) { // if the current region becomes reliable, strand switches
			//echo("this reliable, prev not reliable, or both are unreliable and on different strands\n");
			$arrCurrProblemRegion['start'] = $nPrevEnd+1;
			$arrCurrProblemRegion['end'] = ($arrmRNAonScf[$i+1]['start'] < $arrmRNAonScf[$i+1]['end'])? $arrmRNAonScf[$i+1]['start']-1 : $arrmRNAonScf[$i+1]['end']-1;
			$arrCurrProblemRegion['gene'] = $oCurrGene;
			$arrProblemRegions[] = $arrCurrProblemRegion;

		} else {
			//leave gene.
		}

		$nPrevEnd = ($oCurrGene['end'] > $oCurrGene['start'] )? $oCurrGene['end'] : $oCurrGene['start'];

	}

	print_r($arrProblemRegions);

	//die();

	if (count($arrProblemRegions) == 0 ) {
		echo("No potential problems found in scaffold $sScf, great.\n");
		continue;
	}


	$sRegionBed = "$sWD/region.bed";
	$sRegionBlastxGFF = "$sWD/blastx_region.gff";

	$sAnnotNotes = "$sWD/annotations.notes";
	$sFixedGFF = "$sWD/fixed.gff";

	$hAnnotNotes = fopen($sAnnotNotes , "w");
	$hFixedGFF = fopen($sFixedGFF , "w");



	$nRegionCount = -1;
	foreach($arrProblemRegions as &$oRegion) {

		if ( !array_key_exists('ensemblorthologs' ,$oRegion['gene'] ) || trim(implode("",$oRegion['gene']['ensemblorthologs'])) == '') continue;
		//no ensembl orthologs were found, skip this.
			echo("Ensembl orthologs:".implode(' ', $oRegion['gene']['ensemblorthologs']).PHP_EOL);
			$sQueryProtFa = "$sWD/queryprot.fa";
			$sQueryProtAlnFa = "$sWD/queryprot.aln.fa";
			$sQueryProtHMM = "$sWD/prot.hmm";
			$sGenewiseOut = "$sWD/genewise.out";
			$sFixedGFFThis = "$sWD/fixed.part.gff";
			$hFixedGFFThis = fopen($sFixedGFFThis ,'w');

			$nRegioStart = ($oRegion['start'] < $oRegion['end'])? $oRegion['start'] : $oRegion['end']; //if is first cluster, region start is start of the whole problem region. else, start from the previous end
			$nRegioEnd =  ($oRegion['start'] > $oRegion['end'])? $oRegion['start'] : $oRegion['end']; //if is first cluster, region end is end of the whole problem region. else, end at next

			//$sSubRegionFA = "$sWD/region.fa";

			$hQueryProtFa = fopen($sQueryProtFa , 'w');
			$nProtCount = 0;
			foreach($oRegion['gene']['ensemblorthologs'] as $sEnsemblID) {
				$sProtName = (array_key_exists($sEnsemblID, $arrEnsemblProteinMap))? $arrEnsemblProteinMap[$sEnsemblID] : '';
				if (trim($sProtName) == '' ) continue;
				$nProtCount++;
				$sProtSeq = $oProteinFasta->GetContig($sProtName);
				fwrite($hQueryProtFa  , ">$sProtName\n$sProtSeq\n");
			}


			$sGenewiseInFile = $sQueryProtHMM;
			$sHmmerParam = "";
			$sAlgParam = "623L";

			if ($nProtCount == 0) continue;
			if ($nProtCount > 1) {
				fnExec("mafft --quiet $sQueryProtFa > $sQueryProtAlnFa");
				//produce hmmer profile:
				fnExec("if [ -e $sQueryProtHMM ];then rm $sQueryProtHMM; fi; $HMMER $sQueryProtHMM $sQueryProtAlnFa");
				$sHmmerParam = "-hmmer";

			} else {
				fnExec("cp $sQueryProtFa  $sQueryProtAlnFa");
				$sGenewiseInFile = $sQueryProtAlnFa;
				$sAlgParam = "623";
			}

			//produce genewise alignment using global alignment mode:
			$sOrientArg = ($oRegion['gene']['strand'] == '+')? '-tfor' : '-trev';
			fnExec("genewise $sOrientArg -silent -init global -alg $sAlgParam -gff -kbyte 250000 -u $nRegioStart -v $nRegioEnd $sHmmerParam $sGenewiseInFile $sScfFa  > $sGenewiseOut");
			$smRNAID = fnConvertGenewise($sGenewiseOut , $hFixedGFFThis, $nRegionCount++);

			if (false === $smRNAID ) { //if genewise failed, keep the original fragments
				echo("Genewise failed to produce any model. Keep the original fragments...\n");
				foreach(array($oRegion['gene']) as &$oGene) {
					fwrite($hAnnotNotes , $oGene['id']."\tkeep\n");
				}
				continue;
			}


			//now add quality annotation for the genewise model.
			$oGenewiseModel =  new MappedCDSGFF();
			$oGenewiseModel->LoadGFF($sFixedGFFThis , $sScfFa);
			$sGenewiseSeq = $oGenewiseModel->ExtractmRNASequence('maker', $smRNAID);

			$arrScores = array();

			foreach($oRegion['gene']['ensemblorthologs'] as $sEnsemblID) {
				$sProtName = (array_key_exists($sEnsemblID, $arrEnsemblProteinMap))? $arrEnsemblProteinMap[$sEnsemblID] : '';
				if (trim($sProtName) == '' ) continue;
				$sProtSeq = $oProteinFasta->GetContig($sProtName);
				$oScore = fnScoreExonerateModel($sGenewiseSeq['AA'] , $sProtSeq , $sWD );
				$arrScores[$oScore['score']] = $oScore;
			}

			ksort($arrScores, SORT_NUMERIC);
			$oHighScore = array_pop($arrScores);

			$bKeepGenewiseModel = false;
			//now validates that the obtained genewise model covers the original fragments:
			foreach( array($oRegion['gene']) as $oGene) {
				//check if the original fragments are overlapping with the genewise model.
				$sOrigFragBed  = "$sWD/orig.frag.bed";
				$hOrigFragBed = fopen($sOrigFragBed , 'w');
				fwrite($hOrigFragBed , "$sScf\t".$oGene['start']."\t".$oGene['end']."\n");

				fnExec("bedtools intersect -a $sOrigFragBed  -b $sFixedGFFThis -f 0.01 -F 0.001 | wc -l", $arrRet);
				if (intval($arrRet[0]) >=1 ) {
					//make sure that the model has a lower score than the new genewise model:
					if ($oGene['score'] < $oHighScore['score']) {
						echo("Replace ".$oGene['id']." with new model\n");
						fwrite($hAnnotNotes , $oGene['id']."\tdelete\n");
						$bKeepGenewiseModel = true;
					} else {
						fwrite($hAnnotNotes , $oGene['id']."\tkeep\n");
					}
				} else {
					fwrite($hAnnotNotes , $oGene['id']."\tkeep\n");
				}
				continue;
				
			}

			if ($bKeepGenewiseModel) {
				$hFixedGFFThis = fopen($sFixedGFFThis , 'r');
				while(false !== ($sLn = fgets($hFixedGFFThis) ) ) {
					$sLn = trim($sLn);
					if ($sLn == '') continue;
					$arrF = explode("\t", $sLn);
					if ($arrF[2] == 'mRNA') {
						$oGene['origannot']['OLDID'] = $oGene['origannot']['ID'];
						$oHighScore['ID'] = $smRNAID;
						$oHighScore['Parent'] = 'gene-'.$smRNAID ;
						$arrF[8] = fnArr2Annot( array_merge( $oGene['origannot'] ,$oHighScore) );
						//set id and parent
						

					}
					if ($arrF[2] == 'gene') {
						$oGene['origannot']['OLDID'] = $oGene['origannot']['ID'];
						$arrF[8] = fnArr2Annot( array_merge( array('ID' => 'gene-'.$smRNAID) , sub_array( $oGene['origannot'], array( 'OLDID', 'spgeneid', 'gene' , 'description', 'ensemblorthologs')) ) );

					}

					fwrite($hFixedGFF , implode("\t", $arrF)."\n");
				}
			}
			
		
	}

	//die();

	echo("Scaffold $sScf done.\n");
	exec("touch $sFinishedFlag");
}

function fnExtractRNAInfo( &$arrF) {
	global $sBlastxGFF, $sWD, $oMakerGFF, $oProteinFasta, $sProtFile;
	$arrAnnot =  MappedCDSGFF::fnParseAnnotation($arrF[8]);

	$oRNASeq = $oMakerGFF->ExtractmRNASequence('maker', $arrAnnot['ID']) ;

	$bQualityUnknown = array_key_exists('quality', $arrAnnot)? ($arrAnnot['quality']=='unkown') : false; 
	$sID = $arrAnnot['ID'];
	$sName = array_key_exists('Name', $arrAnnot)? ($arrAnnot['Name']) : $sID; 
	$sParent = $arrAnnot['Parent'];
	$sSpGeneID = $arrAnnot['spgeneid'];
	$arrEnsemblOrthologs = array_key_exists('ensemblorthologs', $arrAnnot)? explode(',', trim($arrAnnot['ensemblorthologs'])) : array(); 

	$bGenewise = ( strpos( strtolower($sID) , 'genewise') !== false ||  strpos( strtolower($sID) , 'sName') !== false );
	$bAugustusAbinit = ( strpos( strtolower($sID) , 'augustus') === 0 );
	$bAugustus = ( strpos( strtolower($sID) , 'augustus') !== false );
	$bEVM = ( strpos( strtolower($sID) , 'evm') !== false );
	$bExonerate = ( strpos( strtolower($sID) , 'exonerate') !== false );

	$nCompleteness = array_key_exists('completeness', $arrAnnot)? $arrAnnot['completeness'] : -1;
	$nScore = array_key_exists('score', $arrAnnot)? $arrAnnot['score'] : 0;
 	$n5PrimeMiss = array_key_exists('5primemiss', $arrAnnot)? $arrAnnot['5primemiss'] : -1;
 	$n3PrimeMiss = array_key_exists('3primemiss', $arrAnnot)? $arrAnnot['3primemiss'] : -1;

	/*if ($bNewAnnoteAvail) {
		print_r($arrAnnot);
		die();
	}*/

	/*if ($sID == 'maker-PLPv12scf55-augustus-gene-0.121-mRNA-1') {
		print_r($arrAnnot);
		die();
	}*/

	return array('qualityunknown' =>$bQualityUnknown, 'id' => $sID, 'geneid' =>  $sParent , 'genewise' => $bGenewise, 'augustus_evid' => ($bAugustus && (!$bAugustusAbinit) ) , 'augustus_abinit' => $bAugustusAbinit, 'evm' => $bEVM, 'exonerate' =>$bExonerate , 'completeness' => $nCompleteness, '5primemiss' => $n5PrimeMiss, '3primemiss' => $n3PrimeMiss, 'start' => $arrF[3], 'end' => $arrF[4], 'strand' => $arrF[6] , 'name' => $sName, 'score' => $nScore, 'spgeneid' => $sSpGeneID , 'ensemblorthologs' => $arrEnsemblOrthologs, 'origannot' => $arrAnnot);

}

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

function fnCoordStartSort($a, $b) {

	if ($a[0]==$b[0]) return 0;
	return ($a[0]<$b[0])? -1:1;

}

function fnConvertGenewise($sGenewiseOut , $hFixedGFF, $nRegionCount) {
	// These genewise models often contains fragments, because these problematic regions probably have assembly error that prevents the gene model to go through completely.
	// This can be seen in the following example:
/*
PLPv12scf1	GeneWise	match	3710701	3710712	6.13	+	.	PLPv12scf1-genewise-prediction-1
PLPv12scf1	GeneWise	cds	3710701	3710712	0.00	+	0	PLPv12scf1-genewise-prediction-1
PLPv12scf1	GeneWise	match	3711060	3735298	200.46	+	.	PLPv12scf1-genewise-prediction-2
PLPv12scf1	GeneWise	cds	3711060	3711230	0.00	+	0	PLPv12scf1-genewise-prediction-2
PLPv12scf1	GeneWise	intron	3711231	3712665	0.00	+	.	PLPv12scf1-genewise-prediction-2
PLPv12scf1	GeneWise	cds	3712666	3712837	0.00	+	0	PLPv12scf1-genewise-prediction-2
PLPv12scf1	GeneWise	intron	3712838	3733759	0.00	+	.	PLPv12scf1-genewise-prediction-2
PLPv12scf1	GeneWise	cds	3733760	3733830	0.00	+	2	PLPv12scf1-genewise-prediction-2
PLPv12scf1	GeneWise	intron	3733831	3733922	0.00	+	.	PLPv12scf1-genewise-prediction-2
PLPv12scf1	GeneWise	cds	3733923	3733962	0.00	+	0	PLPv12scf1-genewise-prediction-2
PLPv12scf1	GeneWise	intron	3733963	3735239	0.00	+	.	PLPv12scf1-genewise-prediction-2
PLPv12scf1	GeneWise	cds	3735240	3735298	0.00	+	2	PLPv12scf1-genewise-prediction-2
PLPv12scf1	GeneWise	match	3735300	3735878	177.09	+	.	PLPv12scf1-genewise-prediction-3
PLPv12scf1	GeneWise	cds	3735300	3735570	0.00	+	0	PLPv12scf1-genewise-prediction-3
PLPv12scf1	GeneWise	intron	3735571	3735735	0.00	+	.	PLPv12scf1-genewise-prediction-3
PLPv12scf1	GeneWise	cds	3735736	3735878	0.00	+	2	PLPv12scf1-genewise-prediction-3
//

On the revere strand the cds are sorted based on the protein:
PLPv12scf1	GeneWise	match	4145078	4082756	10713.60	-	.	PLPv12scf1-genewise-prediction-1
PLPv12scf1	GeneWise	cds	4145078	4145006	0.00	-	0	PLPv12scf1-genewise-prediction-1
PLPv12scf1	GeneWise	intron	4145005	4142855	0.00	-	.	PLPv12scf1-genewise-prediction-1
PLPv12scf1	GeneWise	cds	4142854	4142793	0.00	-	2	PLPv12scf1-genewise-prediction-1
PLPv12scf1	GeneWise	intron	4142792	4140209	0.00	-	.	PLPv12scf1-genewise-prediction-1
PLPv12scf1	GeneWise	cds	4140208	4140113	0.00	-	0	PLPv12scf1-genewise-prediction-1
PLPv12scf1	GeneWise	intron	4140112	4140017	0.00	-	.	PLPv12scf1-genewise-prediction-1
PLPv12scf1	GeneWise	cds	4140016	4139939	0.00	-	0	PLPv12scf1-genewise-prediction-1
PLPv12scf1	GeneWise	intron	4139938	4139837	0.00	-	.	PLPv12scf1-genewise-prediction-1
PLPv12scf1	GeneWise	cds	4139836	4139744	0.00	-	0	PLPv12scf1-genewise-prediction-1
PLPv12scf1	GeneWise	intron	4139743	4139637	0.00	-	.	PLPv12scf1-genewise-prediction-1
PLPv12scf1	GeneWise	cds	4139636	4139440	0.00	-	0	PLPv12scf1-genewise-prediction-1
PLPv12scf1	GeneWise	intron	4139439	4139305	0.00	-	.	PLPv12scf1-genewise-prediction-1
.....
PLPv12scf1	GeneWise	match	4082750	4082475	209.19	-	.	PLPv12scf1-genewise-prediction-2
PLPv12scf1	GeneWise	cds	4082750	4082475	0.00	-	0	PLPv12scf1-genewise-prediction-2
//

*/
	//in the above example, the genewise alignment is broken up to three "matches". with 3 "predictions".
	//convert this into the following format:
/*
PLPv12scf1	maker	gene	3710701	3735878	.	+	.	ID=PLPv12scf1-frag-merged-clust0
PLPv12scf1	maker	mRNA	3710701	3735878	.	+	.	ID=PLPv12scf1-frag-merged-clust0-mRNA-1;Parent=PLPv12scf1-frag-merged-clust0;Fragments=3;FragmentList=3710701-3710712,3711060-3735298,3735300-3735878
PLPv12scf1	maker	exon	3710701	3710712	0.00	+	.	ID=PLPv12scf1-frag-merged-clust0-mRNA-1-frag-1-1;Parent=PLPv12scf1-frag-merged-clust0-mRNA-1;FragmentID=1
PLPv12scf1	maker	CDS	3710701	3710712	0.00	+	0	ID=PLPv12scf1-frag-merged-clust0-mRNA-1-frag-1-1;Parent=PLPv12scf1-frag-merged-clust0-mRNA-1;FragmentID=1
... etc...
PLPv12scf1	GeneWise	cds	3711060	3711230	0.00	+	0	PLPv12scf1-genewise-prediction-2
PLPv12scf1	GeneWise	intron	3711231	3712665	0.00	+	.	PLPv12scf1-genewise-prediction-2
PLPv12scf1	GeneWise	cds	3712666	3712837	0.00	+	0	PLPv12scf1-genewise-prediction-2
PLPv12scf1	GeneWise	intron	3712838	3733759	0.00	+	.	PLPv12scf1-genewise-prediction-2
PLPv12scf1	GeneWise	cds	3733760	3733830	0.00	+	2	PLPv12scf1-genewise-prediction-2
PLPv12scf1	GeneWise	intron	3733831	3733922	0.00	+	.	PLPv12scf1-genewise-prediction-2
PLPv12scf1	GeneWise	cds	3733923	3733962	0.00	+	0	PLPv12scf1-genewise-prediction-2
PLPv12scf1	GeneWise	intron	3733963	3735239	0.00	+	.	PLPv12scf1-genewise-prediction-2
PLPv12scf1	GeneWise	cds	3735240	3735298	0.00	+	2	PLPv12scf1-genewise-prediction-2
PLPv12scf1	GeneWise	match	3735300	3735878	177.09	+	.	PLPv12scf1-genewise-prediction-3
PLPv12scf1	GeneWise	cds	3735300	3735570	0.00	+	0	PLPv12scf1-genewise-prediction-3
PLPv12scf1	GeneWise	intron	3735571	3735735	0.00	+	.	PLPv12scf1-genewise-prediction-3
PLPv12scf1	GeneWise	cds	3735736	3735878	0.00	+	2	PLPv12scf1-genewise-prediction-3
*/

	$hGenewise = fopen($sGenewiseOut  , 'r');
	$nFragmentCount = -1;
	$arrFrags = array(); //key is fragment count, 0 based index
	$nLeftMostCoord = PHP_INT_MAX;
	$nRightMostCoord = 0;
	$sScf = '';
	$sOrient = "";
	$arrFragCoord = array();

	while(false!==($sLn = fgets($hGenewise) )) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		if ($sLn == '//') break;
		$arrF = explode("\t", $sLn);
		$arrF[1] = 'maker';
		$nGFFStart = ($arrF[3] < $arrF[4])? $arrF[3] : $arrF[4];
		$nGFFEnd  = ($arrF[3] < $arrF[4])? $arrF[4] : $arrF[3];
		$arrF[3] = $nGFFStart ;
		$arrF[4] = $nGFFEnd ;
		if ($arrF[2] == 'match') {
			$sScf = $arrF[0];
			$sOrient = $arrF[6];
			$nFragmentCount++;
			$arrFrags[$nFragmentCount] = array('start'=> $nGFFStart, 'end'=>$nGFFEnd , 'cds'=>array());
			$arrFragCoord[] = $nGFFStart."-".$nGFFEnd ;
			if ($arrF[3] < $nLeftMostCoord) $nLeftMostCoord = $nGFFStart;
			if ($arrF[4] > $nRightMostCoord) $nRightMostCoord = $nGFFEnd ;
			continue;
		}
		if ($arrF[2] == 'cds') {
			$arrFrags[$nFragmentCount]['cds'][] = $arrF;
			continue;
		}


	}

	if ($nFragmentCount == -1) {
		return false;
	}

	$sIDBase = "$sScf-frag-improved-region$nRegionCount";
	$smRNAID = "$sIDBase-mRNA-1";
	fwrite($hFixedGFF , "$sScf\tmaker\tgene\t$nLeftMostCoord\t$nRightMostCoord\t.\t$sOrient\t.\tID=$sIDBase\n");
	fwrite($hFixedGFF , "$sScf\tmaker\tmRNA\t$nLeftMostCoord\t$nRightMostCoord\t.\t$sOrient\t.\tID=$smRNAID;Parent=$sIDBase;Fragments=".count($arrFrags).";FragmentList=".implode(",",$arrFragCoord )."\n");
	for($nFrag=0;$nFrag<count($arrFrags);$nFrag++) {
		$oFrag = $arrFrags[$nFrag];
		for($nCDS=0; $nCDS < count($oFrag['cds']); $nCDS++ ) {
			$oCDS = $oFrag['cds'][$nCDS];
			$oExon = $oCDS;
			$oExon[2] = "exon";
			$oExon[7] = ".";
			$oCDS[2] = "CDS"; 
			$sAnnot = "ID=$smRNAID-frag-$nFrag-$nCDS;Parent=$smRNAID;FragmentID=$nFrag";
			$oExon[8] = $oCDS[8] = $sAnnot;
			fwrite($hFixedGFF , implode("\t", $oExon) . "\n" );
			fwrite($hFixedGFF , implode("\t", $oCDS) . "\n" );
		}
	}
	
	return $smRNAID ;
	
}

function fnGetEnsemblProteinNames($sEnsemblProteins) {
	$hIn = fopen($sEnsemblProteins, 'r');

	$arrRet = array();
	while(false !== ($sLn = fgets($hIn) ) ) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		if ($sLn[0] == '>') {
			$arrF = explode("|" , $sLn);
			if (count($arrF) !=3) continue;
			if (substr($arrF[2], 0, 3) == 'ENS') {
				$arrF1 = explode('_' , $arrF[2]);
				$arrRet[$arrF1[0]] = substr($sLn, 1);
			}
		}
	}

	return $arrRet;
}

function sub_array( $haystack,  $needle)
{
    return array_intersect_key($haystack, array_flip($needle));
}

?>
