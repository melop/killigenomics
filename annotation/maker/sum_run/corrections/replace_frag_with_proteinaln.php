<?php
require_once(dirname(__FILE__) . "/lib.php");
require_once(dirname(__FILE__) . "/scores.lib.php");
$sGBlastPath = "/beegfs/group_dv/software/source/genblast/genblast_v139";
putenv("GBLAST_PATH=$sGBlastPath");

//sometimes, a very long protein model is scattered into multiple pieces
// this is often shown as ab initio predictions, such as many augustus "genes".
// when you look at the underlying exonerate and genewise models
// often a more complete model can be identified.
// This program uses bedtools to find genewise models that overlays multiple augustus ab initio genes (gene id starting with "augustus")
// And replaces the augustus fragments with the genewise / exonerate models

$sProtFile = "all.protein.evidence.fa";
$sMakerModels = "maker.refined.removeisoform.missedaddedback.gff";
$sGenome = "query.fa";
$sMaskedGenome = "query.masked.fa";//"../../../../repeatmasker/PLP_v1.0_pilon/scf.fa.masked";
$sOutDir = "./frag_to_aln";
$sFullMakerGFF = 'genome.all.gff';//"../genome_pilon.all.gff";
$sCandidateModels = "split__candidate.gff";

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


$nScfCount = -1;
//$arrScfs = array("PLPv12scf203" => true);

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

	$sBlastxGFF = "$sWD/blastx.gff";
	$sGenewiseGFF = "$sWD/genewise.gff";
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

	fnExec('grep -P "'.$sScf.'\tblastx:\S+Prot\t" '.$sFullMakerGFF." > $sBlastxGFF ");
	fnExec('grep -P "'.$sScf.'\tpred_gff_genewise" '.$sCandidateModels." > $sGenewiseGFF ");

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
	//now scan the gene

	$nPrevEnd = 1;
	$bPrevReliable = true;
	$nPrevStrand = '+';
	$nPrev5primemiss = 0;
	$nPrev3primemiss = 0;
	$bPrevQualUnknown = false;

	$nNumGenes = count($arrmRNAonScf);

	$arrProblemRegions = array(); //all problem regions
	$arrCurrProblemRegion = array('genes'=> array() ); //genes within this problem region

	for($i=0; $i<$nNumGenes; $i++) {
		$oCurrGene = $arrmRNAonScf[$i];
		//print_r($oCurrGene);

		//print_r($oCurrGene);
		//die();
		//$bReliable = ( $oCurrGene['exonerate'] || $oCurrGene['genewise'] || $oCurrGene['completeness'] >=0.90 ); 

		$bReliable = ( $oCurrGene['completeness'] >=0.90 ); 

		//echo($oCurrGene['id'] ."\n");
			//print_r($oCurrGene);
			//if ($bReliable) {echo("reliable\n");}
			//if ($bPrevReliable) {echo("prevreliable\n");}

		

		if ( ($bReliable && (!$bPrevReliable )) || ((( !$bReliable) && (!$bPrevReliable )) && $nPrevStrand!= $oCurrGene['strand'] ) ) { // if the current region becomes reliable, strand switches
			//echo("this reliable, prev not reliable, or both are unreliable and on different strands\n");
			$arrCurrProblemRegion['start'] = $nPrevEnd+1;
			$arrCurrProblemRegion['end'] = $oCurrGene['start']-1;
			if (count( $arrCurrProblemRegion['genes']) >0 ) $arrProblemRegions[] = $arrCurrProblemRegion; // add the current problem region
			$arrCurrProblemRegion = array('genes'=> array() ); //genes within this problem region
			if ((( !$bReliable) && (!$bPrevReliable )) && $nPrevStrand!= $oCurrGene['strand'] ) {
				$arrCurrProblemRegion['genes'][] = $oCurrGene;
			}

		} else if (!$bReliable) {
			//echo("this unreliable\n");
			if ( (!$bPrevReliable) && ($nPrevStrand == $oCurrGene['strand'] ) ) { // join to current region
				//echo("\tprev unreliable and on same strand\n");
				$bDoJoin = true;
				if ( (!$bPrevQualUnknown) && (!$oCurrGene['qualityunknown'] ) ) { // if 5' and 3' missing is known from the previous one check if decrease and increase in 5' and 3' monotonous
					$bIsMonotonousMissing = ($oCurrGene['strand']=='+')? ( $nPrev5primemiss < $oCurrGene['5primemiss'] && $nPrev3primemiss > $oCurrGene['3primemiss'] ) :
											     ( $nPrev5primemiss > $oCurrGene['5primemiss'] && $nPrev3primemiss < $oCurrGene['3primemiss'] ) ;
					$bDoJoin = $bIsMonotonousMissing;

				}

				if ($bDoJoin) {
					//echo("\tjoin\n");
					//print_r($arrCurrProblemRegion['genes']);
					$arrCurrProblemRegion['genes'][] = $oCurrGene;
				} else {
					//echo("\tsave prev. putcurrent in new\n");
					$arrCurrProblemRegion['start'] = $nPrevEnd+1;
					$arrCurrProblemRegion['end'] = $oCurrGene['start']-1;
					if (count( $arrCurrProblemRegion['genes']) >0 ) $arrProblemRegions[] = $arrCurrProblemRegion; // add the current problem region
					$arrCurrProblemRegion = array('genes'=> array( 0 => $oCurrGene) ); //genes within this problem region
					$nPrevEnd = ($i>0)? ($arrmRNAonScf[$i-1]['end']) : 1;
				}
			} else if ($bPrevReliable) { //previous is reliable, this is the first one unreliable
				//echo("\tprev reliable\n");
				$arrCurrProblemRegion['genes'][0] = $oCurrGene;
				//print_r($arrCurrProblemRegion);
				//die();
			} else { //both unreliable, but on different strands.
				//echo("\tprev unreliable, on different strands. tsave prev. putcurrent in new\n");
				$arrCurrProblemRegion['start'] = $nPrevEnd+1;
				$arrCurrProblemRegion['end'] = $oCurrGene['start']-1;
				if (count( $arrCurrProblemRegion['genes']) >0 ) $arrProblemRegions[] = $arrCurrProblemRegion; // add the current problem region
				$arrCurrProblemRegion = array('genes'=> array( 0 => $oCurrGene ) ); //genes within this problem region
				$nPrevEnd = ($i>0)? ($arrmRNAonScf[$i-1]['end']) : 1;
			}
		}

		if ($bReliable ) { 
			$nPrevEnd = $oCurrGene['end'];
		}

		if ( (!$bReliable) && ($nPrevStrand!= $oCurrGene['strand']) ) {
			$nPrevEnd = ($i>0)? ($arrmRNAonScf[$i-1]['end']) : 1;
		}

		$bPrevReliable = $bReliable;
		$nPrevStrand = $oCurrGene['strand'];
		$nPrev5primemiss = $oCurrGene['5primemiss'];
		$nPrev3primemiss = $oCurrGene['3primemiss'];
		$bPrevQualUnknown = $oCurrGene['qualityunknown'];

	}

	//print_r($arrProblemRegions);

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
		$nRegionCount++;
		$bSplitAlready = false;
		$sSplitID = "";
		/*foreach($oRegion['genes'] as &$oGene) {
			if (strpos($oGene['id'], 'SplitChimeric_') === 0) {
				$bSplitAlready = true;
				$sSplitID = $oGene['id'];
				fwrite($hAnnotNotes , "$sSplitID\tkeep\n");
				continue;
			}
		}
		if ($bSplitAlready) {
			echo("Skip region because it contains genes previously known to be misjoined." . $sSplitID.PHP_EOL);
			continue;
		}*/

		//get sub region blastx gff:
		$hRegionBed = fopen($sRegionBed, "w");
		fwrite($hRegionBed , "$sScf\t".($oRegion['start']-1)."\t".$oRegion['end'].PHP_EOL);
		fclose($hRegionBed);

		fnExec("bedtools intersect -a $sBlastxGFF -b $sRegionBed -f 0.5 -F 0.01 | grep -P ".'"\tprotein_match\t"'." > $sRegionBlastxGFF");
		//fnExec("bedtools intersect -a $sGenewiseGFF -b $sRegionBed -f 0.01 -F 0.8 | grep -P ".'"\tmRNA\t"'." > $sRegionGenewiseGFF");
		
		if (filesize($sRegionBlastxGFF) == 0 ) { //&& filesize($sRegionGenewiseGFF)==0 ) {
			echo("Nothing aligns in this region $sScf:".$oRegion['start']."..".$oRegion['end']." skip...\n");
			foreach($oRegion['genes'] as &$oGene) {

				fwrite($hAnnotNotes , $oGene['id']."\tkeep\t".$oGene['newannotation']."\n");
				continue;
				
			}

			continue;
		}

		print_r($oRegion);

		//cluster the blastx hit's start and end positions. make sure that they overlap with the two fragments. This suggests that the two fragments came from one gene.
		$hRegionBlastxGFF = fopen($sRegionBlastxGFF , 'r');

		$sClusterInput = "$sWD/clusterin.txt";
		$sClusterOutput = "$sWD/clusterout.txt";
		$sClusterMedOutput = "$sWD/clustermedout.txt";
		$hClusterInput = fopen($sClusterInput ,'w');
		fwrite($hClusterInput , "Start\tEnd\tID\n");

	
		$nWrittenEvidences = 0;

		while(false !== ($sLn = fgets($hRegionBlastxGFF ) ) ) {
			$sLn = trim($sLn);
			if ($sLn == '') continue;
			$arrF = explode("\t", $sLn);
			$nStart = $arrF[3];
			$nEnd = $arrF[4];
			$sStrand = $arrF[6];

			if ($sStrand != $oRegion['genes'][0]['strand'] ) {
				continue; // match not on the same strand, skip.
			}

			$arrAnnot =  MappedCDSGFF::fnParseAnnotation($arrF[8]);
			$sProtName = $arrAnnot['Name'];
			fwrite($hClusterInput , "$nStart\t$nEnd\t$sProtName\n");
			$nWrittenEvidences++;
		}

		if ($nWrittenEvidences == 0) {
			echo("Nothing aligns in this region $sScf:".$oRegion['start']."..".$oRegion['end']." skip...\n");
			foreach($oRegion['genes'] as &$oGene) {

				fwrite($hAnnotNotes , $oGene['id']."\tkeep\t".$oGene['newannotation']."\n");
				continue;
				
			}

			continue;
		}

		 fnExec("Rscript coordcluster.R $sClusterInput $sClusterOutput $sClusterMedOutput", $arrRetDummy, $nRetVal);
		if ($nRetVal != 0) {
			echo(" warning: Rscript error!\n");
		}

		$hClusterTab = fopen($sClusterOutput , 'r');
		$hClusterMedTab = fopen($sClusterMedOutput , 'r'); //calculated median start and end positions of the cluster

		$arrClusters = array();   
		$arrProteinByCluster = array();

		while( false !== ($sLn = fgets($hClusterMedTab) ) ) {
		        $sLn = trim($sLn);
		        if ($sLn == '') continue;
		        $arrF = explode("\t" , $sLn);
		        if ($arrF[0] == 'clu') continue;
		        $arrClusters[] = array($arrF[1], $arrF[2], $arrF[0] );
		        $arrProteinByCluster[$arrF[0]] = array();
			//echo($arrF[0]);
		}

		while ( false !== ($sLn = fgets($hClusterTab) )  ) {
		        $sLn = trim($sLn);
		        if ($sLn == '') continue;
		        $arrF = explode("\t" , $sLn);
		        if ($arrF[0] == 'Start') continue;

			if ($arrF[4] == 0) continue;
			if (!array_key_exists($arrF[3] , $arrProteinByCluster) ) {
				continue;
			}
			$arrProteinByCluster[$arrF[3]][] = $arrF[2];
			//echo($arrF[3]);
		}

		$nSubFragments = count($arrClusters);
		echo("$nSubFragments clusters of blastx alignments found in the region.\n");
		//print_r($arrProteinByCluster);
		usort($arrClusters , "fnCoordStartSort");
		for($nClust=0; $nClust < count($arrClusters); $nClust++) {
			//process each cluster
			$oClust = $arrClusters[$nClust];
			$nClustID = $oClust[2];
			
			$nRegioStart = ($nClust==0)? $oRegion['start'] : $arrClusters[$nClust-1][1]; //if is first cluster, region start is start of the whole problem region. else, start from the previous end
			$nRegioEnd = ($nClust== (count($arrClusters)-1) )? $oRegion['end'] : $arrClusters[$nClust+1][0]; //if is first cluster, region end is end of the whole problem region. else, end at next start
			$sQueryProtFa = "$sWD/queryprot.fa";
			$sQueryProtAlnFa = "$sWD/queryprot.aln.fa";
			$sQueryProtHMM = "$sWD/prot.hmm";
			$sGenewiseOut = "$sWD/genewise.out";
			$sFixedGFFThis = "$sWD/fixed.part.gff";
			$hFixedGFFThis = fopen($sFixedGFFThis ,'w');

			//$sSubRegionFA = "$sWD/region.fa";

			$hQueryProtFa = fopen($sQueryProtFa , 'w');
			$nProtCount = 0;
			foreach($arrProteinByCluster[$nClustID] as $sProtName) {
				$nProtCount++;
				$sProtSeq = $oProteinFasta->GetContig($sProtName);
				fwrite($hQueryProtFa  , ">$sProtName\n$sProtSeq\n");
			}
			//fnExec('echo ">'."$sScf"."_$nRegioStart"."_$nRegioEnd".'"'." > $sSubRegionFA; fastahack -r $sScf:$nRegioStart..$nRegioEnd $sGenome >> $sSubRegionFA");
			//use mafft to align the proteins

			$sGenewiseInFile = $sQueryProtHMM;
			$sHmmerParam = "";
			$sAlgParam = "623L";

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
			$sOrientArg = ($oRegion['genes'][0]['strand'] == '+')? '-tfor' : '-trev';
			fnExec("genewise $sOrientArg -silent -init global -alg $sAlgParam -gff -kbyte 250000 -u $nRegioStart -v $nRegioEnd $sHmmerParam $sGenewiseInFile $sScfFa  > $sGenewiseOut");
			$smRNAID = fnConvertGenewise($sGenewiseOut , $hFixedGFFThis, $nRegionCount, $nClust);

			if (false === $smRNAID ) { //if genewise failed, keep the original fragments
				echo("Genewise failed to produce any model. Keep the original fragments...\n");
				foreach($oRegion['genes'] as &$oGene) {
					fwrite($hAnnotNotes , $oGene['id']."\tkeep\t".$oGene['newannotation']."\n");
				}
				continue;
			}


			//now add quality annotation for the genewise model.
			$oGenewiseModel =  new MappedCDSGFF();
			$oGenewiseModel->LoadGFF($sFixedGFFThis , $sScfFa);
			$sGenewiseSeq = $oGenewiseModel->ExtractmRNASequence('maker', $smRNAID);

			$arrScores = array();

			foreach($arrProteinByCluster[$nClustID] as $sProtName) {
				$sProtSeq = $oProteinFasta->GetContig($sProtName);
				$oScore = fnScoreExonerateModel($sGenewiseSeq['AA'] , $sProtSeq , $sWD );
				$arrScores[$oScore['score']] = $oScore;
			}

			ksort($arrScores, SORT_NUMERIC);
			$oHighScore = array_pop($arrScores);

			$bKeepGenewiseModel = false;
			//now validates that the obtained genewise model covers the original fragments:
			foreach($oRegion['genes'] as &$oGene) {
				//check if the original fragments are overlapping with the genewise model.
				$sOrigFragBed  = "$sWD/orig.frag.bed";
				$hOrigFragBed = fopen($sOrigFragBed , 'w');
				fwrite($hOrigFragBed , "$sScf\t".$oGene['start']."\t".$oGene['end']."\n");

				fnExec("bedtools intersect -a $sOrigFragBed  -b $sFixedGFFThis -f 0.5 -F 0.001 | wc -l", $arrRet);
				if (intval($arrRet[0]) >=1 ) {
					//make sure that the model has a lower score than the new genewise model:
					if ($oGene['score'] < $oHighScore['score']) {
						echo("Replace ".$oGene['id']." with new model\n");
						fwrite($hAnnotNotes , $oGene['id']."\tdelete\n");
						$bKeepGenewiseModel = true;
					} else {
						fwrite($hAnnotNotes , $oGene['id']."\tkeep\t".$oGene['newannotation']."\n");
					}
				} else {
					fwrite($hAnnotNotes , $oGene['id']."\tkeep\t".$oGene['newannotation']."\n");
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
						$arrF[8] .= ";".fnArr2Annot($oHighScore);

					}
					fwrite($hFixedGFF , implode("\t", $arrF)."\n");
				}
			}
			
		}



		//if ($oRegion['genes'][0]['strand'] == '-') {
			//die();
		//}

	}

	//die();

	echo("Scaffold $sScf done.\n");
	exec("touch $sFinishedFlag");
}

function fnExtractRNAInfo( &$arrF) {
	global $sBlastxGFF, $sWD, $oMakerGFF, $oProteinFasta, $sProtFile;
	$arrAnnot =  MappedCDSGFF::fnParseAnnotation($arrF[8]);

	$oRNASeq = $oMakerGFF->ExtractmRNASequence('maker', $arrAnnot['ID']) ;

	//first, attemp to obtain a score if score does not already exist!
	$bNewAnnoteAvail = false;
	$bQualityUnknown = array_key_exists('quality', $arrAnnot)? ($arrAnnot['quality']=='unkown') : false; 
	if ($bQualityUnknown) {
		$sRegionBed = "$sWD/extractrnainfo.region.bed";
		$hRegionBed = fopen($sRegionBed, "w");
		$sEvidBlastxGFF = "$sWD/extractrnainfo.blastx.gff";
		fwrite($hRegionBed , "$arrF[0]\t".($arrF[3]-1)."\t".$arrF[4].PHP_EOL);
		fclose($hRegionBed);

		fnExec("bedtools intersect -a $sBlastxGFF -b $sRegionBed -f 0.01 -F 0.8 | grep -P ".'"\tprotein_match\t"'." > $sEvidBlastxGFF");
		//fnExec("bedtools intersect -a $sGenewiseGFF -b $sRegionBed -f 0.01 -F 0.8 | grep -P ".'"\tmRNA\t"'." > $sRegionGenewiseGFF");
		$arrEvidenceProtNames = array();
		$hEvidBlastxGFF = fopen($sEvidBlastxGFF, "r");		
		while( false !== ($sLn = fgets($hEvidBlastxGFF) ) ) {
			$arrFx = explode("\t" , trim($sLn));
			if (count($arrFx) != 9) continue;
			$arrAnnotEvidProt = MappedCDSGFF::fnParseAnnotation($arrFx[8]);
			if (array_key_exists('Name', $arrAnnotEvidProt) ) $arrEvidenceProtNames[] = $arrAnnotEvidProt['Name'];
		}

		if (count($arrEvidenceProtNames) ==0 ) {//no blast, try redoing with blast p.
			$sTmpFasta = "$sWD/geneprot.fa";
			$hTmpFasta = fopen($sTmpFasta , 'w');
			$sBlastOut = "$sWD/extractrnainfo.blastp.out";
			fwrite($hTmpFasta , ">prot\n".$oRNASeq['AA']."\n");
			fnExec("blastp -task blastp -db $sProtFile -query $sTmpFasta  -evalue 1e-10 -out $sBlastOut -outfmt '7 sseqid evalue'");
			$hBlastOut = fopen($sBlastOut , "r");
			$nBlastLimit = 10;
			$nBlastCount = 0;
			while( false !== ($sLn = fgets($hBlastOut) ) ) {
				if ($sLn[0] == '#') continue;
				$arrFx = explode("\t" , trim($sLn));
				if (count($arrFx) != 2) continue;
				$arrEvidenceProtNames[] = $arrFx[0];
				//append this protein name into the blastx gff
				$arrToWrite = $arrF;
				$arrToWrite[1] = "blastp:Prot";
				$arrToWrite[2] = "protein_match";
 				$arrToWrite[8] = "Name=".$arrFx[0];
				$hBlastxGFF = fopen($sBlastxGFF , 'a'); //append
				fwrite($hBlastxGFF , implode("\t", $arrToWrite) . "\n" );
				if ( $nBlastLimit <= ++$nBlastCount ) break;
			}

		}

		//now try to annotate
		$arrScores = array();

		foreach($arrEvidenceProtNames as $sProtName) {
			//echo("Evidence: $sProtName");
			if (trim($sProtName) == '') continue;
			$sProtSeq = $oProteinFasta->GetContig($sProtName);
			//echo("Evidence: $sProtSeq");
			$oScore = fnScoreExonerateModel($oRNASeq['AA'] , $sProtSeq , $sWD );
			$arrScores[$oScore['score']] = $oScore;
		}

		if (count($arrScores) > 0) {
			ksort($arrScores, SORT_NUMERIC);
			$oHighScore = array_pop($arrScores);
			$oHighScore['ID'] = $arrAnnot['ID'];
			$oHighScore['Parent'] = $arrAnnot['Parent'];
			$arrAnnot = $oHighScore;
			$bNewAnnoteAvail = true;
		}


	}

	$bQualityUnknown = array_key_exists('quality', $arrAnnot)? ($arrAnnot['quality']=='unkown') : false; 
	$sID = $arrAnnot['ID'];
	$sName = array_key_exists('Name', $arrAnnot)? ($arrAnnot['Name']) : $sID; 
	$sParent = $arrAnnot['Parent'];

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

	return array('qualityunknown' =>$bQualityUnknown, 'id' => $sID, 'geneid' =>  $sParent , 'genewise' => $bGenewise, 'augustus_evid' => ($bAugustus && (!$bAugustusAbinit) ) , 'augustus_abinit' => $bAugustusAbinit, 'evm' => $bEVM, 'exonerate' =>$bExonerate , 'completeness' => $nCompleteness, '5primemiss' => $n5PrimeMiss, '3primemiss' => $n3PrimeMiss, 'start' => $arrF[3], 'end' => $arrF[4], 'strand' => $arrF[6] , 'name' => $sName, 'score' => $nScore, 'newannotationavail' => $bNewAnnoteAvail , 'newannotation' => ( $bNewAnnoteAvail)?  fnArr2Annot($arrAnnot): '.' );

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

function fnConvertGenewise($sGenewiseOut , $hFixedGFF, $nRegionCount, $nClust) {
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

	$sIDBase = "$sScf-frag-merged-region$nRegionCount-clust$nClust";
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


?>
