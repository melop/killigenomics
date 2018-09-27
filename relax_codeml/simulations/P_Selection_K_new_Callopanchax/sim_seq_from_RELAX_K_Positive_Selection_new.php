<?php
//simulate sequence alignments from codeml output (branch-site tests)
//the temporary folder must be saved for this to work. in particular, codeml.3.out and in.fas need to be present in the temp folder.
//the missing data pattern is directly copied from in.fas (with N's and gaps copied to the simulated dataset).


### Change testing HERE ###
$TestRelaxation = false;
$TestPositivSelection = true;
### Change testing HERE ###

$nForegroundForBG0 = 0.2108875;
$nMaxPosOmega = 100;

$cmd='Rscript /beegfs/group_dv/home/LIasi/Simulation_Tools/Get_Omega_K.R ';
// path to table
$RELAXRet= "/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/relax_46_spp/rerun_omega0/sum_Callopanchax.txt";
$sCodemlRet = "/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/codeml_46_spp/sum_Callopanchax.txt";
// variable for creating temporary folder/directory
$sTmp = "tmp/".rand();
$nRandSeed  = 13888;
// call analysis programm
$sEvolver = "/beegfs/group_dv/software/source/paml4.8/src/evolverNSbranchsites";
$sOutFasta = "ret_sim.fasta";
$sOutGapFasta = "ret_sim_gapped.fasta";
$sIsForgroundClade= "Callopanchax";
$processONE=null;
$pipesONE	= array();
$descONE 	= array(
		0	=> array('pipe','r'),
		1	=> array('pipe','w'),
		2	=> array('pipe','w')
	);
$arrAlreadyDone = fnGetPrevResults($sOutFasta, $sOutGapFasta);
$process=null;
$pipes	= array();
$desc 	= array(
		0	=> array('pipe','r'),
		1	=> array('pipe','w'),
		2	=> array('pipe','w')
	);

$arrDrawKParameter= array(1,2,3);
//print_r($arrDrawKParameter);
//print_r($arrOmegaFoldchangeParameter);								
fnOpenfnDrawOmegasK();

// the function opens the Set_Omega.R script in the background
function fnOpenfnDrawOmegasK() {

		global $process, $desc, $pipes, $cmd;
		
		$process = proc_open($cmd, $desc, $pipes);
}
$hOutFasta = fopen($sOutFasta, 'a');
$hOutGapFasta = fopen($sOutGapFasta, 'a');

// activates the function to open the Set_Omega.R script
fnOpenSetOmegaR();

// the function opens the Set_Omega.R script in the background
function fnOpenSetOmegaR() {

		global $processONE, $descONE, $pipesONE;
		
		$processONE = proc_open('Rscript /beegfs/group_dv/home/LIasi/Simulation_Tools/Set_Omega_new.R ', $descONE, $pipesONE);
}

// read in RELAX results (which are cat to in one file with /n between the groups
// read therefore line by line
function fnCallRELAX ($sRELAXRet) {
	

	$arrRet = array();
	
		if (!file_exists($sRELAXRet) ) {
				die("RELAX file not found\n");
			}
		
		$hRELAXRet = fopen($sRELAXRet , 'r');

	while( false !== ($sLn = fgets($hRELAXRet) )) {
		$sLn = trim($sLn);
		

		$arrRELAXRet = explode("\t", $sLn);
		if (count($arrRELAXRet) < 21 ) {
			continue;
		}
		
		


		$arrRet[$arrRELAXRet[0]] =  array_slice($arrRELAXRet,10,12);
	}
	//print_r($arrRet);
	//die();
	return $arrRet;
}
// this function calls a R programm drawing foldchanges for the omega out of a distribution
function fnDrawOmegasK(){
	global $process,$desc,$pipes,$arrDrawKParameter;
	
	$arrRet = array();
	foreach($arrDrawKParameter as $nIndex => $OmegaK){
	//implode joints the given array elements as a string with a "glue" together 
	//($OmegaClassTable[0]= foreground and background values for a specific OmegaClass) by tab
	//echo("write to R\n");
	fwrite($pipes[0], $OmegaK."\n");
	//echo("written\nget out\n");

	$randomK= floatval(fgets($pipes[1]));
	//echo("got results from R\n");

	$arrRet[] = $randomK;
	}
	
	return $arrRet;
}
// this function checks the constrains for the omegas
function fnCheckOmegaConstraints($arrOmegasRELAX) {
	return fnOmega_Zero_One_Two($arrOmegasRELAX[0][0],$arrOmegasRELAX[1][0],$arrOmegasRELAX[2][0]) && fnOmega_Zero_One_Two($arrOmegasRELAX[0][1],$arrOmegasRELAX[1][1],$arrOmegasRELAX[2][1]);
	}

function fnOmega_Zero_One_Two($Omega_Zero,$Omega_One,$Omega_Two) { 
			if (
			$Omega_One <= 1
			&&
			$Omega_Zero <= $Omega_One 
			&&
			$Omega_Two >= 1 )	{
			return true;
			} 
			
			return false;
		}

//read in codeml output 
$hCodemlRet = fopen($sCodemlRet , 'r');
$nGroupCount = 0;
$arrF = array();
	// reads file line wise in and removes white space
	$arrRELAXOutput = fnCallRELAX($RELAXRet);
	//echo("Here \n");
	//print_r($arrRELAXOutput);
	while( false !== ($sLn = fgets($hCodemlRet) )) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		$nGroupCount++;
		
// splitts the strings in $sLn by tab and returns all the obtained substrings in an array
		$arrF = explode("\t", $sLn);
			
	//print_r($arrF);

	
		$sGroup = $arrF[0];
		//print_r($sGroup);

		if (!array_key_exists($sGroup , $arrRELAXOutput) ) {
    	echo "ERROR: no matching RELAX/Codeml files for $sGroup \n ";
    	continue;
		}

		//print_r($arrRELAXOutput[0]);
		$arrRelaxRecord = $arrRELAXOutput[$sGroup];
		
		if ($arrRelaxRecord[2] > 1 || $arrRelaxRecord[0] > $arrRelaxRecord[2] || $arrRelaxRecord[4] < 1 )	{
			die ("Omega values of $arrRELAXOutput does not fit the constrains \n");
			} 
			

		//print_r($arrRelaxRecord);

		$sCodemlPath =  $arrF[7];
		//print_r($sCodemlPath);
		$sInFasta = "$sCodemlPath"."in.fas";
		$sCodeml3Out = "$sCodemlPath"."codeml.3.out";
		print_r($sGroup."\n");
		//print_r($sInFasta."\n");
		//print_r($sCodeml3Out."\n");
		if (!file_exists($sInFasta) ) {
			die("$sInFasta not found\n");
		}

		if (!file_exists($sCodeml3Out) ) {
			die("$sCodeml3Out not found\n");
		}

		$arrSimAln = fnSim($sGroup, $sInFasta, $sCodeml3Out, $arrRelaxRecord);
		if($arrSimAln===false){
		continue;
		}

		if ($arrSimAln === false) {
			echo("$sGroup skipped, simulation error.\n");
			continue;
		}

		foreach($arrSimAln['simaln'] as $sTax => $sSeq ) {
			//>Ortho:0;Sp:AAU;MappedToRef:AAU;SpGeneId:029799;GroupId:Group_0_0
			fwrite($hOutFasta, ">Ortho:$nGroupCount;Sp:$sTax;MappedToRef:$sTax;SpGeneId:000000;GroupId:$sGroup\n");
			fwrite($hOutGapFasta, ">Ortho:$nGroupCount;Sp:$sTax;MappedToRef:$sTax;SpGeneId:000000;GroupId:$sGroup\n");
			fwrite($hOutFasta , $sSeq."\n");
			fwrite($hOutGapFasta , $arrSimAln['gappedaln'][$sTax]."\n");
		}
}

function fnMarkTreeWithOmegas($sTree, $arrOmegasRELAX) {
	global $processONE,$descONE,$pipesONE, $sIsForgroundClade;
	$arrRet = array();
	
	foreach($arrOmegasRELAX as $nIndex => $OmegaClass){
	//implode joints the given array elements as a string with a "glue" together 
	//($OmegaClassTable[0]= foreground and background values for a specific OmegaClass) by tab
	//echo("write to R\n");
	fwrite($pipesONE[0], $sTree."\t".$sIsForgroundClade."\t".implode("\t" , $OmegaClass)."\n");
	//echo("written\nget out\n");

	$TreeWithOmegas= fgets($pipesONE[1]);
	//echo("got results from R\n");

	$arrRet[] = $TreeWithOmegas;
	
	}
	
	return $arrRet;
}

function fnSim($sGroup, $sInFasta, $sCodeml3Out, $arrRelaxRecord) {
	global $sTmp, $nRandSeed , $sUseOmegaType, $sEvolver;
	//call external programm to delet the $sTmp directory
	#exec("rm $sTmp/*");
	// execute an external programm, which in this case crates the directory $sTmp
$sTmp = "tmp/".rand();

while(file_exists($sTmp)) {
	$sTmp = "tmp/".rand();
}
echo("tmpfolder $sTmp\n");
exec("mkdir -p $sTmp");
	$arrNGapPattern = fnFindGapPattern($sInFasta);
	//print_r($arrNGapPattern);
	//echo($sInFasta."\n");
	$sMCCodon = "$sTmp/MCcodonNSbranchsites.dat";
	$hMCCodon = fopen($sMCCodon , 'w');
	fwrite($hMCCodon , "0          * 0,1:seqs or patters in paml format (mc.paml); 2:paup format (mc.nex)\n");
	fwrite($hMCCodon , "$nRandSeed       * random number seed (odd number)\n");
	fwrite($hMCCodon , $arrNGapPattern['ntax']." ".intval($arrNGapPattern['seqlen']/3)." 1   * <# seqs>  <# codons>  <# replicates>\n\n");
	
	$arrSimInfo = fnCollectSimInfo($sCodeml3Out,$arrRelaxRecord, $sGroup);
	if ($arrSimInfo ===false)	{
	return false;
		}
	
	
	//print_r($arrSimInfo);

	//print_r($arrSimInfo);
	fwrite($hMCCodon , $arrSimInfo['treelen']."           * <tree length; see note below; use -1 if tree has absolute branch lengths>\n\n");
	fwrite($hMCCodon , $arrSimInfo['tree']."\n\n");
	fwrite($hMCCodon , count($arrSimInfo['OmegaProportions'])."            * number of site classes, followed by frequencies\n");
	foreach($arrSimInfo['OmegaProportions']as $nCol => $nVal) {
		fwrite($hMCCodon , "  " . number_format($nVal, 6));
		}
	fwrite($hMCCodon , "\n\n");
	fwrite($hMCCodon ,$arrSimInfo['OmegaTrees'][0]."\n");
	fwrite($hMCCodon ,$arrSimInfo['OmegaTrees'][1]."\n");
	fwrite($hMCCodon ,$arrSimInfo['OmegaTrees'][2]."\n");
	
	


	fwrite($hMCCodon , "\n\n".$arrSimInfo['kappa']."     * kappa\n\n");

	foreach($arrSimInfo['codonfreq'] as $nRow => $arrRow) {
		foreach($arrRow as $nCol => $nVal) {
			if ($nCol >0) {
				fwrite($hMCCodon , " ");
			}
			fwrite($hMCCodon , number_format($nVal, 6) );
		}
		fwrite($hMCCodon , "\n" );
	}

	fwrite($hMCCodon , "\n0    * genetic code (0:universal; 1:mammalian mt; 2-10:see below)\n\n//end of file\n" );
	exec("WD=`pwd`; cd $sTmp; $sEvolver; cd \$WD;");

	//added by ray:
	//check the site category files, and make sure that positive selected sites "3" are included. if not, return false;
	$hTestOutput = popen("grep '3' $sTmp/siterates.txt | wc -l", 'r');
	$nLnCount = intval(fgets($hTestOutput));	
	pclose($hTestOutput );
	if ($nLnCount == 0) {
		echo("$sGroup excluded because no positive sites were produced by evolver.\n");
		return false;
	} else {
		echo("$sGroup , $nLnCount sites simulated to be under positive selection\n");
	}
	//now parse the output file
	$hMCPAML = fopen("$sTmp/mc.paml" , 'r');

	$arrSimAlns = fnGetSimAln($hMCPAML , $arrNGapPattern);

	return $arrSimAlns;
}

function fnCollectSimInfo($sCodeml3Out, $arrRelaxRecord, $sGroup) {
	//find codon frequencies from output:
	global $TestRelaxation, $TestPositivSelection, $nForegroundForBG0, $nMaxPosOmega;
	
	$h1 = popen("grep -A21 \"Sums of codon usage counts\" $sCodeml3Out" , 'r');

	$arrCodonFreq = array();
	$nCodonCountSum = 0;
	while( false !== ($sLn = fgets($h1) ) ) {
		preg_match_all("/\d+/" , $sLn , $arrM);
		//print_r($arrM);
		if (count($arrM[0]) != 4 ) continue;
		$arrCodonFreq[] = $arrM[0];
		$nCodonCountSum += array_sum($arrM[0]);
	}

	//print_r($arrCodonFreq);
	if (count($arrCodonFreq)!=16 ) {
		die("Codon frequency table in $sCodeml3Out is of unexpected dimensions.\n");
	}



	foreach($arrCodonFreq as $nRow => $arrRow) {
		foreach($arrRow as $nCol => $nFreq) {
			$arrCodonFreq[$nRow][$nCol] = $nFreq/$nCodonCountSum;
		}
	}

	pclose($h1);

	//print_r($arrCodonFreq);

	//read the tree:
	$h1 = popen("grep -A4 \"tree length =\" $sCodeml3Out" , 'r');
	$nLn = 0;
	$nTreeLen = 0;
	$sTree = '';
	while( false !== ($sLn = fgets($h1) ) ) {
	/*	if (empty($nLn)){
		return(false);
		}
	*/
		$nLn++;
		if ($nLn == 1) {
			$arrF = explode("=" , $sLn);
			if (count($arrF) !=2) {
				die("Cannot parse tree length\n");
			}

			$nTreeLen = floatval(trim($arrF[1]));
		}

		if ($nLn == 5) {
			$sTree = trim($sLn);

		} 

		
	}

	pclose($h1);
	
	if ($sTree == ''){
		echo("$sGroup has not tree \n");
		return false ;
	}

	//read kappa
	$h1 = popen("grep \"kappa\" $sCodeml3Out" , 'r');
	$nKappa = 0;
	while( false !== ($sLn = fgets($h1) ) ) {
		$arrF = explode("=" , $sLn);
		if (count($arrF) !=2) {
			die("Cannot parse kappa\n");
		}

		$nKappa = floatval(trim($arrF[1]));

	}

	pclose($h1);
	$arrOmegasRELAX= array();
	do {
	
	$arrDrawnOmegasK=fnDrawOmegasK();
	//print_r($arrRelaxRecord);
	//print_r($arrDrawnOmegasK);
	//print_r($arrDrawnOmegasFoldchanges);
if($TestRelaxation){
		/*if($arrRelaxRecord[2]==1){
		$arrOmegasRELAX= array(
					array($arrRelaxRecord[0]/$arrDrawnOmegasFoldchanges[0],$arrRelaxRecord[0]),
					array($arrRelaxRecord[2],$arrRelaxRecord[2]),
					array($arrRelaxRecord[4]/1,$arrRelaxRecord[4]));
					} else {
					*/
	$arrOmegasRELAX= array(
					array( pow($arrRelaxRecord[0] , $arrDrawnOmegasK[0]),$arrRelaxRecord[0]),
					array( pow($arrRelaxRecord[2] , $arrDrawnOmegasK[1]),$arrRelaxRecord[2]),
					array($arrRelaxRecord[4],$arrRelaxRecord[4]));
					#}

			if ($arrRelaxRecord[0] == 0) { // if the background is 0, set foreground to 0.2108875 in case of 
				$arrOmegasRELAX[0][0] = $nForegroundForBG0;
			} 
					
				}
if($TestPositivSelection){
			$arrOmegasRELAX= array(
					array($arrRelaxRecord[0],$arrRelaxRecord[0]),
					array($arrRelaxRecord[2],$arrRelaxRecord[2]),
					array( pow($arrRelaxRecord[4] , $arrDrawnOmegasK[2]) ,$arrRelaxRecord[4]));
					}

			if ($arrOmegasRELAX[2][0] > $nMaxPosOmega) {
				if ( $arrOmegasRELAX[2][1] >  $nMaxPosOmega) { // if background is larger than max
					$arrOmegasRELAX[2][1] = $nMaxPosOmega / 1.5;
				}
				$arrOmegasRELAX[2][0] = $arrOmegasRELAX[2][1] * 1.5;
			}
	} while(!fnCheckOmegaConstraints($arrOmegasRELAX));
	print_r($arrOmegasRELAX);	//die();					
							
	$arrOmegaProportionsRELAX= array($arrRelaxRecord[1],$arrRelaxRecord[3],$arrRelaxRecord[5]);
	//print_r($arrOmegaProportionsRELAX);
	$arrMarkedTree = fnMarkTreeWithOmegas($sTree, $arrOmegasRELAX);
	//print_r($arrMarkedTree);
	
	
	return array('codonfreq' => $arrCodonFreq, 'treelen' => $nTreeLen , 'tree' => $sTree , 'kappa' => $nKappa , 'OmegaTrees' => $arrMarkedTree  ,'OmegaProportions' => $arrOmegaProportionsRELAX);
}

function fnFindGapPattern($sInFasta) {
	$hF = fopen($sInFasta , 'r');
	$arrRet = array();
	$sSeqName = '';
	$sSeq = '';

	$nSeqLen = 0;
	do {
		$sLn = fgets($hF);
		if ($sLn === false || $sLn[0] == '>') {
			if ($sSeqName != '') {
				$arrRet[$sSeqName] = $sSeq;
				$nSeqLen = strlen($sSeq);
			}
			$sSeqName = substr(trim($sLn) , 1);
			$sSeq = "";

			if ($sLn === false) break;
			continue;
		}

		$sSeq .= trim($sLn);
	} while(true);

	$arrGapPatterns = array();
	$arrNPatterns = array();

	foreach($arrRet as $sSp => $sSeq) {
		if (strlen($sSeq) != $nSeqLen ) {
			die("Alignment of different lengths!\n");
		}
		preg_match_all('/N+/', $sSeq , $arrM, PREG_OFFSET_CAPTURE);
		$arrNPatterns[$sSp] = array();
		foreach($arrM[0] as $arrMatch) {
			$arrNPatterns[$sSp][] = array( $arrMatch[1] , strlen($arrMatch[0])); //offset, length
		} 

		preg_match_all('/-+/', $sSeq , $arrM, PREG_OFFSET_CAPTURE);
		$arrGapPatterns[$sSp] = array();
		foreach($arrM[0] as $arrMatch) {
			$arrGapPatterns[$sSp][] = array( $arrMatch[1] , strlen($arrMatch[0])); //offset, length
		} 
	}

	return array('gap' => $arrGapPatterns , 'n' =>$arrNPatterns , "seqlen" => $nSeqLen, "ntax" => count($arrRet));
}

function fnGetSimAln($hMCPAML , $arrNGapPattern) {

	$arrSimAln = array();
	$arrSimAlnGapped = array();
	while( false !== ($sLn = fgets($hMCPAML) ) ) {
		$sLn = trim($sLn);
		$arrF = preg_split("/\s{2,}/" , $sLn);
		if (trim($arrF[0]) == '') continue;
		if (trim($arrF[1]) == '') continue;
		if (strlen($arrF[1]) < 12 ) continue;
		$sTax = trim($arrF[0]);
		$sSeq = str_replace(' ' , '', trim($arrF[1]));
		$arrSimAln[$sTax] = $sSeq;

		if (!array_key_exists($sTax , $arrNGapPattern['n'])) {
			echo("Missing data pattern undefined for $sTax\n");
			return false;
		}

		if (!array_key_exists($sTax , $arrNGapPattern['gap'])) {
			echo("Gap pattern undefined for $sTax\n");
			return false;
		}

		foreach($arrNGapPattern['n'][$sTax] as $arrPattern) {
			$sSeq = substr_replace($sSeq, str_repeat('N', $arrPattern[1]), $arrPattern[0] , $arrPattern[1]);
		}
		foreach($arrNGapPattern['gap'][$sTax] as $arrPattern) {
			$sSeq = substr_replace($sSeq, str_repeat('-', $arrPattern[1]), $arrPattern[0] , $arrPattern[1]);
		}

		$arrSimAlnGapped[$sTax] = $sSeq;

	}

	return array('simaln' => $arrSimAln , 'gappedaln' => $arrSimAlnGapped);

	
}

function fnGetPrevResults($sOutFasta, $sOutGapFasta) {
	if (!file_exists($sOutGapFasta)) {
		return array();
	}
	$hPrevID = popen('grep -oP "GroupId:\S+" '. $sOutGapFasta.' | uniq', 'r');
	$arrDone = array();
	while( false !== ($sLn = fgets($hPrevID) )) {
		$sLn = trim($sLn);
		$arrF = explode(':' , $sLn);
		$arrDone[$arrF[1]] = true;
	}

	return $arrDone;
}

?>
