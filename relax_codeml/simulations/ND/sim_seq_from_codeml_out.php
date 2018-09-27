<?php
//simulate sequence alignments from codeml output (branch-site tests)
//the temporary folder must be saved for this to work. in particular, codeml.3.out and in.fas need to be present in the temp folder.
//the missing data pattern is directly copied from in.fas (with N's and gaps copied to the simulated dataset).

$arrCodemlRet = glob("../..//hyphyrelax/codeml_46_spp/Codeml_Nothobranchius/ret_part*.of.200.txt");
$sTmp = "tmp/".rand();
$nRandSeed  = 13888;
$sEvolver = "/software/source/paml4.8/src/evolverNSsites";
$sUseOmegaType = "average"; // can be "average", "foreground", "background";
$sOutFasta = "ret_sim.fasta";
$sOutGapFasta = "ret_sim_gapped.fasta";

$arrAlreadyDone = fnGetPrevResults($sOutFasta, $sOutGapFasta);

$hOutFasta = fopen($sOutFasta, 'a');
$hOutGapFasta = fopen($sOutGapFasta, 'a');

exec("mkdir -p $sTmp");
//fnSim('Group_10_6' , '/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/codeml_46_spp/codemltmp/9227585/in.fas' , '/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/codeml_46_spp/codemltmp/9227585/codeml.3.out');
//die();

foreach($arrCodemlRet as $sCodemlRet) {
	$hCodemlRet = fopen($sCodemlRet , 'r');
	$nGroupCount = 0;
	while( false !== ($sLn = fgets($hCodemlRet) )) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;

		$nGroupCount++;

		$arrF = explode("\t", $sLn);
		$sGroup = $arrF[0];

		if (array_key_exists($sGroup , $arrAlreadyDone)) {
			echo("$sGroup skipped, already done.\n");
			continue;
		}

		if ( trim($arrF[5]) == '' || trim($arrF[6]) == '') {
			echo("$sGroup skipped because codeml did not finish.\n");
			continue;
		}

		$sCodemlPath =  $arrF[7];
		$sInFasta = "$sCodemlPath/in.fas";
		$sCodeml3Out = "$sCodemlPath/codeml.3.out";
		if (!file_exists($sInFasta) ) {
			die("$sInFasta not found\n");
		}

		if (!file_exists($sCodeml3Out) ) {
			die("$sCodeml3Out not found\n");
		}

		$arrSimAln = fnSim($sGroup, $sInFasta, $sCodeml3Out);

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
}

function fnSim($sGroup, $sInFasta, $sCodeml3Out) {
	global $sTmp, $nRandSeed , $sUseOmegaType, $sEvolver;
	exec("rm $sTmp/*");
	$arrNGapPattern = fnFindGapPattern($sInFasta);
	//print_r($arrNGapPattern);
	echo($sInFasta."\n");
	$sMCCodon = "$sTmp/MCcodonNSsites.dat";
	$hMCCodon = fopen($sMCCodon , 'w');
	fwrite($hMCCodon , "0          * 0,1:seqs or patters in paml format (mc.paml); 2:paup format (mc.nex)\n");
	fwrite($hMCCodon , "$nRandSeed       * random number seed (odd number)\n");
	fwrite($hMCCodon , $arrNGapPattern['ntax']." ".intval($arrNGapPattern['seqlen']/3)." 1   * <# seqs>  <# codons>  <# replicates>\n\n");
	
	$arrSimInfo = fnCollectSimInfo($sCodeml3Out);

	//print_r($arrSimInfo);
	fwrite($hMCCodon , $arrSimInfo['treelen']."           * <tree length; see note below; use -1 if tree has absolute branch lengths>\n\n");
	fwrite($hMCCodon , $arrSimInfo['tree']."\n\n");
	fwrite($hMCCodon , count($arrSimInfo['omegatable'][0])."            * number of site classes, followed by frequencies and omega's.\n");

	$nOmegaRow = ($sUseOmegaType == "average")? 3: (($sUseOmegaType == "foreground")? 2:1 );

	foreach(array(0, $nOmegaRow) as $nRow ) {

		foreach($arrSimInfo['omegatable'][$nRow] as $nCol => $nVal) {
			fwrite($hMCCodon , "  " . number_format($nVal, 6));
		}
		fwrite($hMCCodon , "\n");
	}

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

	//now parse the output file
	$hMCPAML = fopen("$sTmp/mc.paml" , 'r');

	$arrSimAlns = fnGetSimAln($hMCPAML , $arrNGapPattern);

	return $arrSimAlns;
}

function fnCollectSimInfo($sCodeml3Out) {
	//find codon frequencies from output:

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

	//read dN/dS categories:
	$h1 = popen("grep -A5 \"dN/dS (w) for site classes\" $sCodeml3Out" , 'r');
	$nLn = 0;
	$arrOmegaTable = array();
	while( false !== ($sLn = fgets($h1) ) ) {
		$nLn++;
		if ($nLn >= 4) {
			preg_match_all("/\d+\.?\d+/", $sLn, $arrF);
			if (count($arrF[0]) !=4) {
				echo($sLn);
				die("Cannot parse omega table\n");
			}

			$arrOmegaTable[] = $arrF[0];
		}

	}

	pclose($h1);

	$arrOmegaTable[3] = $arrOmegaTable[2];
	//compute average omega between background and foreground:

	foreach($arrOmegaTable[3] as $nCat => $nVal) {
		$arrOmegaTable[3][$nCat] = ($arrOmegaTable[1][$nCat] + $arrOmegaTable[2][$nCat]) / 2;
	}
	
	return array('codonfreq' => $arrCodonFreq, 'treelen' => $nTreeLen , 'tree' => $sTree , 'kappa' => $nKappa , 'omegatable' => $arrOmegaTable  );
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
			print_r($arrNGapPattern);
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
