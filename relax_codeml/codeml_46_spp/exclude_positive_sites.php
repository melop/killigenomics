<?php
//sometimes codeml and relax analyses detect the same genes
//this could mean that both positive and relaxed selection operates on the protein,
//or could be that either method is counfounded by positive selection (hyphy relax) / relaxed selection (codeml)
//this script reads in the list of genes that have conflicting results between codeml and relax
//and eliminate the positively selected sites detected by BEB in codeml
//a new DNA alignment file is produced so that it can be used to rerun codeml and relax.
//if the results in codeml and relax are driven by different sites in the alignment, then relax should still show the same result while codeml will no longer be significant
//if both becomes not significant, it means that either relax or codeml has made a mistake.

$nMinPosterior = 0.85; //mean posterior probablility in BEB test in codeml for deleting the site.
$sSp = "Nothobranchius";
$sConflictList = "sum_$sSp.fdr.relax_conflicts.txt";
$sOutDir = "relax_codeml_conflict/$sSp/bed_pp_cutoff_$nMinPosterior";

//print_r(fnParseBEBExclude('/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/codeml_46_spp/codemltmp/7583511/codeml.3.out')); 

//die();
exec("mkdir -p $sOutDir");

$hList = fopen($sConflictList , "r");
$hFasta = fopen($sOutDir."/ret_aln_possite_excluded_DNA.fasta" , 'w');
$nCount = 0;
while( false !== ($sLn = fgets($hList) )) {
	$sLn = trim($sLn);
	if ($sLn == '') continue;
	$arrF = explode("\t" , $sLn);
	$sOrthoID = $arrF[0];
	$sCodemlPath = $arrF[7];
	$sCodemlOut = "$sCodemlPath/codeml.3.out";
	$sInFasta = "$sCodemlPath/in.fas";
	if (!file_exists($sCodemlOut )) {
		echo("$sCodemlOut does not exist any more, cannot do\n");
		continue;
	}

	$arrExcludeSites = fnParseBEBExclude($sCodemlOut);
	$arrFasta = fnParseFasta($sInFasta, $arrExcludeSites);

	foreach($arrFasta as $sSq => $sSeq) {
		$sHeader = "Ortho:$nCount;Sp:$sSq;MappedToRef:$sSq;SpGeneId:000000;GroupId:$sOrthoID";
		fwrite($hFasta , ">$sHeader\n".wordwrap($sSeq, 75, "\n" , true)."\n");
	}

	$nCount++;
}

function fnParseBEBExclude($sCodemlOut) {
	global $nMinPosterior;
	$hIn = popen(" grep -A10000000 \"Bayes Empirical Bayes (BEB) analysis\" \"$sCodemlOut\" | grep -B10000000 \"The grid\" | head -n -3 " , 'r');
	$arrRet = array();; 
	while( false !== ($sLn = fgets($hIn) )) {
		if ($sLn[0]!=' ') continue;
		$sLn = trim($sLn);
		list($nPos, $sAA, $nProb) = preg_split("/\s+/" , $sLn);
		$nProb = floatval($nProb);
		if ($nProb >= $nMinPosterior) {
			$arrRet[$nPos] = array($sAA, $nProb);
		}
	}

	pclose($hIn);

	return $arrRet;
}

function fnParseFasta($sInFasta, $arrExcludeSites) {
	$hIn = fopen($sInFasta , 'r');

        $arrRet = array();
        $sSeqName = "";
        $sSeq = "";
        do {
                $sLn = fgets($hIn);
                if ($sLn === false || $sLn[0] == '>') {
                        if ($sSeqName != '') {
                                $arrRet[$sSeqName] = $sSeq;
                        }

                        if ($sLn === false) break;
                        $sSeqName = substr(trim($sLn), 1);
                        $sSeq = "";
                        continue;
                }

                $sSeq .= trim($sLn);
        } while(true);

	$arrExcRet = array();
        foreach($arrRet as $sSp => $sSeq) {
		$nStartCut = 0;
		$sTrimSeq = "";
		foreach($arrExcludeSites as $nPos => $arrInfo) {
			$nEndCut = $nPos*3 - 3;
			$sTrimSeq .= substr($sSeq , $nStartCut , $nEndCut - $nStartCut);
			$nStartCut = $nPos*3;
		}
		if ($nStartCut < strlen($sSeq)) {
			$sTrimSeq .= substr($sSeq , $nStartCut);
		}
		$arrExcRet[$sSp] = $sTrimSeq;
	}
	
	return $arrExcRet;
}
?>
