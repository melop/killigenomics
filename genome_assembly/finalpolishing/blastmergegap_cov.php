<?php
/*
Blast a scaffold against itself, then look for tandem duplications of very high identity that probably results from overlapping overhangs
This version keeps the junction site:
 <==============Dup1===========>NNNNNN==JUNCTION==<==============Dup2===========>
results in:
 NNNNNN==JUNCTION==<==============Dup2===========>

 <==============Dup1===========>==JUNCTION==NNNNNN<==============Dup2===========>
results in:
 <==============Dup1===========>==JUNCTION==NNNNNN

*/

$sIn = "./release/NOR/1.0/NORv1.0.fa"; #assembly to be fixed
$sOut = "blast.gapmerged.fa"; //output file
$sTmpFolder = "./tmp/";
$sAGP = "overlap.joined.agp"; //agp output
$nLnWidth = 75; //line width in output

$nTotalPart = 1;
$nThisPart = 0;
$nMinBlockIdentity = 90; //min sequence identity between the overhangs
$nMaxHitDistance = 100000; //if hits are more than 100kb away,ignore them, probably just interspersed repeats.
$nMinOverallIdentity = 0.95;
$nMaxGapPerc = 0.15;
$nMaxInterHitDistance = 3000; //if the sub hits are < this bp away, merge and treat as one big block
$nMaxNonNDistanceBetweenDups = 5000; // the two duplicated blocks need to be <5000bp
$nMaxNonNDistanceBetweenDupsPerc = 1; // the two duplicated blocks need to be <1% of duplicated block length

$sBAMFile = "./mapped/mergedPE.bam";
$nAvgGenomeCov = 64.96965;
$sRepeatGFF = "./annotation/repeatmasker/NOR_v1.0/scf.fa.out.gff"; //repeat annotations
$nMinLenAfterRepeatExclusion = 20; //at least 20bp left after masking the repeats in the overlap
$nMinHaploidLODCutoff = log10(5); //0.7 is about 5 times more likely . log10 base, if this is 1, then it means that it's 10 times more likely the region is haploid than diploid
$nMinKeepSegmentCovDiffLOD = log10(1.3); //if abs(log10(coverage left / coverage right))) is less than this cutoff, then keep the part where it joins continuously with the junction. otherwise keep the part with higher coverage.

//echo(fnGetCov("NORscf444" , 19739 , 21685, true, false, true));
//die();

while(count($argv) > 0) {
    $arg = array_shift($argv);
    switch($arg) {
        case '-N':
            $nTotalPart  = trim(array_shift($argv));
            break;
        case '-f':
            $nThisPart  = trim(array_shift($argv));
            break;
        case '-i':
            $sIn  = trim(array_shift($argv));
            break;
        case '-T':
            $sTmpFolder  = trim(array_shift($argv));
            break;
        case '-b':
            $nMinBlockIdentity  = trim(array_shift($argv));
            break;
        case '-h':
            $nMaxHitDistance  = trim(array_shift($argv));
            break;
        case '-O':
            $nMinOverallIdentity  = trim(array_shift($argv));
            break;
 	case '-M':
            $nMaxGapPerc  = trim(array_shift($argv));
            break;
	case '-D':
            $nMaxInterHitDistance  = trim(array_shift($argv));
            break;
	case '-g':
            $nMaxNonNDistanceBetweenDups  = trim(array_shift($argv));
            break;
	case '-G':
            $nMaxNonNDistanceBetweenDupsPerc  = trim(array_shift($argv));
            break;
	case '-B':
            $sBAMFile  = trim(array_shift($argv));
            break;
	case '-C':
            $nAvgGenomeCov  = trim(array_shift($argv));
            break;
	case '-R':
            $sRepeatGFF  = trim(array_shift($argv));
            break;
	case '-r':
            $nMinLenAfterRepeatExclusion  = trim(array_shift($argv));
            break;

    }
}

$hFastaIn = fopen($sIn, "r");

/*
start communication with poisson engine in R

*/

$oDesc = array(
   0 => array("pipe", "r"),  // stdin is a pipe that the child will read from
   1 => array("pipe", "w"),  // stdout is a pipe that the child will write to
   2 => array("file", "/dev/null", "a") // stderr is a file to write to
);

$arrPoissonPipes = array();
$hPoisson = proc_open('Rscript poisson.R', $oDesc, $arrPoissonPipes);
$nAvgHaploidGenomeCov = $nAvgGenomeCov / 2;


$sSeqName = "";
$sSeq = "";
$nSeqCount = -1;
do {
	$sLn = fgets($hFastaIn);
	if ($sLn === false) {
		fnRun($sSeqName , $sSeq);
		break;
	}

	$sLn = trim($sLn);

	if ($sLn == "") continue;

	if ($sLn[0] == ">") {
		fnRun($sSeqName , $sSeq);
		$sSeqName = substr($sLn, 1);
		$sSeq = "";
		continue;
	}

	$sSeq .= $sLn;		
} while(true) ;

function fnRun($sSeqName , $sSeq) {
	global $sAGP, $sOut, $nTotalPart, $nThisPart, $nSeqCount, $sTmpFolder , $nMinBlockIdentity , $nMinOverallIdentity, $nMaxInterHitDistance, $nMaxNonNDistanceBetweenDups, $nMaxNonNDistanceBetweenDupsPerc, $nMaxGapPerc, $nLnWidth, $nMinHaploidLODCutoff, $nMinKeepSegmentCovDiffLOD;
	if ($sSeqName == "" || $sSeq == "") return;


	$nSeqCount++;

	if ( ($nSeqCount % $nTotalPart) != $nThisPart ) return;

	//echo("$nSeqCount / $nTotalPart  / $nThisPart / " . ($nSeqCount % $nTotalPart)."\n");

	echo("part $nThisPart of $nTotalPart\n");

	$sWD = "$sTmpFolder/$sSeqName/";

	fnCMD("mkdir -p $sWD");

	$sFasta = "$sWD/scf.fa";
	$hFasta = fopen($sFasta, "w");
	fwrite($hFasta, ">$sSeqName\n".$sSeq."\n");

	//now makeblastdb
	$sBlastOut = "$sWD/blast.txt";
	$sBlastFlag = "$sWD/blast.ok";

	if (!file_exists($sBlastFlag)) {
		fnCMD("makeblastdb -in $sFasta  -dbtype nucl");
		fnCMD("blastn -task megablast -query $sFasta -db $sFasta -outfmt 7 -evalue 1e-300 > $sBlastOut && touch $sBlastFlag ");
	}

	if (!file_exists($sBlastFlag)) {
		echo("blast failed: $sSeqName\n");
		return;
	}

	$arrAlnMatches = fnGetAlnMatches($sBlastOut , $nMinBlockIdentity);

	//print_r($arrAlnMatches);

	$bCandidateDupStarted = false;
	$nCandidateLeftBlockStart = 0;
	$nCandidateLeftBlockEnd = 0;

	$nCandidateRightBlockStart = 0;
	$nCandidateRightBlockEnd = 0;

	$nAlnLen = 0;
	$nIdenticalLen = 0;
	$nTotalGapLen = 0;

	$arrDupBlocks = array();
	$arrAlnMatchKeys = array_keys($arrAlnMatches);
	for( $i=0; $i<count($arrAlnMatchKeys); $i++) {
		$oMatch = $arrAlnMatches[$arrAlnMatchKeys[$i]];
		$nQueryStart = $oMatch[6];
		$nQueryEnd = $oMatch[7];

		$nHitStart = $oMatch[8];
		$nHitEnd = $oMatch[9];

		$nDivergence = $oMatch[2]/100;
		$nAlnBlockLen = $oMatch[3];

	

		if ($nQueryStart >= $nHitStart) {
			if (!$bCandidateDupStart) { // now the candidate duplicate junction point is found

				$nCandidateLeftBlockStart = $nHitStart;
				$nCandidateLeftBlockEnd = $nHitEnd;
				$nCandidateRightBlockStart = $nQueryStart;
				$nCandidateRightBlockEnd = $nQueryEnd;

				$nAlnLen = $nAlnBlockLen;
				$nIdenticalLen = $nAlnBlockLen * $nDivergence;
				$nTotalGapLen = 0;
					echo("init alnlen: $nAlnLen iden: $nIdenticalLen  totalgap: $nTotalGapLen\n");
				$bCandidateDupStart = true;
			} else { // if the match block has started already
				if ( abs($nQueryStart - $nCandidateRightBlockEnd )+ 1 <= $nMaxInterHitDistance && abs($nHitStart - $nCandidateLeftBlockEnd )+ 1 <= $nMaxInterHitDistance ) { // extend the block
					//sanity check:
					if ( $nHitEnd > $nCandidateRightBlockStart ) { // this hit ends too far down stream, this is complex region, don't attemp to merge
						$bCandidateDupStart = false; //save old block
						$nAlnLen = 0;
						continue;
					}

					$nAlnLen += $nAlnBlockLen;
					$nIdenticalLen += $nAlnBlockLen * $nDivergence;
					$nTotalGapLen += ($nQueryStart - $nCandidateRightBlockEnd >=0)? ($nQueryStart - $nCandidateRightBlockEnd + 1):0;
					$nCandidateLeftBlockEnd = ($nCandidateLeftBlockEnd < $nHitEnd)? $nHitEnd : $nCandidateLeftBlockEnd;
					$nCandidateRightBlockEnd = ($nCandidateRightBlockEnd < $nQueryEnd)? $nQueryEnd : $nCandidateRightBlockEnd ;

					echo("alnlen: $nAlnLen iden: $nIdenticalLen  totalgap: $nTotalGapLen\n");
				} //otherwise, do not add further matches to the block
			}
		} 

		if($nQueryStart < $nHitStart || $i == (count($arrAlnMatches)-1) )
		{
			$bCandidateDupStart = false; //save old block

			if ($nAlnLen >0 ) {
				//compute block stats
				$nFinalIdentity = $nIdenticalLen / $nAlnLen;
				$nPercentGapBetweenHits = $nTotalGapLen / ($nAlnLen + $nTotalGapLen);
				$nJunctionGapLen = $nCandidateRightBlockStart - $nCandidateLeftBlockEnd - 1;
				$nJunctionGapLen = ($nJunctionGapLen<=0)? 0:$nJunctionGapLen;
				$nJunctionNonNPerc = 0;
				$sJunctionGap = "";
				if ($nJunctionGapLen > 0) {
					$sJunctionGap = substr($sSeq , $nCandidateLeftBlockEnd , $nJunctionGapLen);
					$nJunctionNonNCount = strlen(str_replace( 'N', '', $sJunctionGap)) ;
					$nJunctionNonNPerc = $nJunctionNonNCount / ($nAlnLen + $nTotalGapLen);
				}

				$nLeftBlockLen = $nCandidateLeftBlockEnd -$nCandidateLeftBlockStart+1;
				$nRightBlockLen = $nCandidateRightBlockEnd -$nCandidateRightBlockStart+1;

				if ($nLeftBlockLen <= 0 ) die("ERROR: nLeftblocklen <= 0!\n"); 
				if ($nRightBlockLen <= 0 ) die("ERROR: nRightBlockLen <= 0!\n"); 
				echo("alnlen: $nAlnLen ; gapped alnlen ". ($nAlnLen + $nTotalGapLen) . " $nCandidateLeftBlockStart - $nCandidateLeftBlockEnd vs $nCandidateRightBlockStart - $nCandidateRightBlockEnd : $nFinalIdentity , gap $nPercentGapBetweenHits, junction gap: $nJunctionNonNCount ($nJunctionNonNPerc) \n");
					echo("left block: \n".substr($sSeq , $nCandidateLeftBlockStart -1, $nLeftBlockLen ) . "\n");
					echo("junction:\n$sJunctionGap\n");
					echo("right block: \n".substr($sSeq , $nCandidateRightBlockStart -1, $nRightBlockLen ) . "\n");


				$bAcceptMerge = ($nFinalIdentity >= $nMinOverallIdentity && $nPercentGapBetweenHits <= $nMaxGapPerc && $nJunctionNonNCount <= $nMaxNonNDistanceBetweenDups && $nJunctionNonNPerc <= $nMaxNonNDistanceBetweenDupsPerc) ;

				if ($bAcceptMerge) { // now, check the read coverage in these two duplicated regions. They should both be haploid depth
					$nLeftBlockCov = fnGetCov($sSeqName , $nCandidateLeftBlockStart  , $nCandidateLeftBlockEnd, true, false, true);
					if ($nLeftBlockCov === false) {
						$nLeftBlockCov = fnGetCov($sSeqName , $nCandidateLeftBlockStart  , $nCandidateLeftBlockEnd, true, false, false); //try it without masking
					}
					$nRightBlockCov = fnGetCov($sSeqName , $nCandidateRightBlockStart  , $nCandidateRightBlockEnd, true, false, true);
					if ($nRightBlockCov === false) {
						$nRightBlockCov = fnGetCov($sSeqName , $nCandidateRightBlockStart  , $nCandidateRightBlockEnd, true, false, false); //try it without masking
					}

					$nLeftHaploidLOD = fnHaploidLOD($nLeftBlockCov);
					$nRightHaploidLOD = fnHaploidLOD($nRightBlockCov);
					$nAvgHaploidLOD = fnHaploidLOD( ($nLeftBlockCov+$nRightBlockCov) / 2 );
					echo("left cov: $nLeftBlockCov, left haploid LOD: $nLeftHaploidLOD ; right cov: $nRightBlockCov, right haploid LOD: $nRightHaploidLOD ; Average Haploid LOD: $nAvgHaploidLOD\n");

					$bAcceptMerge = ( $bAcceptMerge &&  $nAvgHaploidLOD >= $nMinHaploidLODCutoff ); // must be likely both haploids
				}
				
				if ($bAcceptMerge) {
					echo("accepted\n"); 
					//figure out keep left or keep right, by checking the N patterns in the junction sequence, keep the side with shorter Ns.
					$bKeepLeft = ($nLeftBlockCov > $nRightBlockCov); //by default, saves the segment with better mapping rate
					$nDiffCov = abs(log10($nLeftBlockCov / $nRightBlockCov));
					if (strlen($sJunctionGap)>2 && $nDiffCov<$nMinKeepSegmentCovDiffLOD) { // use the part that doesn't join with the junction with N's
						$sJunctionGap = strtoupper($sJunctionGap);
						preg_match_all('/^N+|N+$/', $sJunctionGap , $arrMatches, PREG_OFFSET_CAPTURE);
						if (count($arrMatches[0]) == 0 || $nJunctionNonNCount == 0) {
							//no N's, just keep the longer segment, plus the junction
						} else {
							//find the longest streth of N
							$nLongestIndex = 0;
							$nLongestNstretch = 0;
							$nLongestNStart = 0;
							for($nNIdx=0; $nNIdx<count($arrMatches[0]); $nNIdx++) {
								$nNLen = strlen($arrMatches[0][$nNIdx][0]);
								$nNStart = $arrMatches[0][$nNIdx][1];
								if ($nNLen > $nLongestNstretch) {
									$nLongestNstretch = $nNLen;
									$nLongestIndex = $nNIdx;
									$nLongestNStart = $nNStart;
								}	
							}

							if ($nLongestNStart < (strlen($sJunctionGap)/2) ) { //if a long streth of N's is on the left side , keep the right side
								$bKeepLeft = false;
							} else {

								$bKeepLeft = true;
							}
						}
					}

					if ($bKeepLeft) {echo("keep left\n");} else {echo("keep right\n");}
					$arrDupBlocks[$nCandidateLeftBlockStart] = array($nCandidateLeftBlockStart ,$nCandidateLeftBlockEnd, $nCandidateRightBlockStart, $nCandidateRightBlockEnd , $bKeepLeft);
				} else {
					echo("rejected\n");
				}
			}

			$nAlnLen = 0;
			continue;

		}
	}



	$hOut = fopen("$sWD/$sOut" , "w");
	$hAGP = fopen("$sWD/$sAGP" , "w"); 

	ksort($arrDupBlocks);

	//from dupblocks to alignment blocks (for agp)
	$arrBlocks = array();
	$nPrevBlockEnd = 0;

	foreach($arrDupBlocks as &$oDupBlock) {
		if ($oDupBlock[0] < $nPrevBlockEnd + 1) {
			continue; //the this block is nested within the previous block, don't need to add.
		}
		if ($oDupBlock[4]) { //keep left
			$arrBlocks[] = array( $nPrevBlockEnd + 1, $oDupBlock[2] -1 ); //keep the left sequence + junction gap
		} else {
			if ( ($oDupBlock[0] -1) - ($nPrevBlockEnd + 1) > 0 ) {
				$arrBlocks[] = array( $nPrevBlockEnd + 1, $oDupBlock[0] -1 ); //keep the junction gap + right sequence
			}
			$arrBlocks[] = array( $oDupBlock[1] + 1, $oDupBlock[3] ); //keep the junction gap + right sequence


		}
		$nPrevBlockEnd = $oDupBlock[3];
	}

	if ($nPrevBlockEnd != strlen($sSeq) ) {
		$arrBlocks[] = array($nPrevBlockEnd + 1, strlen($sSeq)); 
	}

	//print_r($arrDupBlocks);


	$sNewSeq = "";
	$nNewBlockStart = 1; //this is one based
	$nComponentNum = 1;
	foreach($arrBlocks as &$oBlock) {
/*
EG1_scaffold1	1	3043	1	W	AADB02037551.1	1	3043	+
EG1_scaffold2	1	40448	1	W	AADB02037552.1	1	40448	+

*/
		$nBlockLen = $oBlock[1] - $oBlock[0] + 1;
		if ($nBlockLen <= 0) { print_r($arrDupBlocks ); print_r($arrBlocks); die("blocklen <= 0! $nBlockLen\n"); }
		$sNewSeq .= substr($sSeq, $oBlock[0] -1 , $oBlock[1] - $oBlock[0] + 1);  
		fwrite($hAGP, $sSeqName."\t".$nNewBlockStart."\t".strlen($sNewSeq)."\t".$nComponentNum."\tW\t$sSeqName\t".($oBlock[0])."\t".($oBlock[1])."\t+\n" );
		$nNewBlockStart = strlen($sNewSeq) + 1;
		$nComponentNum++;
	}

	fwrite($hOut, ">$sSeqName\n".wordwrap($sNewSeq, $nLnWidth, "\n", true )."\n");

}

function fnGetAlnMatches($sBlastOut , $nMinBlockIdentity) {
	global $nMaxHitDistance;
	$hB = fopen($sBlastOut , 'r');

	$arrAln = array(); //index is query start position
	while(false !== ($sLn = fgets($hB) )) {
		$sLn = trim($sLn);
		if ($sLn == "") continue;
		if ($sLn[0] == '#') continue;
		$arrF = explode("\t" , $sLn);

		if ( $arrF[6] == $arrF[8] && $arrF[7] == $arrF[9] ) continue;

		if ( $arrF[6] > $arrF[7] || $arrF[8] > $arrF[9] ) continue; //don't look at reverse matches

		if ($arrF[2] < $nMinBlockIdentity ) continue;

		if (abs($arrF[6] - $arrF[8]) > $nMaxHitDistance) continue;

		$arrAln[$arrF[6]] = $arrF;
	}

	ksort($arrAln , SORT_NUMERIC);
	return $arrAln;
}

function fnCMD($sCmd, $bEcho=true) {
	if ($bEcho) {echo($sCmd.PHP_EOL);}
	exec($sCmd);
}

function fnGetCov($sScf , $nStart, $nEnd, $bKeepTmpBam = true, $bOnlyMakeBam = false, $bIgnoreRepeat = true) { //$start end 1 based
	global $sBAMFile , $nAvgGenomeCov , $sRepeatGFF, $sTmpFolder, $nMinLenAfterRepeatExclusion;

	$sWD = "$sTmpFolder/$sScf/";
	fnCMD("mkdir -p $sWD", false);


	$sBed = "$sScf\t".($nStart-1)."\t$nEnd\n"; 
	$sRegionBed = "$sWD/tmp.bed";
	$hRegionBed = fopen($sRegionBed , "w");
	fwrite($hRegionBed , $sBed);

	if ( (!$bKeepTmpBam) || ( !file_exists("$sWD/tmp.sub.bam")) ) {
		fnCMD("samtools view -b $sBAMFile $sScf > $sWD/tmp.sub.bam" , false);
	}


	if ($bOnlyMakeBam) return;

	if ($bIgnoreRepeat) {
		if (!file_exists("$sWD/scfrepeats.gff") ) {
			fnCMD("grep -P \"$sScf\t\" $sRepeatGFF  > $sWD/scfrepeats.gff " , false);
		}

		fnCMD("bedtools subtract -a $sRegionBed  -b $sWD/scfrepeats.gff > $sRegionBed.norepeat", false);

		//check the output
		$hBedNoRepeats = fopen("$sRegionBed.norepeat" , "r");
		$nAccLen = 0;
		while( false !== ($sLn = fgets($hBedNoRepeats))) {
			$sLn = trim($sLn);
			$arrF = explode("\t", $sLn);
			if (count($arrF) <3) continue;
			$nAccLen += $arrF[2] - $arrF[1] ;
		}

		if ($nAccLen < $nMinLenAfterRepeatExclusion) {
			return false;
		}

		$sRegionBed = "$sRegionBed.norepeat";
	}


	//get the coverage count:
	fnCMD("bedtools coverage -sorted -a $sRegionBed -b $sWD/tmp.sub.bam -d > $sWD/tmp.cov");
	return fnGetAvgFromCovFile("$sWD/tmp.cov");

}

function fnGetAvgFromCovFile($sFile) {
	$hCov = fopen($sFile , "r");
	$nSumCov = 0;
	$nSites = 0;
	while(false !== ($sLn = fgets($hCov) )) {
		$sLn = trim($sLn);
		$arrF = explode("\t" , $sLn);
		if (count($arrF) != 5) continue;
		$nSumCov += $arrF[4];
		$nSites++;
	}

	return $nSumCov / $nSites;
}

function fnHaploidLOD($nCov) {
	global $nAvgHaploidGenomeCov, $nAvgGenomeCov, $arrPoissonPipes, $hPoisson;

	if (!$hPoisson) {
		die("Error: Communication between R is broken.\n");
	}

	fwrite($arrPoissonPipes[0] , "$nCov\t$nAvgGenomeCov\n");
	$nDiploidProb = floatval(fgets($arrPoissonPipes[1]) );

	fwrite($arrPoissonPipes[0] , "$nCov\t$nAvgHaploidGenomeCov\n");
	$nHaploidProb = floatval(fgets($arrPoissonPipes[1]) );

	return log10($nHaploidProb) - log10($nDiploidProb);
	
}

?>
