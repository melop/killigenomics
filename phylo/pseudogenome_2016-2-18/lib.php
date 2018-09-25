<?php
ini_set('memory_limit','300000M');
require_once(dirname(__FILE__) . "/config.php");
require_once(dirname(__FILE__) . "/codons.php");

ini_set('display_errors','On');
error_reporting(E_ALL|E_STRICT);

	$Mask_Operation = MASK_BLASTX; // do blastn search
	$Mask_Blast_Db = MASK_DB_NR; // Usse nt database
	$flag_mask = $Mask_Operation | $Mask_Blast_Db ;

	

	
	function fnWriteText($file, $content) {
		$hFile = fopen($file.Ext_Writing,"w");
		
		if (!$hFile) {
			return false;
		}
		
		fwrite($hFile, $content);
		
		fclose($hFile);
		
		rename( $file.Ext_Writing, $file);
		
	}
	
	function fnReadText($file) {
	/*
		$fh = fopen($file, 'r');
		if (!$fh) {
			return false;
		}
		
		$theData = fgets($fh);
		fclose($fh);
		return $theData;
		*/
		
		return file_get_contents($file);
	}
	
	function fnLog($content) {
		global $tmp_path;
		
		$log_file = $tmp_path."conductor.log";
		$hFile = fopen($log_file,"a");
		
		if (!$hFile) {
			return false;
		}
		
		fwrite($hFile, $content.PHP_EOL);
		
		fclose($hFile);
		
	}
	
	
	function microtime_float()
	{
		list($usec, $sec) = explode(" ", microtime());
		return ((float)$usec + (float)$sec);
	}

	
class HspUtil
{
	private $arrHsps = array();
	
	public function __construct() { }
	
	public function AddHsp($nLeft, $nRight, $nFrame) {
	
		$nCountArrHsps = count($this->arrHsps); 
		for ($i=0;$i<=$nCountArrHsps-1;$i++) //Compare the new hsp to the current summarized hsps
		{
			$nHspLeft = $this->arrHsps[$i][0];
			$nHspRight = $this->arrHsps[$i][1];
			$nHspFrame = $this->arrHsps[$i][2];
			
			if ($nHspFrame != $nFrame) {
				continue; //frames don't fit, don't need to compare
			}
			
			//Now see if overlaps
			if ( $nRight < $nHspLeft || $nLeft > $nHspRight) //no overlap
			{
				continue;
			}
			
			//They have overlap, update the current Hsp definition instead of appending it
			
			$this->arrHsps[$i][0] = ( $nHspLeft < $nLeft )? $nHspLeft:$nLeft;
			$this->arrHsps[$i][1] = ( $nHspRight > $nRight)? $nHspRight:$nRight;
			$this->arrHsps[$i][3] = ($this->arrHsps[$i][0] + $this->arrHsps[$i][1]) / 2; //average pos
			return;
			
		}
		
		//If this new hsp was not found in the current list, append it
		//echo("appendhsp");
		$this->appendHsp($nLeft, $nRight, $nFrame);
		return;
		
	}
	
	public function GetSummarizedHsps() 
	{
		
		if (count($this->arrHsps) == 0) { //no hsps
			return false;
		}
		
		$this->orderHsps();
		return  $this->arrHsps;
		
		
	}
	
	private function appendHsp($nLeft, $nRight, $nFrame) {
	
		$hsp = array( 0 => $nLeft, 1 => $nRight, 2=> $nFrame, 3 =>  ($nLeft+$nRight)/2 );
		
		array_push($this->arrHsps, $hsp); // push it at the end
		//echo(count($this->arrHsps));
		
	}
	/*
	private function orderHsps() 
	{
		uasort($this->arrHsps, 
		
		function ($Hsp1, $Hsp2) 
			{
				//echo("ordercallback".PHP_EOL);
				if ($Hsp1[2] > 0) { //coding chain on current order
					return $Hsp1[3] > $Hsp2[3]? 1:-1 ;
				}
				else { // coding chain on reverse complement
					return $Hsp1[3] < $Hsp2[3]? 1:-1;
				}
			}
		
		);
		//echo(count($this->arrHsps));
		return;
		
	}*/
	

	
}


class AlignmentUtil {
	
	private $arrLabels = array(); // Labels
	private $arrRawSeq = array(); // Raw sequences
	private $arrTransSeq = array(); // Translated peptide sequence
	private $arrAlnSeq = array(); // Aligned Sequences
	private $arrTransToAlnMap = array(); //mapping of translated sequence to alned sequences
	private $arrAlnToTransMap = array(); //mapping of alned sequence to translated sequences
	private $sConsensusLn = ""; //The consensus line containing *, : and .
	private $arrExcludeMask = array();
	private $arrStrongAreas = array();
	
	private $ScoreMatrix = array("*" => 10, ":"=>4, "." => 1, " " => -2); // highest score is 10!! 10 8 5 0
	private $SumScoreMatrix = array("*" => 10, ":"=>2, "." => 1, " " => -2); // highest score is 10!! 10 8 5 0
	private $nLookLeft = 5; // look at 3 positions to the left of the current position
	private $nLookRight = 5; //look at 3 positions to the right of the current position
	private $nScoreCutoff = 0.5; //cutoff in percentage

	/*
	private $sTmpDir = CLUSTAL_TMP;
	private $sClustalExe = CLUSTAL_PATH;
	*/
	
	private $nInstanceID ;
	private $sTmpFileBase;
	private $bOnlyFilterWeakFlanks=false; // whether filter weak positions on flankings only or scan the whole sequence 
	
	const SEQ_TYPE_DNA=0;
	const SEQ_TYPE_CODING_DNA=2;
	const SEQ_TYPE_PEPTIDE=1;
	
	private $rawSequenceType = AlignmentUtil::SEQ_TYPE_DNA; 
	
	public function __construct() {
      $this->nInstanceID = rand( 1000000 , 9999999) ;
	  $this->sTmpFileBase = CLUSTAL_TMP."/CLUSTALTMP_".$this->nInstanceID."/";
	}
	
	public function TrimWeakFlanksOnly($b) {
		
		if (!isset($b)) {
			return $this->bOnlyFilterWeakFlanks;
		}
		
		$this->bOnlyFilterWeakFlanks = 	$b;
		
	}
	
	public function AddSequence($Label, $s, $type, $bTranslate) {
	
		$s = fnCleanDNASeq($s);
		$s = preg_replace('/\|/', "", $s); //Replace all HSP dividers
		
		array_push( $this->arrLabels , $Label ); 
		
		array_push( $this->arrAlnSeq , "");
		array_push( $this->arrTransToAlnMap, array() ); //initialize map 
		array_push( $this->arrAlnToTransMap, array() );
		array_push( $this->arrStrongAreas, "" );
		
		if ( ($type == AlignmentUtil::SEQ_TYPE_DNA && $bTranslate ) || $type == AlignmentUtil::SEQ_TYPE_CODING_DNA ) {
		
			$this->rawSequenceType = AlignmentUtil::SEQ_TYPE_CODING_DNA; 
		
			$arrTranslated = fnDNA2Peptide($s, true); //throw away any stop codons
			array_push($this->arrTransSeq, $arrTranslated[0] );
			array_push($this->arrRawSeq, $arrTranslated[1] );
			/*
			echo(strlen($this->arrTransSeq[0] ). "\n");
			echo(strlen($this->arrRawSeq[0]));
			*/
			return;
		}
		elseif ($type == AlignmentUtil::SEQ_TYPE_PEPTIDE) {
			$this->rawSequenceType = AlignmentUtil::SEQ_TYPE_PEPTIDE;
			array_push($this->arrTransSeq, $s );
			array_push( $this->arrRawSeq , $s );
			return;
		}
		elseif ($type == AlignmentUtil::SEQ_TYPE_DNA ) { // non-coding DNA
			$this->rawSequenceType = AlignmentUtil::SEQ_TYPE_DNA;
			array_push( $this->arrRawSeq , $s );
			return;
			
		}
		
		array_push( $this->arrRawSeq , $s );
		
		
	}
	
	public function DoAlignment() {
		
		$sFas = ""; //prepare fasta input
		
		if ($this->rawSequenceType == AlignmentUtil::SEQ_TYPE_CODING_DNA) // if is coding DNA, align translated protein sequence
		{
		
			$nCountArrTransSeq = count($this->arrTransSeq);
			for ($i=0;$i<$nCountArrTransSeq ;$i++) {
			
				$sFas .= ">Seq$i".PHP_EOL;
				$sFas .= $this->arrTransSeq[$i].PHP_EOL;
			
			}
		
		}
		else { //align raw sequence.
		
			$nCountArrRawSeq = count($this->arrRawSeq);
			for ($i=0;$i<$nCountArrRawSeq;$i++) {
			
				$sFas .= ">Seq$i".PHP_EOL;
				$sFas .= $this->arrRawSeq[$i].PHP_EOL;
			
			}
		
		}
		
		/*
		echo $this->sTmpFileBase;
		return;*/
		
		mkdir($this->sTmpFileBase, 0777, true);
		fnWriteText($this->sTmpFileBase."in.fas", $sFas);
		
		if ($this->rawSequenceType == AlignmentUtil::SEQ_TYPE_CODING_DNA || $this->rawSequenceType == AlignmentUtil::SEQ_TYPE_PEPTIDE) { //USE clustal omega
		
			$sClustal_Cmd = CLUSTAL_PATH." -i \"".$this->sTmpFileBase."in.fas\" -o \"".$this->sTmpFileBase."out.aln\" --outfmt=clu";
		
		}
		elseif ($this->rawSequenceType == AlignmentUtil::SEQ_TYPE_DNA) { // use clustalw
			
			$sClustal_Cmd = CLUSTALW_PATH." -INFILE=\"".$this->sTmpFileBase."in.fas\" -OUTFILE=\"".$this->sTmpFileBase."out.aln\"";
	
		}
		
		$clustal_ret = exec($sClustal_Cmd);
		
		//echo($clustal_ret );
		$this->parseClustalOutput($this->sTmpFileBase."out.aln");
		
		
		//delete the temp files:
		foreach(glob($this->sTmpFileBase.'*.*') as $v){
			unlink($v);
		}
		
		rmdir($this->sTmpFileBase);
		//finished deleting
		

		$this->mapAlnToRaw();
		
		/*
		print_r($this->arrTransToAlnMap[1]);
		print_r($this->arrAlnToTransMap[1]);
		*/
		
		$this->filterWeakPos();
		
	}
	
	public function GetStrongAreasRaw() {
		
		
		$nCountArrRawSeq = count($this->arrRawSeq);
		
		for($j=0;$j<$nCountArrRawSeq ;$j++) 
		{
		
			
			$sSeq = $this->arrRawSeq[$j];
			
				$nCountArrExcludeMask = count($this->arrExcludeMask);
				for($i=0;$i<$nCountArrExcludeMask ;$i++) 
				{
					if ($this->arrExcludeMask[$i]==0)
					{
						//echo("$i weak\n");
						continue; // weak, exclude
					}
					
					
					if (!isset($this->arrAlnToTransMap[$j][$i])) {
					
					//echo("$j $i \n");
					
						$this->arrStrongAreas[$j] .= ($this->rawSequenceType == AlignmentUtil::SEQ_TYPE_CODING_DNA)? "---":"-";
						continue;
						
					}
					
					//echo("$j $i \n");
					
					$nPos = ($this->rawSequenceType == AlignmentUtil::SEQ_TYPE_CODING_DNA)? $this->arrAlnToTransMap[$j][$i] * 3 : $this->arrAlnToTransMap[$j][$i] ;//- 2;
					
					//$nPos = $nPos==-2? 0:$nPos;
					
					$this->arrStrongAreas[$j] .=  ($this->rawSequenceType == AlignmentUtil::SEQ_TYPE_CODING_DNA)? substr( $sSeq, $nPos, 3) : substr( $sSeq, $nPos, 1);
			

					
				}
		

		}
		
		
		return $this->arrStrongAreas;
		
		
	}
	
	private function parseClustalOutput($sFile) 
	{
		$hFile= @fopen($sFile, "r");
		
		$bFirstLineRead = false;
		$nDataBlockStart = 0;
		
		if ($hFile) {
			while (($sLine = fgets($hFile)) !== false) {
				//echo $sLine;
				//echo(strpos( $sLine, "Seq"));
				
				preg_match('/^([A-Za-z]*)(\d*)\s+(\S*)/', $sLine, $matches );
				
				$bIsSeqLine = ($matches[1] == "Seq");
				$nSeqNum = intval($matches[2]);
				$sSeq = $matches[3];
				//print_r($matches);
				
				/*
				print_r($bIsSeqLine);
				
				continue;
				*/
				
				if ( !$bIsSeqLine && !$bFirstLineRead ) { //haven't reached the first line yet
				
					continue;
				
				}
				
				if ( $bIsSeqLine && !$bFirstLineRead) { //Got to the first line.
				
					$this->arrAlnSeq[$nSeqNum] .= $sSeq;
					
					$nDataBlockStart = strpos($sLine, $sSeq);
					
					//echo($nDataBlockStart);
					
					$bFirstLineRead = true;
					
					continue;
				
				}
				
				if ( $bIsSeqLine && $bFirstLineRead ) {
					
					$this->arrAlnSeq[$nSeqNum] .= $sSeq;
					
					continue;
				}
				
				if ( !$bIsSeqLine && $bFirstLineRead) {
					
					if (strlen($sLine) <= $nDataBlockStart) { // this is an empty line
						continue;
					}
					
					$sCon = substr($sLine, $nDataBlockStart );
					
					$sCon = preg_replace('/\n/', "", $sCon);
					$sCon = preg_replace('/\r/', "", $sCon);
					
					$this->sConsensusLn .= $sCon;
					
				}
				
				
			}
			if (!feof($hFile)) {
				return false;
			}
			fclose($hFile);
		}
		
		/*
		print_r($this->sConsensusLn);
		print_r($this->arrAlnSeq);
		*/
		
	}
	
	private function mapAlnToRaw() 
	{
	
		$nCountArrAlnSeq = count($this->arrAlnSeq);
		for ($i=0;$i<$nCountArrAlnSeq;$i++) { // go over each sequence
		
			$arrAlnSeq = str_split($this->arrAlnSeq[$i]);
			
			if ($this->rawSequenceType == AlignmentUtil::SEQ_TYPE_CODING_DNA )
			{
				$arrTransSeq = str_split($this->arrTransSeq[$i]); // use translated sequence
			}
			else { // non coding dna and plain peptide input
				$arrTransSeq = str_split($this->arrRawSeq[$i]); //use raw sequence
			}
			
			$arrMap = array();
			
			$nIndexTransSeq = 0;
			$nIndexAlnSeq = 0;
			$nLenTransSeq = count($arrTransSeq);
			$nLenAlnSeq = count($arrAlnSeq);
			
			while($nIndexTransSeq < $nLenTransSeq && $nIndexAlnSeq < $nLenAlnSeq ) {
			
				if ($arrTransSeq[$nIndexTransSeq] == $arrAlnSeq[$nIndexAlnSeq] ) 
				{ //if the two characters the same, record index:
				
					$arrMap[$nIndexTransSeq] = $nIndexAlnSeq;
					$nIndexTransSeq++;
					$nIndexAlnSeq++;
					continue;
				
				}
				
				$nIndexAlnSeq++;
				continue;
				
			
			}
			
			//print_r($arrMap);
			$this->arrTransToAlnMap[$i] = $arrMap;
			$this->arrAlnToTransMap[$i] = array_flip($arrMap);
		
		}
		
		
	}

	public function GetAlignmentScore() {
	
		//$ScoreMatrix = array("*" => 10, ":"=>4, "." => 1, " " => -2); // highest score is 10!! 10 8 5 0
	

		$nScore = 0;
		
		$arrConsensus = str_split($this->sConsensusLn);
		$nConsensusLen = count($arrConsensus);
		
		/*
		foreach ($this->SumScoreMatrix as $char => $score) {
			$nScore += substr_count($this->sConsensusLn , $char) * $score ;
		}
		*/
		
		for ($i=0;$i<$nConsensusLen;$i++) {
		
			if ($this->arrExcludeMask[$i] == 0) {
				continue; // this site is excluded, continue;
			}
			
			$nScore += $this->SumScoreMatrix[ $arrConsensus[$i] ];
		
		}
		
		return $nScore;
		
	}
	
	
	private function filterWeakPos() {
	
		$arrConsensus = str_split($this->sConsensusLn);
		$nConsensusLen = count($arrConsensus);
		
		$this->arrExcludeMask = array_fill(0, $nConsensusLen , 1); //fill the mask with 1
		
		for($i=0;$i<$nConsensusLen;$i++)
		{
		
			$nStart = ( ($i - $this->nLookLeft) >=0)? $i - $this->nLookLeft:0;
			$nEnd = ( ($i + $this->nLookRight) < $nConsensusLen)? ($i + $this->nLookRight) : ($nConsensusLen -1);
			
			$nTotalLen = $nEnd - $nStart + 1;
			$nHighestScore = (1 + $this->nLookLeft + $this->nLookRight) * 10;
			$nScore=0;
			
			if (($i - $this->nLookLeft) < 0) {
				$nScore += ($i - $this->nLookLeft) * $this->ScoreMatrix[" "] ; //if not exist assume to be gap
			}
			
			if (($i + $this->nLookRight) > $nConsensusLen) {
				$nScore += (($i + 1 + $this->nLookRight) - $nConsensusLen) * $this->ScoreMatrix[" "];
			}
			
			for ($j=$nStart;$j<=$nEnd;$j++) 
			{
				$nScore += $this->ScoreMatrix[ $arrConsensus[$j] ];
			}
			
			$nPercentScore = $nScore  / $nHighestScore;
			
			//echo($nStart." ".$nEnd." ".$nPercentScore."\t");
			
			if ( $nPercentScore >= $this->nScoreCutoff) {
				$this->arrExcludeMask[$i] = 1;
				
				if ($this->bOnlyFilterWeakFlanks) { // if only remove the weak flanks, don't continue.
					break;
				}
				
			}
			else {
				$this->arrExcludeMask[$i] = 0;
			}
			
			
		
		}
		
		if ( ! $this->bOnlyFilterWeakFlanks) {
			return;
		}
		
		
		for($i=$nConsensusLen-1;$i>=0;$i--) //do it from the other direction
		{
		
			$nStart = ( ($i - $this->nLookLeft) >=0)? $i - $this->nLookLeft:0;
			$nEnd = ( ($i + $this->nLookRight) < $nConsensusLen)? ($i + $this->nLookRight) : ($nConsensusLen -1);
			
			$nTotalLen = $nEnd - $nStart + 1;
			$nHighestScore = (1 + $this->nLookLeft + $this->nLookRight) * 10;
			$nScore=0;
			
			if (($i - $this->nLookLeft) < 0) {
				$nScore += ($i - $this->nLookLeft) * $this->ScoreMatrix[" "] ; //if not exist assume to be gap
			}
			
			if (($i + $this->nLookRight) > $nConsensusLen) {
				$nScore += (($i + 1 + $this->nLookRight) - $nConsensusLen) * $this->ScoreMatrix[" "];
			}
			
			for ($j=$nStart;$j<=$nEnd;$j++) 
			{
				$nScore += $this->ScoreMatrix[ $arrConsensus[$j] ];
			}
			
			$nPercentScore = $nScore  / $nHighestScore;
			
			//echo($nStart." ".$nEnd." ".$nPercentScore."\t");
			
			if ( $nPercentScore >= $this->nScoreCutoff) {
				$this->arrExcludeMask[$i] = 1;
				
				if ($this->bOnlyFilterWeakFlanks) { // if only remove the weak flanks, don't continue.
					break;
				}
				
			}
			else {
				$this->arrExcludeMask[$i] = 0;
			}
			
			
		
		}
		
		//print_r($this->arrExcludeMask);
		
	}
	
	public static function RemoveUnreliableFlankingRegions($sSeq1, $sSeq2) //input must be aligned to correct codons already and of equal length 
	{
		 
		if (strlen($sSeq1)!= strlen($sSeq2)) {
			return false;
		}
		
		if ( strlen($sSeq1) % 3 !=0) { //must be aligned!
			return false;
		}
		
		$nTotalLen = strlen($sSeq1);
		$nNewStart = 0;
		$nNewEnd = $nTotalLen-1 - 2;
		
		for ($i=0;$i<=$nTotalLen-3;$i+=3) 
		{
			$Codon1 = substr($sSeq1 , $i, 3);
			$Codon2 = substr($sSeq2 , $i, 3);
			
			$Codon1Simple = ( substr_count($Codon1, substr($Codon1, 0,1)) == 3 );
			$Codon2Simple = ( substr_count($Codon2, substr($Codon2, 0,1)) == 3 );
			
			if ( $Codon1 == $Codon2 && !$Codon1Simple && !$Codon2Simple)
			{
				//found equal codon
				$nNewStart = $i;
				break;
			}
		}
		
		for ($i=$nTotalLen-3; $i>=0;$i-=3) 
		{
			$Codon1 = substr($sSeq1 , $i, 3);
			$Codon2 = substr($sSeq2 , $i, 3);
			
			$Codon1Simple = ( substr_count($Codon1, substr($Codon1, 0,1)) == 3 );
			$Codon2Simple = ( substr_count($Codon2, substr($Codon2, 0,1)) == 3 );
			
			if ( $Codon1 == $Codon2 && !$Codon1Simple && !$Codon2Simple)
			{
			
				//found equal codon
				$nNewEnd = $i;
				break;
			}
		}
		
		$sNewSeq1 = substr($sSeq1 , $nNewStart , $nNewEnd + 3 - $nNewStart );
		$sNewSeq2 = substr($sSeq2 , $nNewStart , $nNewEnd + 3 - $nNewStart );
		
		return array($sNewSeq1, $sNewSeq2);
	}
	
	
}

function fnDNA2Peptide($DNA, $bNoStop = false) {
	global $CodonMap;

	$DNA = strtoupper($DNA);
	$Peptide = "";
	$DNA2 = "";
	
	/*
	if (strlen($DNA) % 3 !=0 ) {
		echo("sequence error.");
	}
	*/
	
	for ($i=0;$i<=strlen($DNA)-3;$i+=3) {
	
		$Codon = substr($DNA, $i, 3);
		
		if (strlen($Codon)<3) {
			continue;
		}
		
		if (strpos($Codon, "N") !== false) { // if the codon contains uncertain nucleotide
		
			continue; //throw it away!
		}
		
		//echo($Codon." ");
		if(!isset($CodonMap[$Codon]))
		{
			if ($bNoStop) {
				continue;
			}
		}
		
		
		
		$AminoAcid = $CodonMap[$Codon];
		
		if ($bNoStop && $AminoAcid=="*") { // if throw away stop codon
			continue;
		}
		
		$Peptide .= $AminoAcid;
		$DNA2 .= $Codon;
		
	}
	
	return array( $Peptide, $DNA2 ) ;

}

function fnSeqType($peptide) {

	if (strpos($peptide, "*") !== false) {
		return 2; //pseudogene
	}
	else {
		return 1;
	}
}

function fnReverseComplement($sIn) {

	//$sRet = strrev($sIn);
	/*
	$arrIn = array_reverse(str_split($sIn));
	
	
	
	for ($i=0;$i<=count($arrIn)-1;$i++) {
		$sChar = $arrIn[$i];
		switch($sChar) {
			case "A" : $sChar="T"; break;
			case "T" : $sChar="A"; break;
			case "C" : $sChar="G"; break;
			case "G" : $sChar="C"; break;
			case "U" : $sChar="A"; break;
			case "-" : $sChar="-"; break;
			case "N" : $sChar="N"; break;
			
		}
		
		$sRet .= $sChar;
		
		
		
	}
	*/
	/*
	$sRet = preg_replace('/A/i', '1', $sRet);
	$sRet = preg_replace('/T/i', '2', $sRet);
	$sRet = preg_replace('/C/i', '3', $sRet);
	$sRet = preg_replace('/G/i', '4', $sRet);
	$sRet = preg_replace('/1/', 'T', $sRet);
	$sRet = preg_replace('/2/', 'A', $sRet);
	$sRet = preg_replace('/3/', 'G', $sRet);
	$sRet = preg_replace('/4/', 'C', $sRet);
	
	return $sRet;
	*/
	
	return strrev( strtr($sIn, 'ACBDGHKMNSRUTWVYacbdghkmnsrutwvy', 'TGVHCDMKNSYAAWBRTGVHCDMKNSYAAWBR'));
	
	
}

function fnCleanDNASeq($s) {
	$s = strtoupper($s);
	$s = preg_replace('/[\n\r\-]/', "", $s);
	//$s = preg_replace('/\r/', "", $s);
	//$s = preg_replace('/\-/', "", $s);
	
	return $s;
}

function fnDenovoPredictORF($sDNA) {

	$sDNA = fnCleanDNASeq($sDNA);
	$arrFrames = array(1, 2, 3, -1, -2, -3);
	$sRetDNAORF = "";
	$sRetPeptideORF = "";
	
	$arrCandidateORF = array();
	
	$nCountArrFrames = count($arrFrames);
	
	for ($f=0;$f<$nCountArrFrames;$f++) {
	
		$nFrame = $arrFrames[$f];
		$sSeq = $sDNA;
		
		if ($nFrame < 0) {
			$sSeq = fnReverseComplement($sSeq);
			$nStart = (0 - $nFrame) - 1;
		}
		else {
			$nStart = $nFrame - 1;
		}
		
		$sSeq = substr($sSeq, $nStart );
		
		$arrRet = fnDNA2Peptide($sSeq, false);//preserve stop codons
		$sPeptide = $arrRet[0]; 
		$sSeq = $arrRet[1];
		
		$arrStops = fnStringPositions($sPeptide, "*");
		$arrStarts = fnStringPositions($sPeptide, "M");
		
		$arrCandidateORF[$f] = fnFindLongestReadingFrame(strlen($sPeptide), $arrStarts, $arrStops);
		
		$arrCandidateORF[$f][3] = substr($sPeptide, $arrCandidateORF[$f][0], $arrCandidateORF[$f][2] );
		$arrCandidateORF[$f][4] = substr($sSeq, $arrCandidateORF[$f][0]*3, $arrCandidateORF[$f][2] * 3 );

	}	
	
	//print_r($arrCandidateORF);
	
	//Now go over the candidate ORFs, get the one significantly longer than others:
	
	$arrSignficantLongerORF = array();
	
	$nCountArrCandidateORF = count($arrCandidateORF);
	
	for ($i=0;$i<$nCountArrCandidateORF;$i++) {
		
		$arrCurrORF = $arrCandidateORF[$i];
		$arrCandidateORFCopy = $arrCandidateORF;
		unset($arrCandidateORFCopy[$i]); //remove it.
		$arrOtherORFs = array_values($arrCandidateORFCopy);
		
		$arrOtherORFLens = array();
		
		$nCountArrOtherORFs = count($arrOtherORFs);
		for ($j=0;$j<$nCountArrOtherORFs;$j++) {
			$arrOtherORFLens[] = $arrOtherORFs[$j][2];
		}
		
		$nAvg = average($arrOtherORFLens);
		$nStd = standard_deviation($arrOtherORFLens);
		
		if ($arrCurrORF[2] > ($nAvg+$nStd) )
		{
			$arrSignficantLongerORF[] = $arrCurrORF;
		}
		
	}
	
	if (count($arrSignficantLongerORF) == 1) 
	{
		return array($arrSignficantLongerORF[0][3], $arrSignficantLongerORF[0][4] );
	}
	else {
		return false; // ORF not found or several long ORFs found
	}

}

class Bcftools{ //Bcftools wrapper

	private $sBcfFile;
	private $arrBuffer;
	private $sBufferScfld;
	private $nBufferStart;
	private $nBufferEnd;
	private $bIsError;

	public function SetBcf($sPath) {
		
		if (! file_exists($sPath)) {
			die("Cannot find bcf file: $sPath !");
		}

		if (! file_exists("$sPath.csi")) {
			echo("Building bcf index for $sPath ...\n");
			exec(BCFTOOLS_PATH." index $sPath" );
			echo("Done.\n");

		}

		$this->sBcfFile = $sPath;
		$this->arrBuffer = array();
		$this->sBufferScfld = '';
		$this->nBufferStart = -1;
		$this->nBufferEnd = -1;
		$this->bIsError = false;
	}

	public function GetVariantsAt($sChr, $nIndex) { //index is 1 based

		//echo("$sChr: $nIndex\n");
	
		if ($sChr != $this->sBufferScfld || $nIndex < $this->nBufferStart || $nIndex > $this->nBufferEnd) { 
			// update buffer
			$this->arrBuffer = null;
			$this->arrBuffer = array();
			$this->bIsError = false;
/*
			$arrRet = array();
			$nReturnVal = -1;
			$sCmd = BCFTOOLS_PATH." view -N $this->sBcfFile ".$sChr.':'.$nIndex.'-'.($nIndex + BCFTOOLS_BUFFERSIZE - 1);
			$descriptorspec = array(
			   0 => array("pipe", "r"),   // stdin is a pipe that the child will read from
			   1 => array("pipe", "w"),   // stdout is a pipe that the child will write to
			   2 => array("pipe", "w")    // stderr is a pipe that the child will write to
			);
			flush();
			ob_implicit_flush(true);
			ob_end_flush(); 
			$process = proc_open($sCmd, $descriptorspec, $pipes, './', array());
			$bBCFOutputVerified = false;

			if (is_resource($process)) {
			    while ($sLn = fgets($pipes[1])) {

				if ($sErr = fgets($pipes[2])) {
					if (strlen($sErr)>1) {
						echo("Warning: Region $sChr:$nIndex error in bcf ".$this->sBcfFile.". $sErr\n");
						proc_close($process);
						return -1;
					}
				}

				if (substr($sLn, 0, 1) != '#') {
					$arrFields = array();
	

	//		Scfld0	9993	.	A	G	120	.	DP=17;AF1=1;AC1=2;DP4=0,0,9,7;MQ=60;FQ=-60.1	GT:PL:GQ	1/1:138,48,0:93

	

					if (strlen($sLn) > 4000 ) {
						continue;
					}

        				if (0 == preg_match('/^(\S+)\t(\S+)\t\S+\t(\S+)\t(\S+)\t(\S+)\t\S\tDP=([0-9]+);\S/' , $sLn, $arrFields)) { // no match
						continue;
					}
					//print_r($arrFields); die();

					if (!$bBCFOutputVerified) {
						if ($arrFields[1] != $sChr) { // BCFTools will return stuff from the begining if it doesn't find the target region. 
							echo("Warning: Region $sChr:$nIndex not found in bcf ".$this->sBcfFile.".\n");
							return false;
						}
						else {
							$bBCFOutputVerified = true;
						}
					}

					$this->arrBuffer[$arrFields[2]] = array_slice($arrFields, 3);
				}
				flush();
			    }

			    if ( (proc_get_status($process))["exitcode"] != 0) {
				echo("Warning: Error occured when trying to pull out buffer region from bcf file with bcftools. ".$this->sBcfFile." : $sChr : $nIndex\n");
				proc_close($process);
				return -1;
			    }
			    proc_close($process);
			}

			
*/
			
			$arrRet = array();
			$nReturnVal = -1;

			//echo(BCFTOOLS_PATH." view -H $this->sBcfFile ".$sChr.':'.$nIndex.'-'.($nIndex + BCFTOOLS_BUFFERSIZE - 1). " 2>/dev/null | grep \"^[^\\#]\" ".PHP_EOL. memory_get_usage() . PHP_EOL);
			exec(BCFTOOLS_PATH." view -H $this->sBcfFile ".$sChr.':'.$nIndex.'-'.($nIndex + BCFTOOLS_BUFFERSIZE - 1). " 2>/dev/null | grep \"^[^\\#]\" ", $arrRet, $nReturnVal);
			
			if ($nReturnVal != 0 ) {
				echo("Warning: Error occured when trying to pull out buffer region from bcf file with bcftools. ".$this->sBcfFile." : $sChr : $nIndex\n");
				$this->bIsError = true;
				$this->sBufferScfld = $sChr;
				$this->nBufferStart = $nIndex;
				$this->nBufferEnd = $nIndex + BCFTOOLS_BUFFERSIZE - 1;
				return -1;
			}


			$bBCFOutputVerified = false;
			foreach( $arrRet as $sLn ) {
				if (substr($sLn, 0, 1) != '#') {
					$arrFields = array();
	

	//		Scfld0	9993	.	A	G	120	.	DP=17;AF1=1;AC1=2;DP4=0,0,9,7;MQ=60;FQ=-60.1	GT:PL:GQ	1/1:138,48,0:93

	

					if (strlen($sLn) > 6000 ) {
						continue;
					}

        				if (0 == preg_match('/^(\S+)\t(\S+)\t\S+\t(\S+)\t(\S+)\t(\S+)\t\S+\t(\S+)\t(\S*)\t?(\S*)/' , $sLn, $arrFields)) { // no match
						continue;
					}
					//print_r($arrFields); die();

					if (!$bBCFOutputVerified) {
						if ($arrFields[1] != $sChr) { // BCFTools will return stuff from the begining if it doesn't find the target region. 
							echo("Warning: Region $sChr:$nIndex not found in bcf ".$this->sBcfFile.".\n");
							return false;
						}
						else {
							$bBCFOutputVerified = true;
						}
					}
					$arrFieldsParsed = fnParseVCFTags($arrFields); //tags will be appended to the end of the array, for example, DP=1 becomes 'DP' => 1, GT=1/1 becomes 'GT' => 1/1
					$this->arrBuffer[$arrFieldsParsed[2]] = array_slice($arrFieldsParsed, 3);
				}
			}
			

			$this->sBufferScfld = $sChr;
			$this->nBufferStart = $nIndex;
			$this->nBufferEnd = $nIndex + BCFTOOLS_BUFFERSIZE - 1;
		
		}

		if ($this->bIsError) { //if there's an error
			return -1;
		} else if (array_key_exists( $nIndex, $this->arrBuffer)) {
			return $this->arrBuffer[$nIndex];
		}
		else {
			return false;
		}
	
	}

}

function fnParseVCFTags($arrFields) {
	$arrRet = $arrFields;
	$sInfo = (array_key_exists(6, $arrFields))? $arrFields[6]:""; // in format KEY1=VALUE1;KEY1=VALUE1; "=VALUE" is optional
	$sExtraKeys = (array_key_exists(7, $arrFields))? $arrFields[7]:""; // in format KEY1:KEY2:KEY3 
	$sExtraVals = (array_key_exists(8, $arrFields))? $arrFields[8]:""; // in format VALUE1:VALUE2:VALUE3
	
	if ($sInfo != "") {
		$arrInfo = explode(";" , $sInfo);
		foreach($arrInfo as $sInfoField) {
			$arrInfoField = explode("=", $sInfoField);
			if (count($arrInfoField) == 2) {
				$arrRet[$arrInfoField[0]] = $arrInfoField[1];
			}
		}
	}

	if ($sExtraKeys!="" && $sExtraVals!="") {
		$arrExtraKeys = explode(":" , $sExtraKeys);
		$arrExtraVals = explode(":" , $sExtraVals);
		if (count($arrExtraKeys)  == count($arrExtraVals) ) {
			for($k=0;$k<count($arrExtraKeys);$k++) {
				$arrRet[ $arrExtraKeys[$k] ] = $arrExtraVals[$k] ;
			}
		}
	}

	return $arrRet;
}

class FastaHack{ //fastahack wrapper

	private $sFastaDb;
	private $sBuffer;
	private $sBufferScfld;
	private $nBufferStart;
	private $nBufferEnd;
	private $arrContigIndex;
	
	public function SetDb($db) {
	
		if (! file_exists($db)) {
			die("Cannot find fasta database file: $db !");
		}

		
		$this->sFastaDb = $db;
		$this->sBuffer = '';
		$this->sBufferScfld = '';
		$this->nBufferStart = -1;
		$this->nBufferEnd = -1;
	
	}
	
	private function LoadFastaIndex() {
		$hFai = fopen($this->sFastaDb . '.fai' , r);
		
	}

	public function GetCDSRegion($sChr, $nStart, $nEnd) {
	
		if ($nStart > $nEnd) {
			die("Fastahack error: start index > end index");
		}
	
		if ($sChr != $this->sBufferScfld || $nStart < $this->nBufferStart || $nEnd > $this->nBufferEnd) { 
			// update buffer
			$this->sBufferScfld = $sChr;
			$this->nBufferStart = $nStart;
			$this->nBufferEnd = $nStart + FASTAHACK_BUFFERSIZE - 1;
			
			
			$arrRet = array();
			$nReturnVal = -1;
			exec(FASTAHACK_PATH." $this->sFastaDb -r ".$this->sBufferScfld.':'.$this->nBufferStart.'..'.$this->nBufferEnd, $arrRet, $nReturnVal);
			//$sLn = shell_exec(FASTAHACK_PATH." $this->sFastaDb -r ".$this->sBufferScfld.':'.$this->nBufferStart.'..'.$this->nBufferEnd);
			if ($nReturnVal != 0 ) {
				die("Error occured when trying to pull out buffer region from fasta file with fastahack.");
			}
		
			foreach( $arrRet as $sLn ) {
				if (strpos($sLn, "index file") === false) {
					$this->sBuffer = $sLn;
					break;
				}
			}
		
		}
		
		return substr($this->sBuffer, $nStart - $this->nBufferStart, $nEnd - $nStart + 1);
		

	}

	public function GetContig($sChr) {
	

			$arrRet = array();
			$nReturnVal = -1;
			exec(FASTAHACK_PATH." $this->sFastaDb -r ".$sChr, $arrRet, $nReturnVal);
			//$sLn = shell_exec(FASTAHACK_PATH." $this->sFastaDb -r ".$this->sBufferScfld.':'.$this->nBufferStart.'..'.$this->nBufferEnd);
			if ($nReturnVal != 0 ) {
				die("Error occured when trying to pull out buffer region from fasta file with fastahack.");
			}
		
			foreach( $arrRet as $sLn ) {
				if (strpos($sLn, "index file") === false) {
					return $sLn;
				}
			}

		
	}

}

class GenemarkGTF {

	private $sGTF; // gtf file;
	public $arrGenes;
	
	public function LoadGTF($sGTFFile) {
	
		if (! file_exists($sGTFFile)) {
			die("Cannot find GTF file: $sGTFFile");
		}
		
		$this->sGTF = $sGTFFile;
		
		$hFile = fopen($this->sGTF , "r");
		
		if (!$hFile) {
			die("Cannot open GTF file for reading: $sGTFFile");
		}
		
		$this->arrGenes = array();
		
		while (($sLn = fgets($hFile)) !== false) {
			$arrFields = array();
        	preg_match('/(\S+)\t(\S+)\t(\S+)\t([0-9]+)\t([0-9]+)\t(\S+)\t(\S)\t([0-2])\tgene_id \"(\S+)\"/' , $sLn, $arrFields);
        	
        	/*
        	Array
(
    [0] => Scfld0       GeneMark.hmm    start_codon     1173    1175    .       +       0       gene_id "1_g"
    [1] => Scfld0   //chromosome name
    [2] => GeneMark.hmm
    [3] => start_codon  // type
    [4] => 1173 //start
    [5] => 1175 //end
    [6] => .    // dont know what
    [7] => + // orientation
    [8] => 0    // frame
    [9] => 1_g  // gene name
)
        	
        	*/
        	//print_r($arrFields);
        	//die();
        	$sGeneName = $arrFields[9];
        	
        	if (!array_key_exists($sGeneName, $this->arrGenes)) {
        		$this->arrGenes[$sGeneName] = array(); //Array for each CDS.
        	}
        	
        	if ($arrFields[3] != "CDS") {
        		continue; // only process CDS
        	}
        	
        	$arrCDS = array($arrFields[1], $arrFields[4] , $arrFields[5], $arrFields[7], $arrFields[8]); //chr, start index, end index, orientation, frame
    		array_push($this->arrGenes[$sGeneName] , $arrCDS );
    	
    	}
    	
    	print_r($this->arrGenes["100_g"]);
		
		fclose($hFile);
	}

}

class AUTest {
	private $sAlignment; //must be in phylip format
	private $arrTrees; //an array of alternative topologies in newick format
	private $nInstanceID ;
        
        const IN_FILE_NAME = "in.phy";
	const RAXML_OUT = "output";
        const TREE_FILE_NAME = "in.tre";
	const RAXML_RET_PREFIX = "RAxML_perSiteLLs";


	public function SetAlignment($sAlg) {
		$this->sAlignment = $sAlg;
	}

	public function SetTrees(&$arrT) {
		$this->arrTrees = $arrT;
	}

	public function DoTest() {

 		$this->nInstanceID = rand( 1000000 , 9999999) ;

		try {

			if (!file_exists(AU_TEMP_PATH)) {
				mkdir(AU_TEMP_PATH);
			}
		}
		catch(Exception  $e) 
		{}

		$sTempDir = AU_TEMP_PATH . "/". $this->nInstanceID . "/";
		mkdir($sTempDir);
		chdir($sTempDir);

		$hAln = fopen($sTempDir."/".AUTest::IN_FILE_NAME, "w");
		fwrite($hAln , $this->sAlignment);
		fclose($hAln);

		$hTrees = fopen($sTempDir."/".AUTest::TREE_FILE_NAME, "w");

		foreach( $this->arrTrees as $sTree) {
			fwrite($hTrees, $sTree.PHP_EOL);
		}
		fclose($hTrees);

		//Run RAxML

		$sCmd = RAXML_PATH . " -f g -s " . $sTempDir."/".AUTest::IN_FILE_NAME . " -m GTRGAMMA -z ". $sTempDir."/".AUTest::TREE_FILE_NAME ." -n ".AUTest::RAXML_OUT;
		exec( $sCmd );

		$sRaxmlLnOut = AUTest::RAXML_RET_PREFIX.".".AUTest::RAXML_OUT;
		if (!file_exists($sRaxmlLnOut )) {
			echo("RAxML did not generate output $sRaxmlLnOut! ");
			return false;
		}

		//now do consel
		exec(CONSEL_PATH."/seqmt --puzzle RAxML_perSiteLLs.output");
		exec(CONSEL_PATH."/makermt ".AUTest::RAXML_RET_PREFIX);
		exec(CONSEL_PATH."/consel ".AUTest::RAXML_RET_PREFIX);
		exec(CONSEL_PATH."/catpv -s 1 ".AUTest::RAXML_RET_PREFIX." > pvals.txt");
		
		$sPvals = file_get_contents("pvals.txt");
		$arrPvalLns = explode("\n", $sPvals);

		$arrRet = array();

		for ($i=4;$i<count($arrPvalLns); $i++) {

			if (substr($arrPvalLns[$i], 0, 1) != "#") {
				continue;
			}
			$arrMatches = array();
			preg_match('/\#\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+\|\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+\|/', $arrPvalLns[$i], $arrMatches); 
			/*#    2    2  190.9  2e-07  3e-07 |      0  1e-83      0      0      0      0 |*/
			
			array_shift($arrMatches);
			array_push($arrRet , implode("\t", $arrMatches) );
		}

		chdir(dirname(__FILE__) );

		rrmdir($sTempDir);//remove dir

		return $arrRet;

	}
}

class CodemlUtil {

        private $nInstanceID ;
        private $sTmpFileBase;
        private $sTreeContent;
        private $arrLabels = array();
        private $arrSeqs = array();
        private $arrRet = array();
        
        private $bInit = false;
        
        const IN_FILE_NAME = "in.fas";
        const TREE_FILE_NAME = "in.tre";
        const PREFIX_CTRL_FILE1 = "codeml.1";
        const PREFIX_CTRL_FILE2 = "codeml.2";
        const PREFIX_CTRL_FILE3 = "codeml.3";
        
        
        private function initialize() {
        
          if ($this->bInit) {return;};
          
          
      $this->nInstanceID = rand( 1000000 , 9999999) ;
          $this->sTmpFileBase = CODEML_TMP."/".$this->nInstanceID."/";
          
          mkdir($this->sTmpFileBase, 0777, true); // create temp file folder
          
          //Read the template control file:
          $sCtrlTemplate = fnReadText(CODEML_CTRL_TEMPLATE);
          
          
          //Common variables shared by all control files:
          $sCtrlTemplate = str_replace("%SEQ_FILE%", $this->sTmpFileBase . CodemlUtil::IN_FILE_NAME ,  $sCtrlTemplate);
          $sCtrlTemplate = str_replace("%TREE_FILE%",  $this->sTmpFileBase . CodemlUtil::TREE_FILE_NAME  ,  $sCtrlTemplate);
          
          //Control file 1 for just getting plain dN/dS:
          $sCtrl1 = str_replace("%OUT_FILE%",  $this->sTmpFileBase . CodemlUtil::PREFIX_CTRL_FILE1.".out"  ,  $sCtrlTemplate);
          $sCtrl1 = str_replace("%NSSITES%", "0"  ,  $sCtrl1);
          $sCtrl1 = str_replace("%RUN_MODE%", "-2"  ,  $sCtrl1); //Run mode pairwise
          $sCtrl1 = str_replace("%FIX_OMEGA%", "0"  ,  $sCtrl1); //fix omega
          $sCtrl1 = str_replace("%FIX_ALPHA%", "1"  ,  $sCtrl1); //fix gamma
          
          fnWriteText($this->sTmpFileBase .  CodemlUtil::PREFIX_CTRL_FILE1 . ".ctl" , $sCtrl1);
          
          //Control file 2 for M1, M2, M7, M8
          $sCtrl2 = str_replace("%OUT_FILE%",  $this->sTmpFileBase . CodemlUtil::PREFIX_CTRL_FILE2.".out"  ,  $sCtrlTemplate);
          $sCtrl2 = str_replace("%NSSITES%", "1 2 7 8"  ,  $sCtrl2);
          $sCtrl2 = str_replace("%RUN_MODE%", "0"  ,  $sCtrl2); //Run mode user tree
          $sCtrl2 = str_replace("%FIX_OMEGA%", "0"  ,  $sCtrl2); //estimate omega
          $sCtrl2 = str_replace("%FIX_ALPHA%", "1"  ,  $sCtrl2); //fix gamma
          
          fnWriteText($this->sTmpFileBase .  CodemlUtil::PREFIX_CTRL_FILE2 . ".ctl", $sCtrl2);
          
          //Control file 3 for M8a
          $sCtrl3 = str_replace("%OUT_FILE%",  $this->sTmpFileBase . CodemlUtil::PREFIX_CTRL_FILE3.".out"  ,  $sCtrlTemplate);
          $sCtrl3 = str_replace("%NSSITES%", "8"  ,  $sCtrl3);
          $sCtrl3 = str_replace("%RUN_MODE%", "0"  ,  $sCtrl3); //Run mode user tree
          $sCtrl3 = str_replace("%FIX_OMEGA%", "1"  ,  $sCtrl3); //estimate omega
          $sCtrl3 = str_replace("%FIX_ALPHA%", "1"  ,  $sCtrl3); //fix gamma
          
          fnWriteText($this->sTmpFileBase .  CodemlUtil::PREFIX_CTRL_FILE3 . ".ctl", $sCtrl3);
          
          chdir($this->sTmpFileBase ); //change working dir
          
          $this->bInit = true;
          
          
        }
        
        public function AddSequence($Label , $Seq)  //the sequence is already well aligned so no need to do anything else
        {
                $this->arrLabels[] = $Label;
                $this->arrSeqs[] = $Seq;
        }
        
        public function SpecifyTree($sNewickTree)
        {
                $this->sTreeContent = $sNewickTree;
        }
        
        public function RunCodeml($bOnlyGetdNdS) {
        
                if (count($this->arrRet) > 0)
                {
                        return $this->arrRet;
                }
                
                $this->initialize();
                // Prepare the in files:
                $sFas = "";
                for ($i=0;$i<count($this->arrSeqs);$i++) {
                        
                                $sFas .= '>'.trim($this->arrLabels[$i])."\n";
                                $sFas .= trim($this->arrSeqs[$i])."\n";
                                                        
                }

                
                fnWriteText($this->sTmpFileBase . CodemlUtil::IN_FILE_NAME, $sFas);
                
                if (!isset($this->sTreeContent) ) {
                        $sTree = "(";
                        for ($i=0;$i<count($this->arrSeqs);$i++) {
                                                       
                                $sTree .= ($i==0)? $this->arrLabels[$i]: ", ". $this->arrLabels[$i];
                                $sTree .= ($i == (count($this->arrSeqs)-1) )? ");" : "";
                        
                        }
                        fnWriteText($this->sTmpFileBase . CodemlUtil::TREE_FILE_NAME, $sTree );
                }
                else {
                	
                        $sTree = $this->sTreeContent;
                        if (substr(trim($sTree), -1, 1) != ';') {
                        	    $sTree .= ';';             	
                        }
                        fnWriteText($this->sTmpFileBase . CodemlUtil::TREE_FILE_NAME, $sTree );
                }
                
                
                
                //Run Ctrl file 1: ========================================
                
                echo("Running Model 0 ...\n");
                $sCMD = CODEML_PATH . " ". $this->sTmpFileBase .  CodemlUtil::PREFIX_CTRL_FILE1 . ".ctl";
                //echo($sCMD);
                
                //return;
                
                $sRet = exec($sCMD);
                $sF1 = fnReadText($this->sTmpFileBase .  CodemlUtil::PREFIX_CTRL_FILE1 . ".out");
                
                //Parse the file:
                /* $this->arrRet[]; */
                //echo($sF1);
                preg_match_all('/t\=\s*\S+\s+S\=\s*\S+\s*N\=\s*\S+\s*dN\/dS\=\s*([0-9]*\.[0-9]+)\s*dN\=\s*([0-9]*\.[0-9]+)\s*dS\=\s*(([0-9]*\.[0-9]+))/m', $sF1, $matches, PREG_PATTERN_ORDER);
                //t=50.0000  S=    46.9  N=   178.1  dN/dS= 0.0086  dN= 0.6661  dS=77.3772

                //print_r($matches);
                
                $this->arrRet["dN"] = $matches[2][0];
                $this->arrRet["dS"] = $matches[3][0];
                $this->arrRet["dNdS"] = $matches[1][0];
                
                if ($bOnlyGetdNdS) {
                        return $this->arrRet;
                }
                
                
                //Run Ctrl file 2 ======================================
                echo("Running Model 1 2 7 8 ...\n");
                
                $sCMD = CODEML_PATH . " ". $this->sTmpFileBase .  CodemlUtil::PREFIX_CTRL_FILE2 . ".ctl";
                //echo("$sCMD\n");
                $sRet = exec($sCMD);
                $sF2 = fnReadText($this->sTmpFileBase .  CodemlUtil::PREFIX_CTRL_FILE2 . ".out");
                
                //Parse the file:
                /* $this->arrRet[]; */
                //echo($sF2);
                preg_match_all('/lnL\(ntime\:\s*\S*\s*np\:\s*\S*\)\:\s+(\-[0-9]*\.[0-9]*)/m', $sF2, $matches, PREG_PATTERN_ORDER);
                
                //print_r($matches);
        
                $this->arrRet["M1"] = $matches[1][0];
                $this->arrRet["M2"] = $matches[1][1];
                $this->arrRet["M7"] = $matches[1][2];
                $this->arrRet["M8"] = $matches[1][3];
                
                //Run Ctrl file 3 for M8a ================================
                echo("Running Model 8a ...\n");
                
                $sCMD = CODEML_PATH . " ". $this->sTmpFileBase .  CodemlUtil::PREFIX_CTRL_FILE3 . ".ctl";
                //echo("$sCMD\n");
                $sRet = exec($sCMD);
                $sF3 = fnReadText($this->sTmpFileBase .  CodemlUtil::PREFIX_CTRL_FILE3 . ".out");
                
                //Parse the file:
                /* $this->arrRet[]; */
                //echo($sF3);
                
                preg_match_all('/lnL\(ntime\:\s*\S*\s*np\:\s*\S*\)\:\s+(\-[0-9]*\.[0-9]*)/m', $sF3, $matches, PREG_PATTERN_ORDER);
                
                //print_r($matches);
        
                $this->arrRet["M8a"] = $matches[1][0];
                
                // do likelihood ratio tests
                echo("Doing LRT ...\n");
                
                $this->arrRet["M2M1"] = $this->LikelihoodRatio( $this->arrRet["M1"] , $this->arrRet["M2"] , 2);
                $this->arrRet["M8M7"] = $this->LikelihoodRatio( $this->arrRet["M7"] , $this->arrRet["M8"] , 2);
                $this->arrRet["M8M8a"] = $this->LikelihoodRatio( $this->arrRet["M8a"] , $this->arrRet["M8"] , 1);
                
                print_r($this->arrRet);
                
                
                //delete the temp files: 
                foreach(glob($this->sTmpFileBase.'*') as $v){
                        unlink($v);
                }
                
                chdir($this->sTmpFileBase . "../");
                rmdir($this->sTmpFileBase);
                //finished deleting
                
                
                return $this->arrRet;
                
        }
        
        public function LikelihoodRatio($lnNullModel, $lnAlterModel, $df) 
        {
                $nD = ($lnAlterModel - $lnNullModel) * 2 ;
                
                if ($nD < 0) {
                        return 1;
                }
                
                exec(CHI2_PATH . " $df $nD", $arrRet);
                $sRet = implode($arrRet);
                
                preg_match_all('/df\s*\=\s*\S*\s*prob\s*\=\s*([0-9]*\.[0-9]*)/m', $sRet, $matches, PREG_PATTERN_ORDER);
                
                return $matches[1][0];
        }
}

function rrmdir($dir) {
    foreach(glob($dir . '/*') as $file) {
        if(is_dir($file))
            rrmdir($file);
        else
            unlink($file);
    }
    rmdir($dir);
}

function fnExludeDivergentAlignment($arrAln , $bCountNAsDifference = true) {

	$nWindowSize = 21; // must be odd
	$nDivergenceLimit = 0.3;
	$nMinimalIdentical =  $nWindowSize - $nWindowSize * 0.3; // at least this number of characters should be the same
	$nLookEachSide = ($nWindowSize - 1 ) / 2;
	$arrTaxa = array_keys($arrAln);
	$arrSeqs = array_values($arrAln);
	$nTaxa = count($arrTaxa);
	$nLen = strlen($arrSeqs[0]);
	
	$arrRet = array();

	//initialize ret array

	foreach( $arrTaxa as $sTaxa) {
		array_push($arrRet , "");
	}

	// go over each position
	for ($i=0;$i<$nLen;$i++) {

	// first see if they're all equal, if so then don't do complicated things

		$arrCol = array(); // column
		$arrCol[0] = substr($arrSeqs[0], $i, 1);
		$bAllIdentical = true;

		for ($l=1;$l < $nTaxa; $l++ ) {
			$arrCol[$l] = substr($arrSeqs[$l] , $i, 1);
			if ($arrCol[0] != $arrCol[$l]) {
				$bAllIdentical = false;
			}
		}

		$bKeepCol = true; //keep column?
		
		if (!$bAllIdentical) {

			$nWinLeft = ($i - $nLookEachSide < 0)? 0 :  (($i + $nLookEachSide) < $nLen)? $i - $nLookEachSide : ($nLen - 1 - $nLookEachSide*2);

			$sRefSeqInWindow = substr($arrSeqs[0], $nWinLeft, $nWindowSize);
			
			for ($j=1;$j<$nTaxa;$j++) {
				$nMatch = $bCountNAsDifference? similar_text($sRefSeqInWindow , substr($arrSeqs[$j], $nWinLeft, $nWindowSize)) : fnSequenceMatchNAsIdentical($sRefSeqInWindow , substr($arrSeqs[$j], $nWinLeft, $nWindowSize));
				if ( $nMatch < $nMinimalIdentical ) { // doesn't pass filter
					$bKeepCol = false;
					break;
				}
			}

			
		}
		else {
			$bKeepCol = true;
		}

		if ( $bKeepCol ) { // copy column to ret

			for ($k=0; $k < $nTaxa; $k++ ) {
				$arrRet[$k] .= $arrCol[$k];
			}

		}

	}

	return array_combine($arrTaxa , $arrRet);

}

function fnSequenceMatchNAsIdentical($str1, $str2) {

	
	$nLen = strlen($str1);

	$nIdentity = 0;
	for($i=0;$i<$nLen;$i++) {
		$sChar1 = substr($str1 , $i, 1);
		$sChar2 = substr($str2 , $i, 1);
		if ($sChar1 == $sChar2) {
			$nIdentity++;
		}
		else if($sChar1=='N' || $sChar2=='N' ) {
			$nIdentity++;
		}
		else {}
	}

	return $nIdentity;

}

function fnFindLongestReadingFrame($totalLen, $arrStarts, $arrStops) { // the stop codon position is not included

	if (count($arrStops) == 0) { // no stop codon
	
		return array(0, $totalLen-1, $totalLen); // the whole sequence is within a putative ORF
	
	}

	$nBestStart = 0;
	$nBestEnd = 0;
	$nLongestLen = 0;
	
	for ($i=0;$i<count($arrStops);$i++) { // go over each stop codon
	
		$nStopPos = $arrStops[$i];
		$nPreviousStopPos = ($i==0)? -1: $arrStops[$i-1];
		//echo("test$nPreviousStopPos $nStopPos");
		
		$arrStartsInBetween = fnValuesInBetween($arrStarts, $nPreviousStopPos, $nStopPos);

		//print_r($arrStartsInBetween);
		
		$nCurrStart =0;
		$nCurrStop = 0;
		$nCurrLen = 0;
		
		if (count($arrStartsInBetween) ==0) { //if there's no start codon in between the two stops
			if ($nPreviousStopPos == -1) {
				$nCurrStart = 0;
				$nCurrStop = $nStopPos - 1;
				
			}
			else {
				continue; // there is no start codon between two stop codons, can't be ORF
			}
		}
		else {
				$nCurrStart = $arrStartsInBetween[0]; //use the furthest start codon
				$nCurrStop = $nStopPos - 1;
		}
		
		$nCurrLen = $nStopPos - $nCurrStart;
		
		if ($nCurrLen > $nLongestLen) { //update the longest reading frame
			$nBestStart = $nCurrStart;
			$nBestEnd = $nCurrStop;
			$nLongestLen = $nCurrLen;
		}
	
	}
	
	for ($i=0;$i<count($arrStarts);$i++) { // go over each start codon
	
		$nStartPos = $arrStarts[$i];
		$nNextStop = fnFindNextStop($nStartPos, $arrStops);
		
		if ($nNextStop == -1) {
			$nNextStop = $totalLen ;
		}
		
		
		
		$nCurrStart = $nStartPos;
		$nCurrStop = $nNextStop - 1;
		$nCurrLen = $nCurrStop - $nCurrStart + 1;
		
				
		if ($nCurrLen > $nLongestLen) { //update the longest reading frame
			$nBestStart = $nCurrStart;
			$nBestEnd = $nCurrStop;
			$nLongestLen = $nCurrLen;
		}
	
	}
	
	return array($nBestStart, $nBestEnd, $nLongestLen);
	
	
	
}

function fnFindNextStop($nStart, $arrStops) {

	foreach($arrStops as $nStop ) {
		if ($nStop > $nStart) {
			return $nStop;
		}
	}
	
	return -1; //not found
}

function fnValuesInBetween($arr, $nMin, $nMax) {

	$arrRet = array();
	for ($i=0;$i<count($arr);$i++) {
		if ($arr[$i] > $nMin && $arr[$i] < $nMax) {
			$arrRet[] = $arr[$i];
		}
	}
	return $arrRet ;

}

function fnStringPositions($haystack, $needle) {

	$arrRet = array();
	$offset = 0;
	$nPos=0;
	
	while( ($nPos=strpos($haystack , $needle, $offset)) !== false ) 
	{
		$arrRet[] = $nPos;
		$offset = $nPos + 1;
	}
	
	
	return $arrRet;

}

function fnIsStop($sCodon , $sOrientation) {
	$sCodon = strtoupper($sCodon);
	if ($sOrientation == '+')
	{
		return ($sCodon=='TAG' || $sCodon=='TGA' || $sCodon=='TAA');
	}
	else {
		return ($sCodon=='CTA' || $sCodon=='TCA' || $sCodon=='TTA');
		
	}
}

function standard_deviation($aValues, $bSample = false)
{
    $fMean = array_sum($aValues) / count($aValues);
    $fVariance = 0.0;
    foreach ($aValues as $i)
    {
        $fVariance += pow($i - $fMean, 2);
    }
    $fVariance /= ( $bSample ? count($aValues) - 1 : count($aValues) );
    return (float) sqrt($fVariance);
}

function average($aValues)
{
    $fMean = array_sum($aValues) / count($aValues);
    return $fMean;
}

?>
