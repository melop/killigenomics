<?php
ini_set('memory_limit','4048M');
require_once(dirname(__FILE__) . "/config.php");
require_once(dirname(__FILE__) . "/codons.php");

ini_set('display_errors','On');
error_reporting(E_ALL|E_STRICT);

	$Mask_Operation = MASK_BLASTX; // do blastn search
	$Mask_Blast_Db = MASK_DB_NR; // Usse nt database
	$flag_mask = $Mask_Operation | $Mask_Blast_Db ;

	

	/*
	function fnWriteText($file, $content) {
		$hFile = fopen($file.Ext_Writing,"w");
		
		if (!$hFile) {
			return false;
		}
		
		fwrite($hFile, $content);
		
		fclose($hFile);
		
		rename( $file.Ext_Writing, $file);
		
	}
	*/
	
	function fnWriteText($file, $content) {
		$hFile = fopen($file,"w");
		
		if (!$hFile) {
			return false;
		}
		
		fwrite($hFile, $content);
		
		fclose($hFile);
		
		
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
		if (file_exists($file)) {
			return file_get_contents($file);
		}
		else {
			return false;
		}
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
		
		deleteDirectory($this->sTmpFileBase);
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

function fnDNA2Peptide($DNA, $bNoStop = false, $bDiscardN=true) {
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
		$bHasN = (strpos($Codon, "N") !== false);
		
		if (strlen($Codon)<3) {
			continue;
		}
		
		if ($bHasN && $bDiscardN) { // if the codon contains uncertain nucleotide
			
			continue; //throw it away!
		}
		
		//echo($Codon." ");
		if( (!isset($CodonMap[$Codon])) && (!$bHasN))
		{
			if ($bNoStop) {
				continue;
			}
		}
		
		$AminoAcid = "";
		
		if (!$bDiscardN && $bHasN) {
			$AminoAcid="?";
		}
		else {
			$AminoAcid = $CodonMap[$Codon];
		}
		
		if ($bNoStop && $AminoAcid=="*") { // if throw away stop codon
			continue;
		}
		
		$Peptide .= $AminoAcid;
		$DNA2 .= $Codon;
		
	}
	
	return array( $Peptide, $DNA2 ) ;

}

function fnDNA2PeptideWithPosMap($DNA, &$arrMap) {
	global $CodonMap;

	$DNA = strtoupper($DNA);
	$Peptide = "";

	$arrNewPosMap = array(); //key is position of amino acid residue, starting from 1. value is an array of length 3, corresponding to codon 1 2 and 3 positions on the genome.
	$nAminoPos = 0;

	if (strlen($DNA) % 3 !=0 ) {
		echo("DNA not a multiplication of 3!\n");
	}

	
	for ($i=0;$i<=strlen($DNA)-3;$i+=3) {
	
		$Codon = substr($DNA, $i, 3);
		$bHasN = (strpos($Codon, "N") !== false);
		
		if (strlen($Codon)<3) {
			continue;
		}
		
		if ($bHasN || (!isset($CodonMap[$Codon])) ) { // if the codon contains uncertain nucleotide
			
			$AminoAcid = "X";
		} else {
			$AminoAcid = $CodonMap[$Codon];
		}
		
		
		$Peptide .= $AminoAcid;
		$nAminoPos++;
		$nPrevSumPos = ($nAminoPos - 1)*3;
		$arrNewPosMap[$nAminoPos] = array($arrMap[$nPrevSumPos+1] , $arrMap[$nPrevSumPos+2] , $arrMap[$nPrevSumPos+3]);
		
		
	}
	
	return array( $Peptide, $arrNewPosMap ) ;

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

	public function SetBcf($sPath) {
		
		if (! file_exists($sPath)) {
			die("Cannot find bcf file: $sPath !");
		}

		if (! file_exists("$sPath.bci")) {
			echo("Building bcf index for $sPath ...\n");
			exec(BCFTOOLS_PATH." index $sPath" );
			echo("Done.\n");

		}

		$this->sBcfFile = $sPath;
		$this->arrBuffer = array();
		$this->sBufferScfld = '';
		$this->nBufferStart = -1;
		$this->nBufferEnd = -1;
	}

	public function GetVariantsAt($sChr, $nIndex) { //index is 1 based

		//echo("$sChr: $nIndex\n");
	
		if ($sChr != $this->sBufferScfld || $nIndex < $this->nBufferStart || $nIndex > $this->nBufferEnd) { 
			// update buffer
			$this->arrBuffer = null;
			$this->arrBuffer = array();
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
			exec(BCFTOOLS_PATH." view -N $this->sBcfFile ".$sChr.':'.$nIndex.'-'.($nIndex + BCFTOOLS_BUFFERSIZE - 1), $arrRet, $nReturnVal);
		
			if ($nReturnVal != 0 ) {
				echo("Warning: Error occured when trying to pull out buffer region from bcf file with bcftools. ".$this->sBcfFile." : $sChr : $nIndex\n");
				return -1;
			}


			$bBCFOutputVerified = false;
			foreach( $arrRet as $sLn ) {
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
			}
			

			$this->sBufferScfld = $sChr;
			$this->nBufferStart = $nIndex;
			$this->nBufferEnd = $nIndex + BCFTOOLS_BUFFERSIZE - 1;
		
		}

		if (array_key_exists( $nIndex, $this->arrBuffer)) {
			return $this->arrBuffer[$nIndex];
		}
		else {
			return false;
		}
	
	}

}

//this is a new wrapper using stdin and stdout stream, much faster!
class FastaHack{ //fastahack wrapper

	private $sFastaDb;
	private $bContigOK;
	private $sContig;
	private $arrPipes;
	private $hFastaHackProc;
	
	public function SetDb($db) {
	
		if (! file_exists($db)) {
			die("Cannot find fasta database file: $db !");
		}

		
		$this->sFastaDb = $db;

		$this->bContigOK = true;

		$this->sContig = '';

		$oDesc = array(
			   0 => array("pipe", "r"),  // stdin is a pipe that the child will read from
			   1 => array("pipe", "w"),  // stdout is a pipe that the child will write to
			   2 => array("pipe", "w") // stderr is a file to write to
		); //stream descriptor

		$sCMD = FASTAHACK_PATH." -c \"$this->sFastaDb\"";
		$this->hFastaHackProc = proc_open($sCMD , $oDesc, $this->arrPipes);

		if (!is_resource($this->hFastaHackProc))  {
                        die("Cannot query fastahack\n");
                }
	
	}
	
	private function LoadFastaIndex() {
		$hFai = fopen($this->sFastaDb . '.fai' , r);
		
	}

	public function GetCDSRegion($sChr, $nStart, $nEnd) {

		if ($sChr == $this->sContig && $this->bContigOK===false) {
			return false;
		}

		 $this->sContig = $sChr;
		
		if ($nStart > $nEnd) {
			die("Fastahack error: start index > end index");
		}

		if (!is_resource($this->hFastaHackProc))  {
                        die("Cannot query fastahack\n");
                }

		fwrite($this->arrPipes[0] , $sChr.':'.$nStart.'..'.$nEnd."\n");//query
		if (!is_resource($this->hFastaHackProc)) {
			$this->bContigOK = false;
			$this->SetDb($this->sFastaDb); //reload the database;
			return false;
		}

		return trim(fgets($this->arrPipes[1]));

	}

	public function GetContig($sChr) {

		if ($sChr == $this->sContig && $this->bContigOK===false) {
			if (!is_resource($this->hFastaHackProc))  {
				$this->SetDb($this->sFastaDb); //reload the database;
			} else {
				echo("Fastahack error: cannot obtain contig $sChr\n");
				return false;
			}

		}

		 $this->sContig = $sChr;
		

		if (!is_resource($this->hFastaHackProc))  {
                        die("Cannot query fastahack\n");
                }

		fwrite($this->arrPipes[0] , $sChr."\n");//query
		if (!is_resource($this->hFastaHackProc)) {
			$this->bContigOK = false;
			$this->SetDb($this->sFastaDb); //reload the database;
			echo("Fastahack error: cannot obtain contig $sChr\n");
			return false;
		}

		return trim(fgets($this->arrPipes[1]));

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

class XiphoGTF {

	private $sGTF; // gtf file;
	public $arrGenes;
	
	public function LoadGTF($sGTFFile, $sFeatureType="CDS") {
	
		if (! file_exists($sGTFFile)) {
			die("Cannot find GTF file: $sGTFFile");
		}
		
		$this->sGTF = $sGTFFile;
		
		$hFile = fopen($this->sGTF , "r");
		
		if (!$hFile) {
			die("Cannot open GTF file for reading: $sGTFFile");
		}
		
		$this->arrGenes = array();
		$nLnCounter = 0;
		while (($sLn = fgets($hFile)) !== false) {
			$sLn = trim($sLn);
			if ($sLn =="") {
				continue;
			}
			$nLnCounter++;
			$arrFields = array();
        	preg_match('/(\S+)\t(\S+)\t(\S+)\t([0-9]+)\t([0-9]+)\t(\S+)\t(\S)\t([0-2\.])\tgene_id \"(\S+)\".+exon_number \"(\S+)\"/' , $sLn, $arrFields);
        	
        	/*
Array
(
    [0] => JH559423.1	protein_coding	CDS	9637	9748	.	-	0	gene_id "ENSXMAG00000000500"; transcript_id "ENSXMAT00000000496"; exon_number "1"
    [1] => JH559423.1
    [2] => protein_coding
    [3] => CDS
    [4] => 9637
    [5] => 9748
    [6] => .
    [7] => -
    [8] => 0
    [9] => ENSXMAG00000000500
    [10] => 1
)

        	
        	*/
        	//print_r($arrFields);
        	//die();
        	$sGeneName = $arrFields[9];
        	
        	if (!array_key_exists($sGeneName, $this->arrGenes)) {
        		$this->arrGenes[$sGeneName] = array(); //Array for each CDS.
        	}
        	
        	if ($arrFields[3] != $sFeatureType) {
        		continue; // only process CDS
        	}
        	
        	$arrCDS = array( "line" => $nLnCounter, "chr" => $arrFields[1], "start" => $arrFields[4] , "end"=>$arrFields[5], "orientation"=>$arrFields[7], "frame"=>($arrFields[8]=="."? 0:$arrFields[8]) , "exon"=> $arrFields[10]); //chr, start index, end index, orientation, frame
    		array_push($this->arrGenes[$sGeneName] , $arrCDS );
    	
    	}
    	
    	//print_r($this->arrGenes["100_g"]);
		
		fclose($hFile);
	}

}


class NothosGTF {

	private $sGTF; // gtf file;
	public $arrGenes;
	
	public function LoadGTF($sGTFFile, $sFeatureType="CDS") {
	
		if (! file_exists($sGTFFile)) {
			die("Cannot find GTF file: $sGTFFile");
		}
		
		$this->sGTF = $sGTFFile;
		
		$hFile = fopen($this->sGTF , "r");
		
		if (!$hFile) {
			die("Cannot open GTF file for reading: $sGTFFile");
		}
		
		$this->arrGenes = array();
		$nLnCounter = 0;
		while (($sLn = fgets($hFile)) !== false) {
			$sLn = trim($sLn);
			if ($sLn =="") {
				continue;
			}
			$nLnCounter++;
			$arrFields = array();
        	//preg_match('/(\S+)\t(\S+)\t(\S+)\t([0-9]+)\t([0-9]+)\t(\S+)\t(\S)\t([0-2\.])\tID\=\S+;Name=(\S+)\:exon\:([0-9]+);/' , $sLn, $arrFields);
        	preg_match('/(\S+)\t(\S+)\t(\S+)\t([0-9]+)\t([0-9]+)\t(\S+)\t(\S)\t([0-2\.])\tID\=\S+;Name=(\S+)\:cds;Parent/' , $sLn, $arrFields);
        	/*
Array
(
0	=>	GapFilledScaffold_11884	maker	exon	10800	10805	.	+	.	ID=14;Name=TAAR6(6of16)-mRNA-1:exon:0;
1	=>	GapFilledScaffold_11884
2	=>	maker
3	=>	exon
4	=>	10800
5	=>	10805
6	=>	.
7	=>	+
8	=>	.
9	=>	TAAR6(6of16)-mRNA-1
10	=>	0
)

        	
        	*/
        	//print_r($arrFields);
        	//die();

		if (count($arrFields)!=10) {
			continue;
		} 
        	$sGeneName = $arrFields[9];
        	
        	if (!array_key_exists($sGeneName, $this->arrGenes)) {
        		$this->arrGenes[$sGeneName] = array(); //Array for each CDS.
        	}
        	
        	if ($arrFields[3] != $sFeatureType) {
        		continue; // only process CDS
        	}
        	
        	$arrCDS = array( "line" => $nLnCounter, "chr" => $arrFields[1], "start" => $arrFields[4] , "end"=>$arrFields[5], "orientation"=>$arrFields[7], "frame"=> 0 /*($arrFields[8]=="."? 0:$arrFields[8])*/ , "exon"=> count($this->arrGenes[$sGeneName]) ); //chr, start index, end index, orientation, frame  (count($this->arrGenes[$sGeneName])+1)  $arrFields[10]
    		array_push($this->arrGenes[$sGeneName] , $arrCDS );
    	
    	}
    	
    	//print_r($this->arrGenes["100_g"]);
		
		fclose($hFile);
	}

}

class MappedCDSGFF { //parse the CDS annotations in GFF and provide easy ways to access positional information
	private $sGFF; // gtf file;
	private $sRefGenome;
	private $oRefGenome;
	private $arrRefGenomeObjs;
	private $arrRefGenomeLoaded;
	private $arrGenes;
	private $arrRNA2Gene; //mRNA id to parent gene id map.
	private $arrSpGeneId2RNAID;
	private $bFastaHackLoaded;
	private $arrGeneCoordMap; //key is scaffold name => array( start1 => array(end, geneid), start2 => array(end, geneid) ) 


	//this assumes a gene -> mRNA -> exon (CDS) structure.
	
	public function LoadGFF($sGFFFile, $sGenome, $sTreatTypeAsCDS = "CDS") { //treat this type as "cds".
		$bRet = true;
	
		if (! file_exists($sGFFFile)) {
			die("Cannot find GFF file: $sGFFFile");
		}
		
		$this->sGFF = $sGFFFile;

		if (is_array($sGenome) ) {
			$this->sRefGenome = $sGenome[0]; //if is array, set the first element
			foreach($sGenome as $sG) {
				$this->arrRefGenomeObjs[$sG] = new FastaHack();
				$this->arrRefGenomeLoaded[$sG] = false;
			}
		} else {
			$this->arrRefGenomeObjs = false;
			$this->sRefGenome = $sGenome;
			$this->oRefGenome = new FastaHack(); //load
			//$this->oRefGenome->SetDb($this->sRefGenome);	
			$this->bFastaHackLoaded = false;
		}


		$this->arrGeneCoordMap = array("+" => array(), "-" => array());
		
		$hFile = fopen($this->sGFF , "r");
		
		if (!$hFile) {
			die("Cannot open GFF file for reading: $sGFFFile");
		}
		
		$this->arrGenes = array();

		while (($sLn = fgets($hFile)) !== false) {
			$sLn = trim($sLn);
			if ($sLn =="") continue;
			if ($sLn[0] =='#') continue;
			$arrF = explode("\t", $sLn);

			if (count($arrF) != 9) continue;

			if ($arrF[3] > $arrF[4]) {
				die("Error: In Gff file, start cannot be larger than end coordinate!\n");
			}

			$sGenePredictor = $arrF[1];
			$sFeatureType = $arrF[2];

			if (!array_key_exists($sGenePredictor , $this->arrGenes)) {
				$this->arrGenes[$sGenePredictor] = array();
			}

			$arrAnnot = $this->fnParseAnnotation($arrF[8]);

			if ($sFeatureType == 'gene') { // check if gene is there
				$this->fnAssertAnnotField($arrAnnot , array('ID') );
				if (!array_key_exists($arrAnnot['ID'] , $this->arrGenes[$sGenePredictor] )) {
					$this->arrGenes[$sGenePredictor][$arrAnnot['ID']] = array('scf'=>$arrF[0] , 
												'start' =>$arrF[3],
												 'end'=> $arrF[4], 
												'strand' => $arrF[6], 
												'annot' => $arrAnnot, 
												'mRNAs' => array());
					if (!array_key_exists( $arrF[0] , $this->arrGeneCoordMap[ $arrF[6]]) ) {
						$this->arrGeneCoordMap[$arrF[6]][$arrF[0]] = array();
					}

					$this->arrGeneCoordMap[$arrF[6]][$arrF[0]][$arrF[3]] = array($arrF[4] , $arrAnnot['ID']);

				} else {
					die("Error: Gene id ".$arrAnnot['ID']." was already defined. gene id must be unique; \n");
				}
				continue;
			}

			if ($sFeatureType == 'mRNA') { // check if gene is there
				$this->fnAssertAnnotField($arrAnnot , array('ID') );
				if (!array_key_exists('Parent' , $arrAnnot) ) {
					//If there is no definition of parent
					$arrAnnot['Parent'] = "gene.".$arrAnnot['ID']; //make it up
					$this->arrGenes[$sGenePredictor][$arrAnnot['Parent']] = array('scf'=>$arrF[0] , 
												'start' =>$arrF[3],
												 'end'=> $arrF[4],
												 'strand' => $arrF[6], 
												'annot' => array('ID'=>$arrAnnot['Parent']),
												'mRNAs' => array());
				}
				
				if (!array_key_exists($arrAnnot['Parent'] , $this->arrGenes[$sGenePredictor] )) { 
					echo("Warning: Gene not found: ".$arrAnnot['Parent'].". Ignored\n");
					$bRet = false;
					continue;
				} 

				$this->arrGenes[$sGenePredictor][$arrAnnot['Parent']]['mRNAs'][$arrAnnot['ID']] = array('scf'=>$arrF[0] , 
												'start' =>$arrF[3],
												 'end'=> $arrF[4],
												 'strand' => $arrF[6], 
												'annot' => $arrAnnot,
												'CDSs' => array(), 'CDSsorted' => false);

				$this->arrRNA2Gene[$arrAnnot['ID']] = $arrAnnot['Parent'];

				if (array_key_exists('spgeneid', $arrAnnot)) {
					$this->arrSpGeneId2RNAID[$arrAnnot['spgeneid']] = $arrAnnot['ID'];
				}

				continue;
			}

			if ($sFeatureType == $sTreatTypeAsCDS) { // check if gene is there
				$this->fnAssertAnnotField($arrAnnot , array('Parent') );
				if (!array_key_exists($arrAnnot['Parent'] , $this->arrRNA2Gene)) {
					echo("The CDS's parent gene ".$arrAnnot['Parent']." is undefined, ignored\n");
					$bRet = false;
					continue;
				}

				$sGeneID = $this->arrRNA2Gene[$arrAnnot['Parent'] ];
				if ($this->arrGenes[$sGenePredictor][$sGeneID]['mRNAs'][$arrAnnot['Parent']]['strand'] != $arrF[6]) {
					print_r($this->arrGenes);
					print_r($this->arrRNA2Gene);
					print_r($arrF);
					echo("Warning:In mRNA ".$arrAnnot['Parent']." of gene $sGeneID, the CDS strand is different from the mRNA! This gene will be removed.\n" );
					//remove parent mRNA
					unset($this->arrRNA2Gene[$arrAnnot['Parent'] ]);
					unset($this->arrGenes[$sGenePredictor][$sGeneID]);
					$bRet = false;
					continue;
				}
				$this->arrGenes[$sGenePredictor][$sGeneID]['mRNAs'][$arrAnnot['Parent']]['CDSs'][$arrF[3]] = $arrF[4];

				continue;
			}

		}

		foreach(array('+','-') as $sOrient) {
			$arrScf = array_keys($this->arrGeneCoordMap[$sOrient]);
			foreach($arrScf as $sScf) {
				ksort($this->arrGeneCoordMap[$sOrient][$sScf] , SORT_NUMERIC);
			}
		}
		//print_r($this->arrGenes);
		fclose($hFile);

		//print_r($this->arrGeneCoordMap);
		//print_r($this->GetNeighborGeneOnScf('PLPv12scf130', 484631 , '+', 1));
		//print_r($this->GetNeighborGeneOnScf('PLPv12scf130', 484631 , '+', -1));
		//print_r($this->GetNeighborGeneOnScf('PLPv12scf130', 484631 , '+', 100));
		//die();
		return $bRet;
	}

	public function GetNeighborGeneOnScf($sScf, $nStart, $sOrient, $nOffSet = 1) { //give the start position and the scaffold of the gene, returns the id of the NEXT gene
		if (!array_key_exists($sScf, $this->arrGeneCoordMap[$sOrient]) ) {
			return false;
		}


		if (!array_key_exists($nStart, $this->arrGeneCoordMap[$sOrient][$sScf])) {
			return false;
		}

		$arrKeys = array_keys($this->arrGeneCoordMap[$sOrient][$sScf]);
		$arrKeysInv = array_flip($arrKeys );
		$nIdx = $arrKeysInv[$nStart];
		$nNextIdx = $nIdx +  $nOffSet;

		if ( ($nNextIdx < count($arrKeys)) && ($nNextIdx >=0) ) {
			$nNextGeneStart = $arrKeys[$nNextIdx];
			return array( 'geneid' => $this->arrGeneCoordMap[$sOrient][$sScf][$nNextGeneStart][1], 'start'=>$nNextGeneStart , 'end' => $this->arrGeneCoordMap[$sOrient][$sScf][$nNextGeneStart][0] ); // return the gene id
		}

		return false;
	}

	public function &GetGenes($sGenePredictor) {
		if (!array_key_exists($sGenePredictor, $this->arrGenes )) {
			//die("Gene predictor not found in GFF: $sGenePredictor\n");
			return array(); //empty array
		} 

		return $this->arrGenes[$sGenePredictor];
	}

	public function GenePredictors() { // returns the gene predictors
		return array_keys($this->arrGenes);
	}

	public function GetmRNA($sGenePredictor, $smRNAID) { // 
		if (!array_key_exists($smRNAID , $this->arrRNA2Gene) ) {
			die("mRNA ID not found: $smRNAID\n");
		}
		$oRNA = ($this->arrGenes[$sGenePredictor][$this->arrRNA2Gene[$smRNAID]]['mRNAs'][$smRNAID]);
		if (!$oRNA['CDSsorted']) {
			if ($oRNA['strand'] == '+') { // sort CDS ascending
				ksort($oRNA['CDSs'], SORT_NUMERIC);
			} else {
				ksort($oRNA['CDSs'], SORT_NUMERIC); // regardless of strand ksort ascending!
			}

			$oRNA['CDSsorted'] = true;
		}

		return $oRNA;
	}

	public function mRNAIDExistsInPredictor($sGenePredictor, $smRNAID) { // 

		if (!array_key_exists($sGenePredictor , $this->arrGenes) ) {
			return false;
		}

		if (!array_key_exists($smRNAID , $this->arrRNA2Gene) ) {
			die("mRNA ID not found: $smRNAID\n");
		}
		$sGeneID = $this->arrRNA2Gene[$smRNAID];
		return array_key_exists( $sGeneID , $this->arrGenes[$sGenePredictor]);		
	}

	public function ExtractmRNASequence($sGenePredictor, $smRNAID, $sGenome = "") { // 
		if (!array_key_exists($smRNAID , $this->arrRNA2Gene) ) {
			die("mRNA ID not found: $smRNAID\n");
		}
		$oRNA = ($this->arrGenes[$sGenePredictor][$this->arrRNA2Gene[$smRNAID]]['mRNAs'][$smRNAID]);

		return $this->ExtractmRNASequenceByModel($oRNA,  $sGenome);
	}

	public function ExtractmRNASequenceByModel(&$oRNA , $sGenome = "") { // 

		if ($this->arrRefGenomeObjs === false) { //single genome mode
			echo("single genome mode;\n");
			if (! $this->bFastaHackLoaded) {
				$this->oRefGenome->SetDb($this->sRefGenome);	
				$this->bFastaHackLoaded = true;
			}
		} else {
			if (!array_key_exists($sGenome , $this->arrRefGenomeObjs) ) {
				die("Multi-reference mode error: $sGenome doesn't exist\n");
			}
			if (! $this->arrRefGenomeLoaded[$sGenome ]) {
				$this->arrRefGenomeObjs[$sGenome]->SetDb($sGenome);
				$this->arrRefGenomeLoaded[$sGenome ] = true;
			}

			$this->sRefGenome = $sGenome;
			$this->bFastaHackLoaded = true;
			$this->oRefGenome = $this->arrRefGenomeObjs[$sGenome];

		}

		if (!$oRNA['CDSsorted']) {
			if ($oRNA['strand'] == '+') { // sort CDS ascending
				ksort($oRNA['CDSs'], SORT_NUMERIC);
			} else {
				ksort($oRNA['CDSs'], SORT_NUMERIC); // regardless of strand ksort ascending!
			}

			$oRNA['CDSsorted'] = true;
		}

		//now go over all CDS
		$sDNA = "";
		$arrPosMap = array(); //key is position on extracted DNA, value is position on the genome. both 1-based.
		foreach($oRNA['CDSs'] as $nStart => $nEnd) {
			$sPartDNA = $this->oRefGenome->GetCDSRegion($oRNA['scf'], $nStart, $nEnd);
			$nCurrLen = strlen($sDNA);
			$nPartDNA = strlen($sPartDNA);
			if ($nPartDNA != (abs( $nEnd - $nStart) + 1) ) {
				die("Something went wrong when extracting sequencing in range: ".$oRNA['scf'].":$nStart..$nEnd ; len: $nPartDNA \n$sPartDNA");
			}

			$nPosEnd = $nCurrLen+$nPartDNA;
			for($nPos = $nCurrLen + 1;$nPos<=$nPosEnd;$nPos++) {
				$arrPosMap[$nPos] = $nStart + ($nPos - $nCurrLen - 1);
			}
			$sDNA .= $sPartDNA;
		}

		//assertion:
		if (count($arrPosMap) != strlen($sDNA) ) {
			die("Something is wrong when computing position map.\n");
		}

		if ($oRNA['strand'] == '-') {
			//print_r($arrPosMap);
			//if on reverse strand, reverse complement DNA sequence
			$sDNA = fnReverseComplement($sDNA);
			//reverse positional map:
			$arrRevValues = array_values($arrPosMap);
			rsort($arrRevValues);
			$arrPosMap = array_combine(array_keys($arrPosMap ) ,  $arrRevValues);

			//echo($sDNA."\n");
			//print_r($arrPosMap);
		}

		$arrTranslation = fnDNA2PeptideWithPosMap($sDNA, $arrPosMap);


		return array('mRNA' => $sDNA, 'RNA2GenomeMap' => $arrPosMap , 'AA' => $arrTranslation[0] , 'Peptide2GenomeMap' => $arrTranslation[1] );
	}

	public function fnParseAnnotation($s) {
		$arrMap = array();
		$arrF = explode(';' , $s);
		foreach($arrF as $v) {
			$arrPair = explode('=' , $v);
			if (count($arrPair) !=2) continue;
			$sKey = trim($arrPair[0]);
			$sValue = preg_replace('/^"|"$/', '', trim($arrPair[1]));
			$arrMap[$sKey] = $sValue;
		}
		return $arrMap;
	}

	public function fnSpGeneId2RNAID($s) {
		if ( !array_key_exists( $s, $this->arrSpGeneId2RNAID)) return false;
		return $this->arrSpGeneId2RNAID[$s];
	}
	private function fnAssertAnnotField(&$arrAnnot, $arrAssert) {
		foreach($arrAssert as &$sAss) {
			if (!array_key_exists($sAss , $arrAnnot) ) {
				die("Expected to find annotation in GFF $sAss, but failed\nDie\n");
			}
		}
	}
}

function fnTranslatorX($arrSeqs , $nPMinSeqConsv=0.75, $nPMinSeqFlank = 0.85, $nMinContNonConsv=2, $nMinBlockLen=5) {
	$sTmp = TRANSLATORX_TEMP_PATH . "/".rand( 1000000 , 9999999) ;
	exec("mkdir -p $sTmp");
	$hIn = fopen($sTmp."/in.fa", "w");
	foreach($arrSeqs as $sName => $sSeq) {
		fwrite($hIn, ">".$sName."\n".$sSeq."\n");
	}

	$nSEQs = count($arrSeqs);
	$nMinSeqConsv=intval($nSEQs * $nPMinSeqConsv + 1);
	$nMinSeqFlank=$nSEQs * $nPMinSeqFlank;


	exec("cd $sTmp; perl ".TRANSLATORX_PATH." -i in.fa -o tranX -p F -c 1 -g \"-b1=$nMinSeqConsv -b2=$nMinSeqFlank -b3=$nMinContNonConsv -b4=$nMinBlockLen -b5=h\" > tr.log 2>&1; cd ../../");

	//print_r(glob($sTmp."/*"));
	if ( (!file_exists($sTmp."/tranX.nt_ali.fasta")) || (!file_exists($sTmp."/tranX.nt_cleanali.fasta")) || (!file_exists($sTmp."/tranX.aa_ali.fasta")) || (!file_exists($sTmp."/tranX.aa_cleanali.fasta"))) {
		return false;
	}

	$sAln = file_get_contents($sTmp."/tranX.nt_ali.fasta");
	$sCleanedAln = file_get_contents($sTmp."/tranX.nt_cleanali.fasta");
	$sAlnAA = file_get_contents($sTmp."/tranX.aa_ali.fasta");
	$sCleanedAA = file_get_contents($sTmp."/tranX.aa_cleanali.fasta");

	exec("rm -R $sTmp");
	return(array($sAln , $sCleanedAln, $sAlnAA, $sCleanedAA ));
}


class BlastTextRet {

	private $sBlastRet; // gtf file;
	public $arrQueries;
	
	public function LoadText($sBlastRetFile) {
	
		if (! file_exists($sBlastRetFile)) {
			die("Cannot find Blast result file: $sBlastRetFile");
		}
		
		$this->sBlastRet = $sBlastRetFile;
		
		$hFile = fopen($this->sBlastRet , "r");
		
		if (!$hFile) {
			die("Cannot open Blast result file for reading: $sBlastRetFile");
		}
		
		$this->arrQueries = array();
		
		while (($sLn = fgets($hFile)) !== false) {
			$sLn = trim($sLn);
			if ($sLn == "") {
				continue;
			}
			
			$arrFields = explode("\t", $sLn); 
        	$sQueryName = $arrFields[0];
        	
        	if (!array_key_exists($sQueryName, $this->arrQueries)) {
        		$this->arrQueries[$sQueryName] = array(); //Array for each CDS.
        	}
        	
/*
Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score
*/
        	$arrCDS = array( "subject" => $arrFields[1], "identity" => $arrFields[2] , "aln_len"=>$arrFields[3], "mismatches"=>$arrFields[4], "gaps"=>$arrFields[5] , "qstart"=> $arrFields[6], "qend"=> $arrFields[7], "sstart"=> $arrFields[8], "send"=> $arrFields[9], "e"=> $arrFields[10], "bitscore"=> $arrFields[11]); //chr, start index, end index, orientation, frame
    		array_push($this->arrQueries[$sQueryName] , $arrCDS );
    	
    	}
    	
    	//print_r($this->arrQueries["100_g"]);
		
		fclose($hFile);
	}

}


class AUTest {
	private $sAlignment; //must be in phylip format
	private $arrTrees; //an array of alternative topologies in newick format
	private $nInstanceID ;
	private $sTempDir;
        
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
		exec("mkdir -p " . AU_TEMP_PATH);

		$this->sTempDir = realpath( AU_TEMP_PATH  ) ."/" . $this->nInstanceID . "/";
		mkdir($this->sTempDir);
		chdir($this->sTempDir);

		$hAln = fopen($this->sTempDir."/".AUTest::IN_FILE_NAME, "w");
		fwrite($hAln , $this->sAlignment);
		fclose($hAln);

		$hTrees = fopen($this->sTempDir."/".AUTest::TREE_FILE_NAME, "w");

		foreach( $this->arrTrees as $sTree) {
			fwrite($hTrees, $sTree.PHP_EOL);
		}
		fclose($hTrees);

		//Run RAxML

		$sCmd = RAXML_PATH . " -f g -s " . AUTest::IN_FILE_NAME . " -m GTRGAMMA -z ". AUTest::TREE_FILE_NAME ." -n ".AUTest::RAXML_OUT;
		echo($sCmd."\n");
		exec( $sCmd );

		$sRaxmlLnOut = $this->sTempDir."/".AUTest::RAXML_RET_PREFIX.".".AUTest::RAXML_OUT;
		if (!file_exists($sRaxmlLnOut )) {
			echo("RAxML did not generate output $sRaxmlLnOut! ");
			chdir(dirname(__FILE__) );
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

		//rrmdir($this->sTempDir);//remove dir

		return $arrRet;

	}

	public function GetTempDir() {
			return $this->sTempDir;
	}

	public function RemoveTemp() {
		if (file_exists($this->sTempDir)) {
			rrmdir($this->sTempDir);
			return $this->sTempDir;
		}
	}
}

class RAXML {
	private $sAlignment; //must be in phylip format
	private $arrTrees; //an array of alternative topologies in newick format
	private $nInstanceID ;
	private $sConstraintTree;
	private $sTempDir;
        
        const IN_FILE_NAME = "in.phy";
	const RAXML_OUT = "output";
        const TREE_FILE_NAME = "in.tre";
	const RAXML_RET_PREFIX = "RAxML_bestTree";


	public function SetAlignment($sAlg) {
		$this->sAlignment = $sAlg;
		$this->sConstraintTree = '';
	}

	public function SetConstraint($sTree) {
		$this->sConstraintTree = $sTree;
	}


	public function GetBestTree($bKeepTmpFolder = true) {

 		$this->nInstanceID = rand( 1000000 , 9999999) ;


		$sTreeRet = "";
		try {

			if (!file_exists(AU_TEMP_PATH)) {
				mkdir(AU_TEMP_PATH);
			}
		}
		catch(Exception  $e) 
		{}

		$this->sTempDir = realpath(AU_TEMP_PATH ). "/". $this->nInstanceID . "/";
		mkdir($this->sTempDir);
		chdir($this->sTempDir);

		$hAln = fopen(RAXML::IN_FILE_NAME, "w");
		fwrite($hAln , $this->sAlignment);
		fclose($hAln);

		$sConstraintParam = '';
		if ($this->sConstraintTree != '') {
			$hConstraint = fopen(RAXML::IN_FILE_NAME . ".constraint.tre", "w");
			fwrite($hConstraint , $this->sConstraintTree);
			fclose($hConstraint);
			$sConstraintParam = " -g \"". RAXML::IN_FILE_NAME . ".constraint.tre" ."\" ";
		}

		//Run RAxML
		//raxmlHPC -s in.phy -m GTRGAMMA -n out -f a -x 23 -N 100

		$sCmd = RAXML_PATH . " $sConstraintParam -T 4 -f a -s \"". RAXML::IN_FILE_NAME . "\" -p 23333 -m GTRGAMMA -x 23 -N 10 -n ".RAXML::RAXML_OUT;
		echo("$sCmd \n");
		exec( $sCmd );

		$sRaxmlLnOut = RAXML::RAXML_RET_PREFIX.".".RAXML::RAXML_OUT;
		if (!file_exists($sRaxmlLnOut )) {
			echo("RAxML did not generate output $sRaxmlLnOut! ");
			chdir(dirname(__FILE__) );
			return array("Success"=>false, "Tree"=>"" , "tmpfolder" => $this->sTempDir);
		}
		else {
			$sTreeRet = fnReadText($sRaxmlLnOut);
		}
		
		
		chdir(dirname(__FILE__) );

		if (!$bKeepTmpFolder) {
			rrmdir($this->sTempDir);//remove dir
			}

		return array("Success"=>true, "Tree"=>$sTreeRet , "tmpfolder" => $this->sTempDir);

	}

	public function GetTempDir() {
			return $this->sTempDir;
	}

	public function RemoveTemp() {
		if (file_exists($this->sTempDir)) {
			rrmdir($this->sTempDir);
			return $this->sTempDir;
		}
	}

}

function fnAln2Phylip($arrAln) {
	$arrKeys = array_keys($arrAln);
	$s = " ".count($arrAln)." ";
	$s .= strlen($arrAln[$arrKeys[0]])."\n";
	foreach($arrAln as $sTaxon => $sSeq) {
		$s .= $sTaxon."\t".$sSeq."\n";
	}
	return $s;
}

class CodemlUtil {

        private $nInstanceID ;
        private $sTmpFileBase;
        private $sTreeContent;
        private $arrLabels = array();
        private $arrSeqs = array();
        private $arrRet = array();
	private $sCurrDir = "?";
        
        private $bInit = false;
        
        const IN_FILE_NAME = "in.fas";
        const TREE_FILE_NAME = "in.tre";
        const PREFIX_CTRL_FILE1 = "codeml.1";
        const PREFIX_CTRL_FILE2 = "codeml.2";
        const PREFIX_CTRL_FILE3 = "codeml.3";
        
        
        private function initialize() {
        
          if ($this->bInit) {return;};
          
          $this->sCurrDir = ($this->sCurrDir=="?")? getcwd() ."/":$this->sCurrDir;
	  do {
      	  $this->nInstanceID = rand( 10000000 , 99999999) ;
          $this->sTmpFileBase = CODEML_TMP."/".$this->nInstanceID."/";
	} while (file_exists($this->sTmpFileBase));

	  chdir($this->sCurrDir);
          mkdir($this->sTmpFileBase, 0777, true); // create temp file folder
          
          $this->sTmpFileBase = realpath(CODEML_TMP."/".$this->nInstanceID) ."/";
          echo("WorkingPath: ".$this->sCurrDir . PHP_EOL);
          echo("TmpPath: ".$this->sTmpFileBase . PHP_EOL);

          //Read the template control file:
          $sCtrlTemplate = fnReadText(CODEML_CTRL_TEMPLATE);
          
          
          //Common variables shared by all control files:
          $sCtrlTemplate = str_replace("%SEQ_FILE%",  CodemlUtil::IN_FILE_NAME ,  $sCtrlTemplate);
          $sCtrlTemplate = str_replace("%TREE_FILE%",   CodemlUtil::TREE_FILE_NAME  ,  $sCtrlTemplate);
          
          //Control file 1 for just getting plain dN/dS:
          $sCtrl1 = str_replace("%OUT_FILE%",   CodemlUtil::PREFIX_CTRL_FILE1.".out"  ,  $sCtrlTemplate);
          $sCtrl1 = str_replace("%NSSITES%", "0"  ,  $sCtrl1);
          $sCtrl1 = str_replace("%RUN_MODE%", "-2"  ,  $sCtrl1); //Run mode pairwise
          $sCtrl1 = str_replace("%FIX_OMEGA%", "0"  ,  $sCtrl1); //fix omega
          $sCtrl1 = str_replace("%FIX_ALPHA%", "1"  ,  $sCtrl1); //fix gamma
          
          fnWriteText($this->sTmpFileBase .  CodemlUtil::PREFIX_CTRL_FILE1 . ".ctl" , $sCtrl1);
          
          //Control file 2 for M1, M2, M7, M8
          $sCtrl2 = str_replace("%OUT_FILE%",   CodemlUtil::PREFIX_CTRL_FILE2.".out"  ,  $sCtrlTemplate);
          $sCtrl2 = str_replace("%NSSITES%", "1 2 7 8"  ,  $sCtrl2);
          $sCtrl2 = str_replace("%RUN_MODE%", "0"  ,  $sCtrl2); //Run mode user tree
          $sCtrl2 = str_replace("%FIX_OMEGA%", "0"  ,  $sCtrl2); //estimate omega
          $sCtrl2 = str_replace("%FIX_ALPHA%", "1"  ,  $sCtrl2); //fix gamma
          
          fnWriteText($this->sTmpFileBase .  CodemlUtil::PREFIX_CTRL_FILE2 . ".ctl", $sCtrl2);
          
          //Control file 3 for M8a
          $sCtrl3 = str_replace("%OUT_FILE%",   CodemlUtil::PREFIX_CTRL_FILE3.".out"  ,  $sCtrlTemplate);
          $sCtrl3 = str_replace("%NSSITES%", "8"  ,  $sCtrl3);
          $sCtrl3 = str_replace("%RUN_MODE%", "0"  ,  $sCtrl3); //Run mode user tree
          $sCtrl3 = str_replace("%FIX_OMEGA%", "1"  ,  $sCtrl3); //estimate omega
          $sCtrl3 = str_replace("%FIX_ALPHA%", "1"  ,  $sCtrl3); //fix gamma
          
          fnWriteText($this->sTmpFileBase .  CodemlUtil::PREFIX_CTRL_FILE3 . ".ctl", $sCtrl3);
          
          //chdir($this->sTmpFileBase ); //change working dir
          
          $this->bInit = true;
          
          
        }

private function initializeBranchSite() {
        
          if ($this->bInit) {return;};
          
          $this->sCurrDir = ($this->sCurrDir=="?")? getcwd() ."/":$this->sCurrDir;
      $this->nInstanceID = rand( 1000000 , 9999999) ;
          $this->sTmpFileBase = CODEML_TMP."/".$this->nInstanceID."/";
	  chdir($this->sCurrDir);
          mkdir($this->sTmpFileBase, 0777, true); // create temp file folder
          
          $this->sTmpFileBase = realpath(CODEML_TMP."/".$this->nInstanceID) ."/";
          echo("WorkingPath: ".$this->sCurrDir . PHP_EOL);
          echo("TmpPath: ".$this->sTmpFileBase . PHP_EOL);

          //Read the template control file:
          $sCtrlTemplate = fnReadText(CODEML_CTRL_TEMPLATE);
          
          
          //Common variables shared by all control files:
          $sCtrlTemplate = str_replace("%SEQ_FILE%",  CodemlUtil::IN_FILE_NAME ,  $sCtrlTemplate);
          $sCtrlTemplate = str_replace("%TREE_FILE%",   CodemlUtil::TREE_FILE_NAME  ,  $sCtrlTemplate);
	  $sCtrlTemplate = str_replace("%MODEL%", "2"  ,  $sCtrlTemplate);
          $sCtrlTemplate = str_replace("%NSSITES%", "2"  ,  $sCtrlTemplate);
          $sCtrlTemplate = str_replace("%FIX_ALPHA%", "1"  ,  $sCtrlTemplate); //fix gamma
          $sCtrlTemplate = str_replace("%RUN_MODE%", "0"  ,  $sCtrlTemplate); //Run mode user tree

          //Control file 2 for Model A null
          $sCtrl2 = str_replace("%OUT_FILE%",   CodemlUtil::PREFIX_CTRL_FILE2.".out"  ,  $sCtrlTemplate);
          $sCtrl2 = str_replace("%FIX_OMEGA%", "1"  ,  $sCtrl2); //estimate omega
          $sCtrl2 = str_replace("%INIT_OMEGA%", "1"  ,  $sCtrl2); //estimate omega
          
          fnWriteText($this->sTmpFileBase .  CodemlUtil::PREFIX_CTRL_FILE2 . ".ctl", $sCtrl2);
          
          //Control file 3 for Model A omega > 1
          $sCtrl3 = str_replace("%OUT_FILE%",   CodemlUtil::PREFIX_CTRL_FILE3.".out"  ,  $sCtrlTemplate);
          $sCtrl3 = str_replace("%FIX_OMEGA%", "0"  ,  $sCtrl3); //estimate omega
          $sCtrl3 = str_replace("%INIT_OMEGA%", "1.5"  ,  $sCtrl3); //estimate omega
          fnWriteText($this->sTmpFileBase .  CodemlUtil::PREFIX_CTRL_FILE3 . ".ctl", $sCtrl3);
          
          //chdir($this->sTmpFileBase ); //change working dir
          
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
        
        public function RunCodeml($bOnlyGetdNdS, $bKeepTmpFolder=false) {
        
                if (count($this->arrRet) > 0)
                {
                        return $this->arrRet;
                }
                
                $this->initialize();
                // Prepare the in files:
                $sFas = "";
		chdir($this->sCurrDir); //always make sure to start from the working directory.
				 
                for ($i=0;$i<count($this->arrSeqs);$i++) {
                        
                                $sFas .= '>'.trim($this->arrLabels[$i])."\n";
                                $sFas .= trim($this->arrSeqs[$i])."\n";
                                                        
                }

                
                fnWriteText( $this->sTmpFileBase . CodemlUtil::IN_FILE_NAME, $sFas);
				$this->arrRet["tmpfolder"] = $this->sTmpFileBase;
                $this->arrRet["Success"] = true;
				
                if (!isset($this->sTreeContent) ) {
                        $sTree = "(";
                        for ($i=0;$i<count($this->arrSeqs);$i++) {
                                                       
                                $sTree .= ($i==0)? $this->arrLabels[$i]: ", ". $this->arrLabels[$i];
                                $sTree .= ($i == (count($this->arrSeqs)-1) )? ");" : "";
                        
                        }
                        fnWriteText( $this->sTmpFileBase . CodemlUtil::TREE_FILE_NAME, $sTree );
                }
                else {
                	
                        $sTree = $this->sTreeContent;
                        if (substr(trim($sTree), -1, 1) != ';') {
                        	    $sTree .= ';';             	
                        }
                        fnWriteText( $this->sTmpFileBase . CodemlUtil::TREE_FILE_NAME, $sTree );
                }
                
                
                
                //Run Ctrl file 1: ========================================
                chdir($this->sTmpFileBase ); //go to temp directory
                echo("Running Model 0 ...\n");

                $sCMD = CODEML_PATH . " ".  CodemlUtil::PREFIX_CTRL_FILE1 . ".ctl";
                //echo($sCMD);
                
                //return;
                
                $sRet = exec($sCMD);
                $sF1 = fnReadText(  CodemlUtil::PREFIX_CTRL_FILE1 . ".out");
                
                //Parse the file:
                /* $this->arrRet[]; */
                //echo($sF1);
                // works for paml4.4 preg_match_all('/t\=\s*\S+\s+S\=\s*\S+\s*N\=\s*\S+\s*dN\/dS\=\s*([0-9]*\.[0-9]+)\s*dN\=\s*([0-9]*\.[0-9]+)\s*dS\=\s*(([0-9]*\.[0-9]+))/m', $sF1, $matches, PREG_PATTERN_ORDER);
                if ($sF1===false) {
					$this->arrRet["S"] = "";
					$this->arrRet["N"] = "";
					$this->arrRet["dN"] = "";
					$this->arrRet["dS"] = "";
					$this->arrRet["dNdS"] = "";
					//chdir($this->sTmpFileBase . "../");
					 $this->arrRet["Success"] = false;
					return $this->arrRet;
				}
				
				preg_match_all('/t\=\s*\S+\s+S\s*\=\s*(\S+)\s*N\s*\=\s*(\S+)\s*dN\/dS\=\s*([0-9]*\.[0-9]+)\s*dN\s*\=\s*([0-9]*\.[0-9]+)\s*dS\s*\=\s*([0-9]*\.[0-9]+)/m', $sF1, $matches, PREG_PATTERN_ORDER);
                
				
				//t=50.0000  S=    46.9  N=   178.1  dN/dS= 0.0086  dN= 0.6661  dS=77.3772

                //print_r($matches);
                $this->arrRet["S"] = $matches[1]; // an array with the pairwise dN dS etc.
                $this->arrRet["N"] = $matches[2];
                $this->arrRet["dN"] = $matches[4];
                $this->arrRet["dS"] = $matches[5];
                $this->arrRet["dNdS"] = $matches[3];
                
                if ($bOnlyGetdNdS) {
						chdir($this->sCurrDir); //always make sure to start from the working directory.
						if (!$bKeepTmpFolder) {
							deleteDirectory($this->sTmpFileBase);
							}
                        return $this->arrRet;
                }
                
                chdir($this->sTmpFileBase ); //go to temp directory
                //Run Ctrl file 2 ======================================
                echo("Running Model 1 2 7 8 ...\n");
                
                $sCMD = CODEML_PATH . " ".  CodemlUtil::PREFIX_CTRL_FILE2 . ".ctl";
                //echo("$sCMD\n");
                $sRet = exec($sCMD);
                $sF2 = fnReadText(  CodemlUtil::PREFIX_CTRL_FILE2 . ".out");
                
                //Parse the file:
                /* $this->arrRet[]; */
                //echo($sF2);
                preg_match_all('/lnL\(ntime\:\s*\S*\s*np\:\s*\S*\)\:\s+(\-[0-9]*\.[0-9]*)/m', $sF2, $matches, PREG_PATTERN_ORDER);
                
                //print_r($matches);
        
                $this->arrRet["M1"] = $matches[1][0];
                $this->arrRet["M2"] = $matches[1][1];
                $this->arrRet["M7"] = $matches[1][2];
                $this->arrRet["M8"] = $matches[1][3];
                
                chdir($this->sTmpFileBase ); //go to temp directory
                //Run Ctrl file 3 for M8a ================================
                echo("Running Model 8a ...\n");
                
                $sCMD = CODEML_PATH . " ".  CodemlUtil::PREFIX_CTRL_FILE3 . ".ctl";
                //echo("$sCMD\n");
                $sRet = exec($sCMD);
                $sF3 = fnReadText( CodemlUtil::PREFIX_CTRL_FILE3 . ".out");
                
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
                
                //print_r($this->arrRet);
                
                
                
                chdir($this->sCurrDir); //always make sure to start from the working directory.
                if (!$bKeepTmpFolder) {
		//delete the temp files: 
	                foreach(glob($this->sTmpFileBase.'*') as $v){
	                        unlink($v);
	                }
			deleteDirectory($this->sTmpFileBase);
		}
                
                
                return $this->arrRet;
                
        }
        
        public function RunCodemlBranchSite($bKeepTmpFolder=false) {
        
                if (count($this->arrRet) > 0)
                {
                        return $this->arrRet;
                }
                
                $this->initializeBranchSite();
                // Prepare the in files:
                $sFas = "";
		chdir($this->sCurrDir); //always make sure to start from the working directory.
				 
                for ($i=0;$i<count($this->arrSeqs);$i++) {
                        
                                $sFas .= '>'.trim($this->arrLabels[$i])."\n";
                                $sFas .= trim($this->arrSeqs[$i])."\n";
                                                        
                }

                
                fnWriteText( $this->sTmpFileBase . CodemlUtil::IN_FILE_NAME, $sFas);
				$this->arrRet["tmpfolder"] = $this->sTmpFileBase;
                $this->arrRet["Success"] = true;
				
                if (!isset($this->sTreeContent) ) {
                   die("User-specified tree is required for branch-site models. Because you need to label the foreground clades.\n");
                }
 		 else {
                	
                        $sTree = $this->sTreeContent;
                        if (substr(trim($sTree), -1, 1) != ';') {
                        	    $sTree .= ';';             	
                        }
                        fnWriteText( $this->sTmpFileBase . CodemlUtil::TREE_FILE_NAME, $sTree );
                }
                
                
                
 
                chdir($this->sTmpFileBase ); //go to temp directory
                //Run Ctrl file 2 ======================================
                echo("Running Model A with fixed omega (null) ...\n");
                
                $sCMD = CODEML_PATH . " ".  CodemlUtil::PREFIX_CTRL_FILE2 . ".ctl";
                //echo("$sCMD\n");
                $sRet = exec($sCMD);
                $sF2 = fnReadText(  CodemlUtil::PREFIX_CTRL_FILE2 . ".out");
                
                //Parse the file:
                /* $this->arrRet[]; */
                //echo($sF2);
                preg_match_all('/lnL\(ntime\:\s*\S*\s*np\:\s*\S*\)\:\s+(\-[0-9]*\.[0-9]*)/m', $sF2, $matches, PREG_PATTERN_ORDER);
                
                //print_r($matches);
        
                $this->arrRet["MA_null"] = $matches[1][0];
                
                chdir($this->sTmpFileBase ); //go to temp directory
                //Run Ctrl file 3 for M8a ================================
                echo("Running Model A with estimated omega for foreground ...\n");
                
                $sCMD = CODEML_PATH . " ".  CodemlUtil::PREFIX_CTRL_FILE3 . ".ctl";
                //echo("$sCMD\n");
                $sRet = exec($sCMD);
                $sF3 = fnReadText( CodemlUtil::PREFIX_CTRL_FILE3 . ".out");
                
                //Parse the file:
                /* $this->arrRet[]; */
                //echo($sF3);
                
                preg_match_all('/lnL\(ntime\:\s*\S*\s*np\:\s*\S*\)\:\s+(\-[0-9]*\.[0-9]*)/m', $sF3, $matches, PREG_PATTERN_ORDER);
                
                //print_r($matches);
        
                $this->arrRet["MA"] = $matches[1][0];
                
                // do likelihood ratio tests
                echo("Doing LRT ...\n");
                
                $this->arrRet["MA-MAnull"] = $this->LikelihoodRatio( $this->arrRet["MA_null"] , $this->arrRet["MA"] , 1);
                
                //print_r($this->arrRet);
                
                
                
                chdir($this->sCurrDir); //always make sure to start from the working directory.
                if (!$bKeepTmpFolder) {
		//delete the temp files: 
	                foreach(glob($this->sTmpFileBase.'*') as $v){
	                        unlink($v);
	                }
			deleteDirectory($this->sTmpFileBase);
		}
                
                
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


class HyPhyUtil {

        private $nInstanceID ;
        private $sTmpFileBase;
        private $sTreeContent;
        private $arrLabels = array();
        private $arrSeqs = array();
        private $arrRet = array();
	private $sCurrDir = "?";
        
        private $bInit = false;
        
        const IN_FILE_NAME = "in.nex";
        const PREFIX_CTRL_FILE1 = "hyphyrelax";
        
        
        private function initialize() {
        
          if ($this->bInit) {return;};
          
          $this->sCurrDir = ($this->sCurrDir=="?")? getcwd() ."/":$this->sCurrDir;
      $this->nInstanceID = rand( 1000000 , 9999999) ;
          $this->sTmpFileBase = HYPHY_TEMP_PATH."/".$this->nInstanceID."/";
	  chdir($this->sCurrDir);
          mkdir($this->sTmpFileBase, 0777, true); // create temp file folder
          
          $this->sTmpFileBase = realpath(HYPHY_TEMP_PATH."/".$this->nInstanceID) ."/";
          echo("WorkingPath: ".$this->sCurrDir . PHP_EOL);
          echo("TmpPath: ".$this->sTmpFileBase . PHP_EOL);

          //Read the template control file:
          $sCtrlTemplate = fnReadText(HYPHY_CTRL_TEMPLATE);
          
          
          //Common variables shared by all control files:
          $sCtrlTemplate = str_replace("%NEXUS%",  $this->sTmpFileBase . HyPhyUtil::IN_FILE_NAME ,  $sCtrlTemplate);
          $sCtrlTemplate = str_replace("%HYPHY_PATH%",   HYPHY_PATH  ,  $sCtrlTemplate);
          
          //Control file 1 for just getting plain dN/dS:
          $sCtrl1 = $sCtrlTemplate ;
          
          fnWriteText($this->sTmpFileBase .  HyPhyUtil::PREFIX_CTRL_FILE1 . ".bf" , $sCtrl1);
          
  
          
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
        
        public function RunRELAX($bKeepTmpFolder=true) {
        
                if (count($this->arrRet) > 0)
                {
                        return $this->arrRet;
                }
                
                $this->initialize();
                // Prepare the in files:
                $sNexus = "#NEXUS\n";
		$sNexus .= "BEGIN DATA;\n";
		$sNexus .= "DIMENSIONS  NTAX=".count($this->arrSeqs)." NCHAR=".strlen($this->arrSeqs[0]).";\n";
		$sNexus .= "FORMAT DATATYPE = DNA GAP = - MISSING = N;\n";
		$sNexus .= "MATRIX\n";
		chdir($this->sCurrDir); //always make sure to start from the working directory.
				 
                for ($i=0;$i<count($this->arrSeqs);$i++) {
                        
                                $sNexus .= trim($this->arrLabels[$i])."\n";
                                $sNexus .= trim($this->arrSeqs[$i])."\n\n";
                                                        
                }

		$sNexus .= "\n;\nEND;\n\n";

                $sNexus .= "begin trees;\n";
		$sNexus .= "tree tree1 = ".$this->sTreeContent.";\n";
		$sNexus .= "END;\n";

                fnWriteText( $this->sTmpFileBase . HyPhyUtil::IN_FILE_NAME , $sNexus);
				$this->arrRet["tmpfolder"] = $this->sTmpFileBase;
                $this->arrRet["Success"] = true;
				             
                
                
                //Run Ctrl file 1: ========================================
                chdir($this->sTmpFileBase ); //go to temp directory
                echo("Running RELAX ...\n");

                $sCMD = HYPHY_PATH . "/HYPHYMP ".  HyPhyUtil::PREFIX_CTRL_FILE1 . ".bf";
                echo($sCMD . "\n");
                
                //return;

                $sRet = exec($sCMD);

                $sF1 = fnReadText(  $this->sTmpFileBase . HyPhyUtil::IN_FILE_NAME . ".RELAX.json"); //read output.
                
                //Parse the file:
                /* $this->arrRet[]; */
                //echo($sF1);
                // works for paml4.4 preg_match_all('/t\=\s*\S+\s+S\=\s*\S+\s*N\=\s*\S+\s*dN\/dS\=\s*([0-9]*\.[0-9]+)\s*dN\=\s*([0-9]*\.[0-9]+)\s*dS\=\s*(([0-9]*\.[0-9]+))/m', $sF1, $matches, PREG_PATTERN_ORDER);

		$arrJSONret = json_decode($sF1 , true);
                if (is_null($arrJSONret) ) {
					$this->arrRet["MG94xREV"] = ""; //loglik of MG94XREV
					$this->arrRet["NULL"] = ""; //Loglik of null model assuming K=1 (no relax or intensification of selection)
					$this->arrRet["Alternative"] = ""; //loglik of alt model, k estimated.
					$this->arrRet["K"] = ""; //estimated k
					$this->arrRet["LR"] = ""; //likelihood ratio
					$this->arrRet["P"] = ""; //likelihood ratio
					$this->arrRet["R_omega0"] = ""; //omega 0 (<1) of reference branches.
					$this->arrRet["R_omega0_prop"] = ""; // proportion of omega 0 (<1) of reference branches.
					$this->arrRet["R_omega1"] = ""; //omega 1 (near 1) of reference branches.
					$this->arrRet["R_omega1_prop"] = ""; // proportion of omega 1 (near 1) of reference branches.
					$this->arrRet["R_omega2"] = ""; //omega 2 (>1) of reference branches.
					$this->arrRet["R_omega2_prop"] = ""; // proportion of omega 2 (>1) of reference branches.
					$this->arrRet["T_omega0"] = ""; //omega 0 (<1) of Test branches.
					$this->arrRet["T_omega0_prop"] = ""; // proportion of omega 0 (<1) of Test branches.
					$this->arrRet["T_omega1"] = ""; //omega 1 (near 1) of Test branches.
					$this->arrRet["T_omega1_prop"] = ""; // proportion of omega 1 (near 1) of Test branches.
					$this->arrRet["T_omega2"] = ""; //omega 2 (>1) of Test branches.
					$this->arrRet["T_omega2_prop"] = ""; // proportion of omega 2 (>1) of Test branches.
					 $this->arrRet["Success"] = false;
                			chdir($this->sCurrDir); //always make sure to start from the working directory.
					return $this->arrRet;
				}
				

		


		$this->arrRet["MG94xREV"] = $arrJSONret["fits"]["Partitioned MG94xREV"]["log-likelihood"]; //loglik of MG94XREV
		$this->arrRet["NULL"] = $arrJSONret["fits"]["Null"]["log-likelihood"]; //Loglik of null model assuming K=1 (no relax or intensification of selection)
		$this->arrRet["Alternative"] = $arrJSONret["fits"]["Alternative"]["log-likelihood"]; //loglik of alt model, k estimated.
		$this->arrRet["K"] = $arrJSONret["fits"]["Alternative"]["K"]; //estimated k
		$this->arrRet["LR"] = $arrJSONret["relaxation-test"]["LR"]; //likelihood ratio
		$this->arrRet["P"] = $arrJSONret["relaxation-test"]["p"]; //p value with 1 degree of freedom. chi sq.
		$this->arrRet["R_omega0"] = $arrJSONret["fits"]["Alternative"]["rate-distributions"]["Reference"][0][0]; //omega 0 (<1) of reference branches.
		$this->arrRet["R_omega0_prop"] = $arrJSONret["fits"]["Alternative"]["rate-distributions"]["Reference"][0][1]; // proportion of omega 0 (<1) of reference branches.
		$this->arrRet["R_omega1"] = $arrJSONret["fits"]["Alternative"]["rate-distributions"]["Reference"][1][0]; //omega 1 (near 1) of reference branches.
		$this->arrRet["R_omega1_prop"] = $arrJSONret["fits"]["Alternative"]["rate-distributions"]["Reference"][1][1]; // proportion of omega 1 (near 1) of reference branches.
		$this->arrRet["R_omega2"] = $arrJSONret["fits"]["Alternative"]["rate-distributions"]["Reference"][2][0]; //omega 2 (>1) of reference branches.
		$this->arrRet["R_omega2_prop"] = $arrJSONret["fits"]["Alternative"]["rate-distributions"]["Reference"][2][1]; // proportion of omega 2 (>1) of reference branches.
		$this->arrRet["T_omega0"] = $arrJSONret["fits"]["Alternative"]["rate-distributions"]["Test"][0][0]; //omega 0 (<1) of Test branches.
		$this->arrRet["T_omega0_prop"] = $arrJSONret["fits"]["Alternative"]["rate-distributions"]["Test"][0][1]; // proportion of omega 0 (<1) of Test branches.
		$this->arrRet["T_omega1"] = $arrJSONret["fits"]["Alternative"]["rate-distributions"]["Test"][1][0]; //omega 1 (near 1) of Test branches.
		$this->arrRet["T_omega1_prop"] = $arrJSONret["fits"]["Alternative"]["rate-distributions"]["Test"][1][1]; // proportion of omega 1 (near 1) of Test branches.
		$this->arrRet["T_omega2"] = $arrJSONret["fits"]["Alternative"]["rate-distributions"]["Test"][2][0]; //omega 2 (>1) of Test branches.
		$this->arrRet["T_omega2_prop"] = $arrJSONret["fits"]["Alternative"]["rate-distributions"]["Test"][2][1]; // proportion of omega 2 (>1) of Test branches.
                
                      

                
                
                chdir($this->sCurrDir); //always make sure to start from the working directory.
                if (!$bKeepTmpFolder) {
		//delete the temp files: 
	                foreach(glob($this->sTmpFileBase.'*') as $v){
	                        unlink($v);
	                }
			deleteDirectory($this->sTmpFileBase);
		}
                
                
                return $this->arrRet;
                
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

function fnExludeCodonsWithN($arrAln ) {


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
	for ($i=0;$i<$nLen;$i+=3) {

	// first see if they're all equal, if so then don't do complicated things
		//print_r($i);
		$arrCol = array(); // column
		$bKeepCol = true; //keep column?

		for ($l=0;$l < $nTaxa; $l++ ) {
			$arrCol[$l] = strtoupper(substr($arrSeqs[$l] , $i, 3));
			
			if (strpos($arrCol[$l] , 'N')!== false) {
				//print_r($i);
				$bKeepCol = false;
				break;
			}
		}		

		if ( $bKeepCol ) { // copy column to ret
			
			for ($k=0; $k < $nTaxa; $k++ ) {
				$arrRet[$k] .= $arrCol[$k];
			}

		}

	}
	
	

	return array_combine($arrTaxa , $arrRet);

}

function fnCountdNdS($arrAln ) {


	$arrTaxa = array_keys($arrAln);
	$arrSeqs = array_values($arrAln);
	$nTaxa = count($arrTaxa);
	$nLen = strlen($arrSeqs[0]);
	$nTotalCodons = 0;
	$nIdenticalCodons = 0;
	$nMoreThanOneMutate = 0;
	$nSynonymousChange = 0;
	$nNonSynChange = 0;
	
	$nS1 = 0; //syn change at codon pos 1
	$nS2 = 0;
	$nS3 = 0;
	$nN1 = 0; //nsyn change at codon pos 1
	$nN2 = 0;
	$nN3 = 0;
	
	$nM1 = 0; //codon pos 1 mutated in a codon pair with >1 substitutions
	$nM2 = 0;
	$nM3 = 0;
	
	$bStopFound = 0; // 0- not found, 1 - found in 1, 2 - found in 2, 3 found in both.
	
	$arrProteinSeqs = array("", "");
	$arrDNASeqs = array("", "");
	$sDNAChanges = "";
	$sProteinChanges = "";
	
	if ($nTaxa!=2) {
		die("Count dNdS only works for two taxa.\n");
	}	
	// go over each position
	for ($i=0;$i<$nLen;$i+=3) {
		$nTotalCodons++;
	// first see if they're all equal, if so then don't do complicated things
		//print_r($i);
		$arrCol = array(); // column
		$bKeepCol = true; //keep column?
		
				


		for ($l=0;$l < $nTaxa; $l++ ) {
			$arrCol[$l] = strtoupper(substr($arrSeqs[$l] , $i, 3));
		}		
		
		if (strlen($arrCol[0])!=3) {
			$nTotalCodons--;
			continue;
		}
		
		$arrAA1 = fnDNA2Peptide($arrCol[0]);
		$arrAA2 = fnDNA2Peptide($arrCol[1]);
		$sAA1 = $arrAA1[0] ;
		$sAA2 = $arrAA2[0] ;
		$arrProteinSeqs[0] .= $sAA1;
		$arrProteinSeqs[1] .= $sAA2;
		$arrDNASeqs[0] .= $arrCol[0];
		$arrDNASeqs[1] .= $arrCol[1];
		
		if ($sAA1 == '*') {
			$bStopFound = $bStopFound | 1;
		}
		
		if ($sAA2 == '*') {
			$bStopFound = $bStopFound | 2;
		}
		
		if ($arrCol[0] == $arrCol[1]) {
			$nIdenticalCodons++;
			$sProteinChanges .= ".";
			$sDNAChanges .= "...";
			continue;
		}
		
		if (levenshtein($arrCol[0] , $arrCol[1]) > 1) { //more than one base change
			$nMoreThanOneMutate++;
			$sProteinChanges .= "!";
			
			if (substr($arrCol[0], 0,1) != substr($arrCol[1], 0,1)) {
				$nM1++;
				$sDNAChanges .= "A";
			} 
			else {
				$sDNAChanges .= ".";
			}
			
			if (substr($arrCol[0], 1,1) != substr($arrCol[1], 1,1)) {
				$nM2++;
				$sDNAChanges .= "B";
			} 
			else {
				$sDNAChanges .= ".";
			}
			if (substr($arrCol[0], 2,1) != substr($arrCol[1], 2,1)) {
				$nM3++;
				$sDNAChanges .= "C";
			} 
			else {
				$sDNAChanges .= ".";
			}
			
			continue;
		}

		
		if ( $sAA1 == $sAA2) {
			$nSynonymousChange++;
			$sProteinChanges .= "S";
			if (substr($arrCol[0], 0,1) != substr($arrCol[1], 0,1)) {
				$nS1++;
				$sDNAChanges .= "1..";
			} 
			else if (substr($arrCol[0], 2,1) != substr($arrCol[1], 2,1)) {
				$nS3++;
				$sDNAChanges .= "..3";
			} else  { // this is impossible, but still here for the sake of completness
				$nS2++;
				$sDNAChanges .= ".2.";
			}
			continue;
		}
		else {
			$nNonSynChange++;
			$sProteinChanges .= "N";
			if (substr($arrCol[0], 1,1) != substr($arrCol[1], 1,1)) {
				$nN2++;
				$sDNAChanges .= ".b.";
			} 
			else if (substr($arrCol[0], 0,1) != substr($arrCol[1], 0,1)) {
				$nN1++;
				$sDNAChanges .= "a..";
			} else  { 
				$nN3++;
				$sDNAChanges .= "..c";
			}
			

			continue;
		}


	}
	
	

	return array("stopcodon" => $bStopFound, "M1" => $nM1, "M2"=>$nM2, "M3"=>$nM3, "S1" => $nS1, "S2" => $nS2 , "S3" => $nS3, "N1"=>$nN1, "N2" => $nN2, "N3" =>$nN3, "DNAs" => $arrDNASeqs, "DNA_changes" => $sDNAChanges, "proteins" => $arrProteinSeqs, "aa_changes" => $sProteinChanges, "total_codons"=>$nTotalCodons , "identical_codons" => $nIdenticalCodons, "morethanonebase_change" =>$nMoreThanOneMutate, "syn"=>$nSynonymousChange, "nsyn"=>$nNonSynChange);

}

function fnExcludeLastStop($arrFullCodingSeq) {

	$bDoExclude = false; //if exclude for one taxon, also need to exclude others to make sure seq same length 

	foreach($arrFullCodingSeq as $sTaxon => $sSeq) {
		$sLastCodon = substr($sSeq , -3);
		$sAAlast = fnDNA2Peptide($sLastCodon, false, false)[0];
		if ($sAAlast == "*" ) {
			$bDoExclude = true;
		}
	}

	if ($bDoExclude) {
	foreach($arrFullCodingSeq as $sTaxon => $sSeq) {

		
			$arrFullCodingSeq[$sTaxon] = substr($sSeq , 0, strlen($sSeq)-3);
	}
	}
	
	return $arrFullCodingSeq;

}

function fnTranslateAlignment($arrAln, $bExcludeLastStop=false) {
	
	$arrRet = array();
	$bStopFound = false;
	$arrTaxaWithStop=array();
	foreach($arrAln as $sTaxon => $sSeq) {
		$sPeptide = fnDNA2Peptide($sSeq, false, false)[0];
		if (strpos($sPeptide, "*") !== false) {
			if (strpos($sPeptide, "*") == (strlen($sPeptide)-1)) { //this is at last of its position which is normal
				if ($bExcludeLastStop) {
					$sPeptide = substr($sPeptide , 0, strlen($sPeptide)-1 );
		
				} else {

					$bStopFound = true;
					$arrTaxaWithStop[] = $sTaxon;
				}
			} else {
			
				$bStopFound = true;
				$arrTaxaWithStop[] = $sTaxon;
			}
		}
		$arrRet[$sTaxon] = $sPeptide;
	}
	
	return array("stopcodon" => $bStopFound, "stoptaxa" => $arrTaxaWithStop, "proteins"=>$arrRet);
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

function fnIs4FoldDegen($sCodon) {
	global $Codon4FoldDegeneracyMap;
	$sCodon = strtoupper($sCodon);

	preg_match_all('/[^ATCG]/', $sCodon, $arrM, PREG_OFFSET_CAPTURE); //check for ambiguous code

	$arrM = $arrM[0];
	if (count($arrM) > 1 ) {
		return -1; //if more than 1 site is ambiguous, do not try to solve the ambiguity
	} else if (count($arrM) == 1) {
		$sAmbBase = $arrM[0][0];
		$sAmbPos = $arrM[0][1];

		if ($sAmbPos != 2) {
			return -1; // the ambiguity base is not on the third position, do not try to resolve.
		} else {
			$sCodon = substr($sCodon, 0, 2) . "T"; // arbitrarily put T, if it is 4 fold will be fine
		}
	}


	if (!array_key_exists($sCodon , $Codon4FoldDegeneracyMap)) {
		return -1;
	}

	return $Codon4FoldDegeneracyMap[$sCodon];
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

function deleteDirectory($dir) {
    if (!file_exists($dir)) return true;
    if (!is_dir($dir)) return unlink($dir);
    foreach (scandir($dir) as $item) {
        if ($item == '.' || $item == '..') continue;
        if (!deleteDirectory($dir.DIRECTORY_SEPARATOR.$item)) return false;
    }
    return rmdir($dir);
}

?>
