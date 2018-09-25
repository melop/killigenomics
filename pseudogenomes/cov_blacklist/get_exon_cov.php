<?php
ini_set('memory_limit', -1);

require_once(dirname(__FILE__) . "/lib.php");

$sMakerGFF = "/beegfs/group_dv/home/RCui/killifish_genomes/annotation/maker/AAU_final_1/run_blastmerge_cov_exonerateOn/corrections/maker.finalannot.improved.gff";
$sGenome = "/beegfs/group_dv/home/RCui/killifish_genomes/pseudogenomes/refs/AAU.fa";
$sOutDIR = "AAU";
$sSp = "AKN1";
$sBAM = "/beegfs/group_dv/home/RCui/killifish_genomes/pseudogenomes/out/bams/AAU/AKN1/sorted.bam";

$nMinMapQ = 0; //do not use mapq cutoff when checking coverage


while(count($argv) > 0) {
    $arg = array_shift($argv);
    switch($arg) {
        case '-G':
            $sMakerGFF = trim(array_shift($argv));
            break;
        case '-g':
            $sGenome  = trim(array_shift($argv));
            break;
        case '-o':
            $sOutDIR  = trim(array_shift($argv));
            break;
        case '-S':
            $sSp  = trim(array_shift($argv));
            break;
        case '-b':
            $sBAM  = trim(array_shift($argv));
            break;

    }
}




$sOutSpDIR = "$sOutDIR/$sSp/";

exec("mkdir -p $sOutSpDIR");

$oGFF = new MappedCDSGFF();
$oGFF->LoadGFF($sMakerGFF , $sGenome);

$arrGenes = $oGFF->GetGenes('maker');

$sExonBed = "$sOutSpDIR/exon.bed";
$sExonSortedBed = "$sOutSpDIR/exon.sorted.bed";
$hExonBed = fopen($sExonBed, 'w');
$arrSpGeneId2Coord = array();

foreach($arrGenes as $sGeneId => $oGene) {
	foreach($oGene['mRNAs'] as $sRNAId => $oRNA) {

		$sScf = $oRNA['scf'];

		if (!array_key_exists('spgeneid' , $oRNA['annot'])) { // if there is no species gene id, skip
			continue;
		}

		$sSpGeneId = $oRNA['annot']['spgeneid'];

		$arrSpGeneId2Coord[$sSpGeneId ] = array($sScf , $oRNA['start'] , $oRNA['end'], $oRNA['strand']);

		foreach($oRNA['CDSs'] as $nStart => $nEnd) {
			$nLeft = ($nStart < $nEnd)? ($nStart-1) : ($nEnd-1);
			$nRight = ($nStart > $nEnd)? $nStart : $nEnd;
			fwrite($hExonBed , "$sScf\t$nLeft\t$nRight\t$sSpGeneId\n");
		}
	}
}

$sHistogram = "$sOutSpDIR/ret.txt";
$sHistogramDone = "$sOutSpDIR/bedtools.done";
exec("awk -v OFS='\t' {'print $1,$2'}  $sGenome.fai > $sOutSpDIR/genomefile.txt");
exec('scfprefix=`head -n 1 '. $sExonBed .' | cut -f 1,1 | sed -r "s/[0-9]+$//"`; sed -r "s/\S+scf//" '.$sExonBed.' | sort -k1,1n -k2,2n | awk -v pf="$scfprefix" \'{ print pf$0 }\' > '.$sExonSortedBed);
exec("if [ -e $sHistogramDone ]; then echo $sHistogramDone exists, skip; else samtools1.2 view -u -q $nMinMapQ $sBAM | bedtools coverage -hist -a $sExonSortedBed -b stdin -sorted -g $sOutSpDIR/genomefile.txt > $sHistogram && touch $sHistogramDone; fi");

if (!file_exists($sHistogramDone) ) {
	die("Failed to count depth using bedtools!\n");
}

//now compute average coverage per exon
$hHist = fopen($sHistogram , 'r');
$sPrevExonID = "";
$sExonID = "";
$nExonLen = 0;
$nExonCov = 0;
$sGeneId = "";
$sPrevGeneId = "";
$nCov = 0;
$nCovBaseCount = 0;
$nExonLen = 0;
$arrF = array();

$arrGenes = array();

do {
	$sLn = fgets($hHist);

	if ($sLn !== false) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;
		$arrF = explode("\t" , $sLn);
		$sExonID = $arrF[0]."_".$arrF[1]."_".$arrF[2];
	}

	if ($sPrevExonID != $sExonID || $sLn === false || $arrF[0] == 'all' ) {
		//echo("Trigger write\n $sPrevExonID != $sExonID ");
		if ($sPrevExonID != '') {
			if (!array_key_exists($sGeneId, $arrGenes) ) {
				$arrGenes[$sGeneId] = array();
			}
			$arrGenes[$sGeneId][] = array( $nCovBaseCount , $nExonLen);
			$nCov = 0;
			$nCovBaseCount = 0;
			$nExonLen = 0;

		}

		if ($sLn === false || $arrF[0] == 'all') {
			break;
		}

	}

	if (count($arrF) >= 7) {
		$sGeneId = $arrF[3];
		$nCov = $arrF[4];
		$nCovBaseCount += $arrF[5] * $nCov;
		$nExonLen = $arrF[6];
	}

	$sPrevExonID = $sExonID;
} while (true);

//print_r($arrGenes);
$hOutTable = fopen($sOutSpDIR."/out.txt" , 'w');
foreach($arrGenes as $sGeneId => $arrExon) {

	$arrInfo = $arrSpGeneId2Coord[$sGeneId];
	$nTotalCov = 0;
	$nTotalLen = 0;
	foreach($arrExon as $oExon) {
		$nTotalCov += $oExon[0];
		$nTotalLen += $oExon[1];
	}

	$nDepth = $nTotalCov / $nTotalLen;
	fwrite($hOutTable , $sGeneId."\t".implode("\t" , $arrInfo). "\t$nTotalCov\t$nTotalLen\t$nDepth\n" );
	
}

?>
