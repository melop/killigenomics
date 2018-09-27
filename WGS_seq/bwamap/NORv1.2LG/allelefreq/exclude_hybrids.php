<?php
//structure analysis reveals some samples to be hybrids or mislabeled (n=4)
//exclude them from the My_in_GFE.txt file in each pop folder
//and need to rerun GFE to estimate allele freq 

$arrExcludeID = array('ORTDRY_414_NO-414-S7_D502_D705', 'RACWET_92_NR-92-R11_S503_D707', 'RACWET_96_NR-96-Y4_S503_D709', 'RACDRY_215_NR-215-UP_S507_N701');
$arrPopFolders = array('ORTDRY', 'RACWET', 'RACDRY');

$arrExcludeID = array_flip($arrExcludeID);
$arrExcludeCol = array(); //column index to exclude

foreach($arrPopFolders as $sFolder) {
	$hF = fopen($sFolder."/My_in_GFE.txt" , 'r');
	$sOutDIR = "$sFolder.excludehybrid";
	exec("mkdir -p $sOutDIR");
	$hO = fopen($sOutDIR."/My_in_GFE.txt", 'w');
	while( false !== ($sLn = fgets($hF) ) ) {
		$sLn = trim($sLn);
		if ($sLn == '') continue;

		$arrF = explode("\t" , $sLn);
		
		if ($arrF[0] == 'scaffold') {
			foreach($arrF as $nCol => &$sVal) {
				if (array_key_exists($sVal , $arrExcludeID) ) {
					$arrExcludeCol[] = $nCol;
					echo("Column $nCol / $sVal will be exlucded from $sFolder\n");
				}
			}
		}

		foreach($arrExcludeCol as $nExCol) {
			unset($arrF[$nExCol]);
		}

		fwrite($hO , implode("\t", $arrF) . "\n");
	}
}

?>
