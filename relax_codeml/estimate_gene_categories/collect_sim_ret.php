<?php
//collect simulation results and makes the sim_ret.txt table
/*
simulated	relax_k<1	relax_k>1	codeml	intersect(relax_k<1 , codeml)
pos		a		f		k		p
relax		b		g		l		q
intensified	c		h		m		r
neutral		d		i		n		s
relax_pos	e		j		o		t

*/

$sGenus = "Nothobranchius";
$arrSimFolders = array( 'pos' => "../sim_nothos_pos" , //positive selection sims.
			'relax' => "../sim_nothos_relaxed",
			'intensified' => "../sim_nothos_intensified",
			'neutral' => "../sim_nothos_neutral",
			'relax_pos' => "../sim_nothos_relax_pos",
			'intense_pos' => "../sim_nothos_intense_pos"
		);

/*$sGenus = "Aphyosemion";
$arrSimFolders = array( 'pos' => "../sim_aphyosemion_pos" , //positive selection sims.
			'relax' => "../sim_aphyosemion_relaxed",
			'intensified' => "../sim_aphyosemion_intensified",
			'neutral' => "../sim_aphyosemion_neutral",
			'relax_pos' => "../sim_aphyosemion_relax_pos",
			'intense_pos' => "../sim_aphyosemion_intense_pos"
		);*/

/*$sGenus = "Scriptaphyosemion";
$arrSimFolders = array( 'pos' => "../sim_scriptaphyosemion_pos" , //positive selection sims.
			'relax' => "../sim_scriptaphyosemion_relaxed",
			'intensified' => "../sim_scriptaphyosemion_intensified",
			'neutral' => "../sim_scriptaphyosemion_neutral",
			'relax_pos' => "../sim_scriptaphyosemion_relax_pos",
			'intense_pos' => "../sim_scriptaphyosemion_intense_pos"
		);*/

/*$sGenus = "Callopanchax";
$arrSimFolders = array( 'pos' => "../sim_callopanchax_pos" , //positive selection sims.
			'relax' => "../sim_callopanchax_relaxed",
			'intensified' => "../sim_callopanchax_intensified",
			'neutral' => "../sim_callopanchax_neutral",
			'relax_pos' => "../sim_callopanchax_relax_pos",
			'intense_pos' => "../sim_callopanchax_intense_pos"
		);*/


$nPCutoffRelax = 0.05;
$nPCutoffCodeml = 0.05;
$bOnlyReal = false;



//$arrRealFiles = array("../relax_46_spp/rerun_omega0/sum_Nothobranchius.txt", "../codeml_46_spp/sum_Nothobranchius.txt");

$arrRealFiles = array("../relax_46_spp/rerun_omega0/Relax_$sGenus/ret*", "../codeml_46_spp/sum_$sGenus.txt");

$sOutTableReal = "real_ret_$sGenus.txt";
$hOutTableReal = fopen($sOutTableReal , 'w'); 


$oReal = fnCount("cat ".$arrRealFiles[0], "cat ".$arrRealFiles[1] );
fwrite($hOutTableReal , implode("\t" , array_keys($oReal) )."\n".implode("\t" , array_values($oReal) )."\n");

if ($bOnlyReal) die();


$sOutTable = "sum_ret_$sGenus.txt";


$hOutTable = fopen($sOutTable , 'w'); 


fwrite($hOutTable , "simulation\trelax_k_less_1\trelax_k_larger_1\tcodeml\tintersect_relax_codeml\tintersect_intense_codeml\n");
foreach($arrSimFolders as $sSimName => $sSimFolder) {
	echo("=============== $sSimName ====================\n");
	$oSimCounts = fnCount("cat $sSimFolder/relax/rerun_omega0/Relax_*/ret*", "cat $sSimFolder/codeml/Codeml_*/ret*");
	$nTotalExcl = $oSimCounts['relaxtotal'] - $oSimCounts['relaxexl'] - $oSimCounts['intenseexl'];
	fwrite($hOutTable , "$sSimName\t".($oSimCounts['relaxsig'] / $nTotalExcl)."\t".($oSimCounts['intensesig'] / $nTotalExcl) ."\t".($oSimCounts['codemlsig']/$oSimCounts['codemltotal'])."\t".($oSimCounts['overlapposrelax']/$oSimCounts['overlapall'])."\t".($oSimCounts['overlapposintense']/$oSimCounts['overlapall'])."\n");
}

function fnCount($sRelaxOutCmd, $sCodemlOutCmd) {
	global $nPCutoffRelax, $nPCutoffCodeml;

	/*
	$nRelaxTotal = fnGetOutput("$sRelaxOutCmd | awk '{if ($2==\"Success\") print $0; }' | cut -f1,1 | sort | uniq > relax_all.txt; cat relax_all.txt | wc -l ;");
	$nRelaxSig = fnGetOutput("$sRelaxOutCmd | awk '{if ( $2==\"Success\" && $11!=0 && $9 < 1 && $5 <$nPCutoffRelax ) print $0; }' | cut -f1,1 | sort | uniq > relax.txt; cat relax.txt | wc -l ;");
	$nRelaxExcluded = fnGetOutput("$sRelaxOutCmd | awk '{if (  $2==\"Success\" && $11==0 && $9 < 1 && $5 <$nPCutoffRelax ) print $0; }' | cut -f1,1 | sort | uniq | wc -l");
	$nIntenseSig = fnGetOutput("$sRelaxOutCmd | awk '{if (  $2==\"Success\" && $11!=0 && $9 > 1 && $5 <$nPCutoffRelax ) print $0; }'  | cut -f1,1 | sort | uniq > intense.txt; cat intense.txt | wc -l ;");
	$nIntenseExcluded = fnGetOutput("$sRelaxOutCmd | awk '{if (  $2==\"Success\" && $11==0 && $9 > 1 && $5 <$nPCutoffRelax ) print $0; }'  | cut -f1,1 | sort | uniq | wc -l");
	*/

	$nRelaxTotal = fnGetOutput("$sRelaxOutCmd | awk '{if ($2==\"Success\") print $0; }' | cut -f1,1 | sort | uniq > relax_all.txt; cat relax_all.txt | wc -l ;");
	$nRelaxSig = fnGetOutput("$sRelaxOutCmd | awk '{if ( $2==\"Success\" &&  $9 < 1 && $5 <$nPCutoffRelax ) print $0; }' | cut -f1,1 | sort | uniq > relax.txt; cat relax.txt | wc -l ;");
	$nRelaxExcluded = 0;
	$nIntenseSig = fnGetOutput("$sRelaxOutCmd | awk '{if (  $2==\"Success\" &&  $9 > 1 && $5 <$nPCutoffRelax ) print $0; }'  | cut -f1,1 | sort | uniq > intense.txt; cat intense.txt | wc -l ;");
	$nIntenseExcluded = 0;

	echo("relax sig (after excl) $nRelaxSig, excl $nRelaxExcluded, total(before excl) $nRelaxTotal\n");
	echo("inten sig (after excl) $nIntenseSig, excl $nIntenseExcluded, total(before excl) $nRelaxTotal\n");


	$nCodemlTotal = fnGetOutput("$sCodemlOutCmd | awk '{if ($2==\"Success\") print $0; }' | cut -f1,1 | sort | uniq  > pos_all.txt; cat pos_all.txt | wc -l ;");
	$nCodemlSig = fnGetOutput("$sCodemlOutCmd | awk '{if (  $2==\"Success\" && $5 <$nPCutoffCodeml ) print $0; }' | cut -f1,1 | sort | uniq > pos.txt; cat pos.txt | wc -l ;");
	echo("codeml: $nCodemlSig / $nCodemlTotal = ".($nCodemlSig / $nCodemlTotal).PHP_EOL);

	$nOverlapAll = fnGetOutput("  comm -12 pos_all.txt relax_all.txt | wc -l");
	$nOverlapPosRelax = fnGetOutput(" comm -12 pos.txt relax.txt | wc -l");
	$nOverlapPosIntense = fnGetOutput(" comm -12 pos.txt intense.txt | wc -l");

	echo("overlap: $nOverlapPosRelax / $nOverlapAll\n");

	return array('relaxsig'=>$nRelaxSig, 'relaxexl' => $nRelaxExcluded, 'relaxtotal' => $nRelaxTotal, 'intensesig' => $nIntenseSig , 'intenseexl'=>$nIntenseExcluded, 'codemlsig' => $nCodemlSig, 'codemltotal' => $nCodemlTotal, 'overlapposrelax' =>$nOverlapPosRelax, 'overlapposintense' => $nOverlapPosIntense ,  'overlapall'=>$nOverlapAll);

	//return array('relaxsig'=>$nRelaxSig+$nRelaxExcluded, 'relaxexl' => 0, 'relaxtotal' => $nRelaxTotal, 'intensesig' => $nIntenseSig + $nIntenseExcluded , 'intenseexl'=>0, 'codemlsig' => $nCodemlSig, 'codemltotal' => $nCodemlTotal, 'overlapposrelax' =>$nOverlapPosRelax, 'overlapall'=>$nOverlapAll);

	//return array('relaxsig'=>$nRelaxSig+$nIntenseExcluded, 'relaxexl' => 0, 'relaxtotal' => $nRelaxTotal, 'intensesig' => $nIntenseSig + $nRelaxExcluded , 'intenseexl'=>0, 'codemlsig' => $nCodemlSig, 'codemltotal' => $nCodemlTotal, 'overlapposrelax' =>$nOverlapPosRelax, 'overlapall'=>$nOverlapAll);

}

function fnGetOutput($sCmd) {
	$p = popen($sCmd , 'r');
	$s = fgets($p);
	pclose($p);
	return intval($s);
}

?>
