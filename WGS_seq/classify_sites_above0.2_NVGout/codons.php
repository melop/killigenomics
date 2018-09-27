<?php
	
	$CodonMap = array(
	"TTT" => "F",
	"TTC" => "F",
	"TTA" => "L",
	"TTG" => "L",
	
	"CTT" => "L",
	"CTC" => "L",
	"CTA" => "L",
	"CTG" => "L",
	
	"ATT" => "I",
	"ATC" => "I",
	"ATA" => "I",
	"ATG" => "M",
	
	"GTT" => "V",
	"GTC" => "V",
	"GTA" => "V",
	"GTG" => "V",
	
	"TCT" => "S",
	"TCC" => "S",
	"TCA" => "S",
	"TCG" => "S",
	
	"CCT" => "P",
	"CCC" => "P",
	"CCA" => "P",
	"CCG" => "P",
	
	"ACT" => "T",
	"ACC" => "T",
	"ACA" => "T",
	"ACG" => "T",
	
	"GCT" => "A",
	"GCC" => "A",
	"GCA" => "A",
	"GCG" => "A",
	
	"TAT" => "Y",
	"TAC" => "Y",
	"TAA" => "*",
	"TAG" => "*",
	
	"CAT" => "H",
	"CAC" => "H",
	"CAA" => "Q",
	"CAG" => "Q",
	
	"AAT" => "N",
	"AAC" => "N",
	"AAA" => "K",
	"AAG" => "K",
	
	"GAT" => "D",
	"GAC" => "D",
	"GAA" => "E",
	"GAG" => "E",
	
	"TGT" => "C",
	"TGC" => "C",
	"TGA" => "*",
	"TGG" => "W",
	
	"CGT" => "R",
	"CGC" => "R",
	"CGA" => "R",
	"CGG" => "R",
	
	"AGT" => "S",
	"AGC" => "S",
	"AGA" => "R",
	"AGG" => "R",
	
	"GGT" => "G",
	"GGC" => "G",
	"GGA" => "G",
	"GGG" => "G",
	
	"NNN" => "?",
	
	"|||" => "|" //this is a hsp divider
	);


	$Codon4FoldDegeneracyMap = array(
	"TTT" => false,
	"TTC" => false,
	"TTA" => false,
	"TTG" => false,
	
	"CTT" => true,
	"CTC" => true,
	"CTA" => true,
	"CTG" => true,
	
	"ATT" => false,
	"ATC" => false,
	"ATA" => false,
	"ATG" => false,
	
	"GTT" => true,
	"GTC" => true,
	"GTA" => true,
	"GTG" => true,
	
	"TCT" => true,
	"TCC" => true,
	"TCA" => true,
	"TCG" => true,
	
	"CCT" => true,
	"CCC" => true,
	"CCA" => true,
	"CCG" => true,
	
	"ACT" => true,
	"ACC" => true,
	"ACA" => true,
	"ACG" => true,
	
	"GCT" => true,
	"GCC" => true,
	"GCA" => true,
	"GCG" => true,
	
	"TAT" => false,
	"TAC" => false,
	"TAA" => false,
	"TAG" => false,
	
	"CAT" => false,
	"CAC" => false,
	"CAA" => false,
	"CAG" => false,
	
	"AAT" => false,
	"AAC" => false,
	"AAA" => false,
	"AAG" => false,
	
	"GAT" => false,
	"GAC" => false,
	"GAA" => false,
	"GAG" => false,
	
	"TGT" => false,
	"TGC" => false,
	"TGA" => false,
	"TGG" => false,
	
	"CGT" => true,
	"CGC" => true,
	"CGA" => true,
	"CGG" => true,
	
	"AGT" => false,
	"AGC" => false,
	"AGA" => false,
	"AGG" => false,
	
	"GGT" => true,
	"GGC" => true,
	"GGA" => true,
	"GGG" => true,
	);


	$CodonNonSyn1stBaseChange = array( //if 1st base changes, is this codon non syn? 
	"TTT" => true,
	"TTC" => true,
	"TTA" => -1, //ambiguous
	"TTG" => -1, 
	
	"CTT" => true,
	"CTC" => true,
	"CTA" => -1,
	"CTG" => -1,
	
	"ATT" => true,
	"ATC" => true,
	"ATA" => true,
	"ATG" => true,
	
	"GTT" => true,
	"GTC" => true,
	"GTA" => true,
	"GTG" => true,
	
	"TCT" => true,
	"TCC" => true,
	"TCA" => true,
	"TCG" => true,
	
	"CCT" => true,
	"CCC" => true,
	"CCA" => true,
	"CCG" => true,
	
	"ACT" => true,
	"ACC" => true,
	"ACA" => true,
	"ACG" => true,
	
	"GCT" => true,
	"GCC" => true,
	"GCA" => true,
	"GCG" => true,
	
	"TAT" => true,
	"TAC" => true,
	"TAA" => true,
	"TAG" => true,
	
	"CAT" => true,
	"CAC" => true,
	"CAA" => true,
	"CAG" => true,
	
	"AAT" => true,
	"AAC" => true,
	"AAA" => true,
	"AAG" => true,
	
	"GAT" => true,
	"GAC" => true,
	"GAA" => true,
	"GAG" => true,
	
	"TGT" => true,
	"TGC" => true,
	"TGA" => true,
	"TGG" => true,
	
	"CGT" => true,
	"CGC" => true,
	"CGA" => true,
	"CGG" => true,
	
	"AGT" => true,
	"AGC" => true,
	"AGA" => true,
	"AGG" => true,
	
	"GGT" => true,
	"GGC" => true,
	"GGA" => true,
	"GGG" => true,
	);

	$CodonNonSyn2ndBaseChange = array( //if 2nd base changes, is this codon non syn? 
	"TTT" => true,
	"TTC" => true,
	"TTA" => true, //ambiguous
	"TTG" => true, 
	
	"CTT" => true,
	"CTC" => true,
	"CTA" => true,
	"CTG" => true,
	
	"ATT" => true,
	"ATC" => true,
	"ATA" => true,
	"ATG" => true,
	
	"GTT" => true,
	"GTC" => true,
	"GTA" => true,
	"GTG" => true,
	
	"TCT" => true,
	"TCC" => true,
	"TCA" => true,
	"TCG" => true,
	
	"CCT" => true,
	"CCC" => true,
	"CCA" => true,
	"CCG" => true,
	
	"ACT" => true,
	"ACC" => true,
	"ACA" => true,
	"ACG" => true,
	
	"GCT" => true,
	"GCC" => true,
	"GCA" => true,
	"GCG" => true,
	
	"TAT" => true,
	"TAC" => true,
	"TAA" => true,
	"TAG" => true,
	
	"CAT" => true,
	"CAC" => true,
	"CAA" => true,
	"CAG" => true,
	
	"AAT" => true,
	"AAC" => true,
	"AAA" => true,
	"AAG" => true,
	
	"GAT" => true,
	"GAC" => true,
	"GAA" => true,
	"GAG" => true,
	
	"TGT" => true,
	"TGC" => true,
	"TGA" => true,
	"TGG" => true,
	
	"CGT" => true,
	"CGC" => true,
	"CGA" => true,
	"CGG" => true,
	
	"AGT" => true,
	"AGC" => true,
	"AGA" => true,
	"AGG" => true,
	
	"GGT" => true,
	"GGC" => true,
	"GGA" => true,
	"GGG" => true,
	);


	$CodonNonSyn3rdBaseChange = array( //if 3rd base changes, is this codon non syn? 
	"TTT" => -1,
	"TTC" => -1,
	"TTA" => -1, //ambiguous
	"TTG" => -1, 
	
	"CTT" => false,
	"CTC" => false,
	"CTA" => false,
	"CTG" => false,
	
	"ATT" => -1,
	"ATC" => -1,
	"ATA" => -1,
	"ATG" => true,
	
	"GTT" => false,
	"GTC" => false,
	"GTA" => false,
	"GTG" => false,
	
	"TCT" => false,
	"TCC" => false,
	"TCA" => false,
	"TCG" => false,
	
	"CCT" => false,
	"CCC" => false,
	"CCA" => false,
	"CCG" => false,
	
	"ACT" => false,
	"ACC" => false,
	"ACA" => false,
	"ACG" => false,
	
	"GCT" => false,
	"GCC" => false,
	"GCA" => false,
	"GCG" => false,
	
	"TAT" => -1,
	"TAC" => -1,
	"TAA" => -1,
	"TAG" => -1,
	
	"CAT" => -1,
	"CAC" => -1,
	"CAA" => -1,
	"CAG" => -1,
	
	"AAT" => -1,
	"AAC" => -1,
	"AAA" => -1,
	"AAG" => -1,
	
	"GAT" => -1,
	"GAC" => -1,
	"GAA" => -1,
	"GAG" => -1,
	
	"TGT" => -1,
	"TGC" => -1,
	"TGA" => true,
	"TGG" => true,
	
	"CGT" => false,
	"CGC" => false,
	"CGA" => false,
	"CGG" => false,
	
	"AGT" => -1,
	"AGC" => -1,
	"AGA" => -1,
	"AGG" => -1,
	
	"GGT" => false,
	"GGC" => false,
	"GGA" => false,
	"GGG" => false,
	);

?>
