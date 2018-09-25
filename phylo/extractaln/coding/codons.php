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
?>
