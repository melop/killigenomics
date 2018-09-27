#!/usr/bin/perl

# Updated on 05/19/15
# Ext_Ref_Nuc.pl to
# extract the nucleotide at
# every position from the
# reference sequence

use strict; 
use warnings;

# Declare and initialize variables
my $in_file_name = '';
my $out_file_name = '';
my $line = '';
my $num_scaffolds = '';
my $pos_gt = '';
my $scaf_num = '';
my @scaf = ();
my @seq = ();
my $line_num = '';
my $i = '';
my $scaffold = '';
my $j = '';
my $pos = '';
my $ref_nuc = '';

die "Usage: Ext_Ref_Nuc.pl <in_file_name> <out_file_name>\n" if (@ARGV != 2);

# Read the specified settings
$in_file_name = $ARGV[0];
$out_file_name = $ARGV[1];

# Open the input file of the reference sequence
open(IN, "<$in_file_name")  ||  die(" There is no $in_file_name ");

# Read the data from the file in a while loop,
# remove newlines, and add the line to the
# output file as it is read
$scaf_num = 0;
while( $line = <IN> ) {
	if ($line =~ /^>/) {
		# Remove the new line at the end of $line
                chomp $line;
		
		$pos_gt = index($line, '>');
		$scaf_num = $scaf_num + 1;
		$scaf[$scaf_num] = substr($line, $pos_gt+1);
	} else { 
		# Remove the new line at the end of $line
		chomp $line;
	
		# Add the line to the sequence
		$seq[$scaf_num] .= $line;
	}
}
# Close the input file
close IN;

# Open the output file
open(OUT, ">$out_file_name") || die("Cannot make $out_file_name");

# Print out the field names
print OUT "scaffold\tsite\tref\n";
# print "scaffold\tsite\tref\n";

for ($i = 1; $i <= $scaf_num; $i++) {
	$scaffold = $scaf[$i];

	# Extract the nucleotide from the reference sequence at each position
	
	for ($j = 0; $j < length($seq[$i]); $j++) {;
		$ref_nuc = substr($seq[$i], $j, 1);
		$pos = $j + 1; 
		print OUT "$scaffold\t$pos\t$ref_nuc\n";
		# print "$scaffold\t$pos\t$ref_nuc\n"; 
	}
}

# Close the output file
close OUT;
