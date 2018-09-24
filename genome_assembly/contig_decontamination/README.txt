1. perform blast to identify contaminants using blastcontigs.sh, output: contig.megablast.txt
2. Run getoriginalcontigname_by_number.sh to obtain original_congitnames.txt
3. Run php idcontainate_contig.php to output a list of likely contaminated contigs, obtain contaminated.list.txt
4. Run php filtercontigs.php to output the filtered contig file to obtain two files: contig.fa.decontaminated.fa (decontaminated) and contig.fa.contaminants.fa (contaminants)


