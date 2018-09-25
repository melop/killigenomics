1. Run getallsites_fasta.php to produce a list containing all protein-coding sites.
2. Optionally, run get4foldsites_fasta.php to produce a list containing only 4-fold sites.
3. Use exportalignments.R to convert the above lists to fasta alignments, with all sites for each species concatenated into one big sequence.
4. Proceed to raxml_allcodons to run RAxML
