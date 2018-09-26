This folder contains scripts to repeat the phylogenetic reconstruction for only the Aplocheiliae and Nothobranchiidae, excluding the more distant outgroups.
This process reuses alignments used for the CodeML scan (~13,000 genes).

1. Run linkfastas.php to link to the fasta files.
2. Run make_alignments.php to concatenate all fasta files into one file aln/all.aln.fas
3. Run 4foldsites.php to extract 4-fold sites from aln/all.aln.fas and saved to aln/all.4fold.aln.phy
4. Proceed to raxml/ to run RAxML on the full alignment
5. Proceed to raxml_4fold/ to run RAxML on only 4-fold sites, and perform molecular clock test
