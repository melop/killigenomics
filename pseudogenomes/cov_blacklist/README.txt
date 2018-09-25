Scan the bam files at the exonic regions (using the GFF annotations), 
if average coverage at the exonic regions of a gene is significantly higher than genomic background, produce a black list
these regions are likely due to mismapping of different paraglogous reads.

These files are used by ../../hyphyrelax/export_protein_multiref

1. Run hhvm driver.php to obtain coverage statistics at each protein coding gene for every species. The output folder format is [REF_SP]/[SP]
2. Run Rscript identify_highcov.R to calculate p values for each gene, saved into [REF_SP]/[SP]/out_Pvalues.txt. 
3. Optionally produce a summary list with orthogroup_blacklist.php 
