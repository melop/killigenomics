Calculate codon usage bias (codon adaptive index) given a list of highly expressed genes as the usage reference.

1. Run hhvm countCodonUsage.php to count codon usage per gene.
2. Run Rscript sum_codon_usage.R to summarize the results.
3. Rscript rscu.R, calculate codon usage for the highly expressed genes
4. Rscript CAI.R, compares the usage for each gene to the reference list
5. correlate_CAI_exp.R and correlateCAI_k.R performs correlation.
