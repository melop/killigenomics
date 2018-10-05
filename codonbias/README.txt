Map RNAseq reads to the reference genome with STAR

1. Link files with ln.sh (example)
2. Run makeref.sh to prepare the reference genome for mapping.
3. Run mapall.sh to perform mapping.
4. Proceed to mappedPE/featurecounts/ to perform read counts
5. Proceed to gene_expression_levels/ to get a list of most highly expressed genes. Also contains scripts for correlating expression and relax k.
6. Proceed to codon_usage to compute codon adaptive index
