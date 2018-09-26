These scripts perform CodeML tests, without exclusion of CMDs.

1. Run Codeml_*.sbatch on SLURM to do codeml tests for all genes.
2. Edit and run summary.sh to get sum_xxxxx.txt tables.
Output format:
* GeneName
* Success or not
* GeneSymbol
* Codon Length
* Log Likelihood ratio
* Likelihood MA
* Likelihood MA_null
* Temporary folder
* GeneFullName
3. The utility exclude_positive_sites.php allows masking of the sites inferred to be positively selected.
4. make_pretty_table.R produces better-looking tables
