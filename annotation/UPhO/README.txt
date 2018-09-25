1. Use prepareProteins.sh to preprocess the killifish proteins for UPhO. This outputs several *.fst files, and a concatenated file allkillifishproteins.fa
2. Use the same method to preprocess Ensembl proteins. Concatenate with the above file into allproteins.fa using joinproteins.sh.
3. Run UPhO.sbatch to perform all-against-all blastp
4. Proceed to the folder mclclusters/

