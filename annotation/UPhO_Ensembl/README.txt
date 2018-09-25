This folder contains scripts used to prepare HMM profiles for Genewise alignments.
For the scripts used for clustering Killifish genes, check ../UPhO/
1. Use prepareProteins.sh to preprocess the Ensembl proteins for UPhO. This outputs several *.fst files, and a concatenated file allproteins.fa
2. Run UPhO.sbatch to perform all-against-all blastp, use mclcluster.sh to convert blast results to MCL input.
3. Proceed to the folder mclclusters/

