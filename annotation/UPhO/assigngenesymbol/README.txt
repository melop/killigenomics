1. Link ../mclclusters/UPhO_orthogroups.csv to the current directory
2. Run make_genesymbol_table_ncbigff.php to make a map between mRNA id and gene symbol/gene names, for NCBI GFF3 format.
3. Run make_genesymbol_table_ensemblgff.php to make a map between mRNA id and gene symbol/gene names, for Ensembl GFF3 format.
4. Run assigngenesymbol.php to add gene symbol and gene names for UPhO_orthogroups.csv. Outputs UPhO_orthogroups.genesymbol.txt.
5. Run killi_orthologs_table.php to assign gene symbols to killifish genes, output killi_orthologs_table.ensemblinformed.txt

