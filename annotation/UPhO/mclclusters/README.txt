1. Perform MCL clustering and output fasta sequences for each gene cluster into ClusteRs/. 
2. evaluate_cluster.sh outputs basic statistics for the run.
3. Run aln_orthogroups.sh to perform sequence alignment and cleaning for each gene cluster.
4. Run UPhO_mcl_raxml.sbatch on SLURM to get gene trees for each gene cluster.
5. Run UPhO_mcl_raxml_fast.sbatch on SLURM (or fasttree_difficult.sh on single computer) for large clusters using fasttree. 
6. Use cptrees.sh to copy trees to trees/
7. Run getOrthoGroup.sh to obtain orthologous groups
8. Run getexactorthologs.php to obtain orthologous groups with a complete representation of all fish species
9. Optionally HMM profiles can be built for each gene family with makehmm.sh
