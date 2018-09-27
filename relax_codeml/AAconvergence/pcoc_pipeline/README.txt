1. Run the pcoc pipeline with pcoc.sbatch, run excludeCMD.R to exclude CMD sites.
2. Run the native pcoc simulations with pcoc_sim.sbatch
3. Run the sim with the ../simulations/ND alignments, run pcoc_codeml_sim.sbatch
4. Use sum_compare.R and sum_compare_excludeCMD.R to check for enrichments of convergent AAs in the relaxed and positive genes.
