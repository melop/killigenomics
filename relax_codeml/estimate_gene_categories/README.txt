1. After performing simualtions using ../simulations/ and running the simulated results with codeml and relax like the real data, edit the paths in collect_sim_ret.php.
2. Run collect_sim_ret.php to output a table of test performance. Outputs sum_ret_XXXXXXX.txt and real_ret_XXXXXXXX.txt
3. Edit and run solve_for_gene_numbers.ML.R to get an estimate of the different gene categories in the real dataset based on simulated performance.
4. compare_between_clades.R allows comparing between gene numbers in different clades.
5. excludeCMD/ is identical with this folder, only with real data sourced from the ../relax_46_spp_excludeCMD and the ../codeml_46_spp_excludeCMD folders.

