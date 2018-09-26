1. Run Relax_xxxxxx.sbatch on SLURM to redo some of the problematic cases excluding CMDs. The results from the parent directory will be automatically read, and if it is identified as having convergence issues, the test is rerun using a newer version of hyphy with multiple CPUs for at most 20 times.
2. The results are updated and written into Relax_xxxxxx/
3. Edit and run summary.sh
4. Run Rscript make_pretty_table.R to make the results more readable. 
