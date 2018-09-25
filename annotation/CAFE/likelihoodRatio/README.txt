1. Copy error models into the errmodels folder from the parent folder.
2. Simulate data with cafe_simulation_lambdamu.sh. Simulated datasets are in genfamily_sim/.
3. Fit models by running computesims.sbatch on slurm.
4. Run likelihood ratio tests by dolh.sh. 
5. Summarize likelihood ratio tests by summarize_sims.php. Output sim_ln_diff.txt.
6. The distribution of logLK ratio can be compared to the observed value to obtain a p value.
