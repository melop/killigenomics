1. This pipeline assumes that you have mapped the reads to the reference genome (in LG) and made the bam files
2. Run countdepth.sh to calculate the average sequencing depth for each sample. Write down the depth in bamcaller.sh.
3. Run bamcaller.sh to prepare input for msmc2.
4. Run runmsmc2.sh to run PSMC' on the input files.
5. Run bootstrap.sh to produce the bootstrap pseudoreplicates.
6. Run bootstrap_run.sh to analyze the bootstrap replicates
7. Plot the results with the bootstraps with plotting.R
8. The utility script smc2ms_command.sh converts the PSMC' output into ms style commands for simulation.
9. Optional: the smc2ms_command_recent_crash.sh script also adds a recent bottleneck. 
