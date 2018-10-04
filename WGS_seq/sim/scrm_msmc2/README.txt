This folder contains simulation scripts using scrm.

1. Run ORT.sh and RAC.sh to simulate genomes given the demography.
2. Run wgsim_lowcov_err*.sh to simulate reads at low coverage at the given sequencing error rates.
3. Run wgsim_hicov_err*.sh to simulate reads at high coverage for a single individual from each population at given error rates.
4. Proceed to bwaalign/ to analyze the performance of PSMC' using the simulated reads.
5. Proceed to msmc2/ to analyze the performance of PSMC' using the true simulated genotypes.
