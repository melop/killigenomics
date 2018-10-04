1. Run match_sample_names.php to link the simulated reads into file name formats that match the real data.
2. Run map.local.sbatch to use BWA to map the simulated reads
3. Use makevcf.sh to call variants, followed by reveel.sh to impute genotypes.

