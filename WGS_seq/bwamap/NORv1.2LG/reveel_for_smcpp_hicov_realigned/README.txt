This folder contains scripts used to prepare input files for smc++.

1. Use ln_bams.sh to link bam files (mapped to the reference genome).
2. Run makevcf.sh to call variants for each population. Run makevcf_twopops.sh to jointly call variants for sister populations.
3. Run makemasks.varqual.sh and makemasks.2pops.varqual.sh to produce quality masks.
4. Run raw2smc_varqual_mask.sh to convert the vcfs and the masks to smc++ input.raw2smc_2pops_*.sh can be used to convert the two-population vcfs, but later analysis shows too much RAM requirement. Thus the reveel imputation is used instead.
5. Run Reveel to impute genotypes with reveel.sh.
6. Run reveel2smc_2pops_ort.sh to convert the reveel output to smc++ format. 
7. Perform bootstrapping with the scripts in bootstrap/ for single populations, and bootstrap_2pops/ for paired populations
