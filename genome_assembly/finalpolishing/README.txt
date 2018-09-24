1. Use mapreads.sh to map pair-end reads to the assembly
2. Run BUSCO on the assembly to be fixed if you haven't done so, convert the coordinates of the single-copy BUSCOs into a bed file
3. Run getBUSCOcov.sh to get coverage at the single-copy BUSCOs, to be used as an estimate of expected coverage, change the nAvgGenomeCov variable in blastmerge_cov.sh
4. run blastmerge_cov.sh to output into the temporary directory fasta and agp files. Redirect the output to a log file. This step can be parallelized on SLURM
5. Summarize the run by summarize.sh, this should output the modified fasta, the agp file, and a lifted version of the masked fasta file.
6. You can use the command in getHaploidLODdistr.sh to summarize the LOD scores (between haploid versus diploid coverage at the overhang part), and use plotDupHaploidLOD.R to visualize.
7. stats.sh contains a command to run quast to obtain the new N50 statistics.


