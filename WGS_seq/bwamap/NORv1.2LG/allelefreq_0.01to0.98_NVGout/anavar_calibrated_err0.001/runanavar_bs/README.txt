1. Run makebsreplicates.R to produce multinomial resampling of the original SFS, write to bs/
2. Run *.sbatch on SLURM to run Anavar on every bootstrap replicate
3. Edit sumresults_bs.R to adjust the range of the DFE to visualize, then run it to produce a figure.
