These scripts runs RELAX analyses given a set of gene alignments and labeled gene trees, excluding CMDs  (see ../export_protein_multiref).

1. Run any of the Relax_*.sbatch on SLURM to obtain results of the RELAX analysis.
	Output in Relax_XXXXXXXXX/ret*.txt, format: 

* GeneName
* Success or not
* GeneSymbol
* CodonCount
* P Value
* MG94xREV
* LogLk NULL
* LogLk Alternative
* K
* Likelihood Ratio
* Background omega0
* Background omega0_prop
* Background omega1
* Background omega1_prop
* Background omega2
* Background omega2_prop
* Foreground omega0
* Foreground omega0_prop
* Foreground omega1
* Foreground omega1_prop
* Foreground omega2
* Foreground omega2_prop
* GeneFullName
* tmpfolder

2. Edit and run summary.sh to summarize results into one file sum_xxxxx.txt
3. These results are not final, proceed to rerun_omega0 to rerun some of the difficult cases where the likelihood search may have stucked.

