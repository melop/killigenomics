 grep "Average Haploid LOD:" slurm-*.out | sed -r 's/.+Average Haploid LOD: //'  > HaploidLODdistr.txt

