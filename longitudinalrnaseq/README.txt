A reanalysis of a previously published RNAseq dataset.

1. Go to download/ to automatically load data from SRA. The cited datasets can be seen in download/datalist.txt.
2. Run buildstarref.sh to prepare the reference genome for STAR.
3. Run mapall.sh to map all reads.
4. Perform featureCount on each sample using countall.sh.
5. Proceed to differential expression analysis using DEseq2/

