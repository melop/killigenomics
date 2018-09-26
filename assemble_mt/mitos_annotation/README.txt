1. mt/ contains the assembled mt genomes. 
2. Run submitall.sh to submit mt genomes for MITOS annotation, output in out/. Run checkall.sh periodically to see if the process is done.
3. Run adjust_contig_start.sh to adjust the mt genomes to the start of 12s rDNA, output in mt_pos_adjusted/. Manually check the annotation, correct small assembly errors such as frame shifts that break annotations.
4. Run submitall_adjust.sh to submit the adjusted mt genomes for MITOS annotation, output in out_adjusted/. Run checkall_adjusted.sh periodically to see if the process is done.
5. Proceed to alignments/ to perform further analyses.
