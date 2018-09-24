1. Use ln.sh to concatenate all transcoder processed transcripts into a a single .mRNA and .gff3 file
2. linebreakfasta.php is provided as a utility to reformat fasta files with the desired line width
3. Edit and run cdhitest.sh to cluster highly similar transcripts together, input Trinity.fasta.transdecoder.mRNA, output Trinity.cdhitout.fasta
4. Run filtertranscript.php to keep only the longest and complete transcript within each cluster. Output: Trinity.cdhitout.filtered.fa, Trinity.fasta.transdecoder.filtered.annot



