grep -P "\tgene\t" finalrun_out3.gff > genecounts.gff
wc -l genecounts.gff
bedtools merge -s -i genecounts.gff > genecounts.merged.gff
wc -l genecounts.merged.gff
