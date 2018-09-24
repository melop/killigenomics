grep -P "\tmRNA\t" split__candidate.gff > candidate.mRNA.gff
grep -P "\tmRNA\t" maker.refined.removeisoform.gff | grep -v "trna" > currentmakermodel.mRNA.gff

bedtools intersect -a candidate.mRNA.gff -b currentmakermodel.mRNA.gff -v -s > nonoverlap_candidate.gff 
