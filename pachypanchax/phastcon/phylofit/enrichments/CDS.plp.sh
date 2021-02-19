sort -k1,1 -k4,4n ../plp.nfzref.acc.gff > plp.nfzref.acc.sorted.gff
bedtools intersect -a plp.nfzref.acc.sorted.gff -b annualconsv/cds.sorted.gff -sorted -wa -wb > plp.nfzref.acc.CDS.gff
