sort -k1,1 -k4,4n ../plp.nfzref.acc.gff > plp.nfzref.acc.sorted.gff
bedtools intersect -a plp.nfzref.acc.sorted.gff -b annualconsv/cds.sorted.gff -sorted -wa -v > plp.nfzref.acc.nonCDS.gff
bedtools intersect -a plp.nfzref.acc.nonCDS.gff -b annualconsv/nfz.repeats.sorted.bed -sorted -wa -v > plp.nfzref.acc.nonCDS.norepeat.gff
