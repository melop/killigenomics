cat ../../chr*_cds.gff | awk '{print substr($0, 5) }' > cds.gff
sort -k1,1 -k4,4n cds.gff > cds.sorted.gff
sort -k1,1 -k4,4n ../../annual.mostconserved.gff > annual.mostconserved.sorted.gff
bedtools intersect -a annual.mostconserved.sorted.gff -b cds.sorted.gff -sorted -wa -v > annual.consv.nonCDS.gff
sort -k1,1 -k2,2n nfz.repeats.bed > nfz.repeats.sorted.bed
bedtools subtract -a annual.consv.nonCDS.gff -b nfz.repeats.sorted.bed -sorted > annual.consv.nonCDS.norepeat.gff
