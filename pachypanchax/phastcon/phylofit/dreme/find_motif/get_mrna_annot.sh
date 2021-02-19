cat /beegfs/group_dv/home/RCui/killifish_genomes/phastcon/cactus/NFZ2.0/NFZ2.0.gff | awk '{if ($3=="mRNA") print $0}' | sort -k1,1 -k4,4n > mRNA.gff
bedtools closest -a TAATTA.brief.gff -b mRNA.gff -D a -fd > TAATTA.brief.gff.mRNA.gff
bedtools closest -a CAGCTG.brief.gff -b mRNA.gff -D a -fd > CAGCTG.brief.gff.mRNA.gff
