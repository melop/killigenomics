ln -sf /beegfs/group_dv/home/RCui/killifish_genomes/WGS_seq/sim/ray_scrm2/chr1.withN.fa ./ref.fa
bwa index ref.fa
samtools1.2 faidx ref.fa
