cat ../../cactus/NFZ2.0/fixgff/cds.gff | awk '{if ($1=="NFZ.chr19") print $0}' > chr19_cds.gff
msa_view --in-format MAF --seqs NFZ.chr19  ../../cactus/NFZ_ref_5spp.maf --4d --features chr19_cds.gff > 4d-codons.ss
