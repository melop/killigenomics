GENOME=/beegfs/group_dv/home/RCui/killifish_genomes/annotation/repeatmasker/NOR_v1.0_pilon/scf.fa

mkdir -p genome
ln -sf `realpath $GENOME` ./genome/genome_pilon.fa

mkdir ESTs #link conspecifc or congeneric ESTs
ln -s /beegfs/group_dv/home/RCui/killifish_genomes/annotation/transcripts/NOR/Trinity.fasta.transdecoder.mRNA ./ESTs/NOR.trinity.mRNA
ln -s /beegfs/group_dv/home/RCui/killifish_genomes/annotation/transcripts/NFZ/Trinity.fasta.transdecoder.mRNA ./ESTs/NFZ.trinity.mRNA
ln -s /beegfs/group_dv/home/RCui/killifish_genomes/annotation/transcripts/NKH/Trinity.fasta.transdecoder.mRNA ./ESTs/NKH.trinity.mRNA

mkdir proteins
ln -s /beegfs/group_dv/home/RCui/killifish_genomes/annotation/transcripts/CTO/Trinity.fasta.transdecoder.pep proteins/CTO.trinity.pep
ln -s /beegfs/group_dv/home/RCui/killifish_genomes/annotation/transcripts/A.striatum/Trinity.fasta.transdecoder.pep proteins/AST.trinity.pep
ln -s /beegfs/group_dv/home/RCui/killifish_genomes/annotation/transcripts/AAU/Trinity.fasta.transdecoder.pep proteins/AAU.trinity.pep
ln -s /beegfs/group_dv/home/RCui/killifish_genomes/annotation/transcripts/NFZ/Trinity.fasta.transdecoder.pep proteins/NFZ.trinity.pep
ln -s /beegfs/group_dv/home/RCui/killifish_genomes/annotation/transcripts/NKH/Trinity.fasta.transdecoder.pep proteins/NKH.trinity.pep
ln -s /beegfs/group_dv/home/RCui/killifish_genomes/annotation/transcripts/NOR/Trinity.fasta.transdecoder.pep proteins/NOR.trinity.pep
ln -s /beegfs/group_dv/home/RCui/killifish_genomes/annotation/transcripts/PLP/Trinity.fasta.transdecoder.pep proteins/PLP.trinity.pep
ln -s /beegfs/group_dv/home/RCui/killifish_genomes/annotation/uniprot/NCBI_refseq_prot_teleostei_nomt.fasta proteins/


ln -s `realpath ../../uniprot/uniprot-taxonomyAtherinomorphaeNOTCyprinodontiformes.fasta` proteins/
ln -s `realpath ../../uniprot/uniprot-taxonomyCyprinodontiformes.fasta` proteins/
ln -s `realpath ../../uniprot/uniprot-taxonomyTeleosteiNOTAtherinomorphae.fasta` proteins/
