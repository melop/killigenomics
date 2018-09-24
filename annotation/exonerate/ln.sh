#join the NFZ NOR assembled transcripts into a single file

 

sed 's/TRINITY_/NKH_T_/g' /beegfs/group_dv/home/RCui/killifish_genomes/annotation/transcripts/NKH/Trinity.fasta.transdecoder.mRNA > Trinity.fasta.transdecoder.mRNA
sed 's/TRINITY_/NKH_T_/g' /beegfs/group_dv/home/RCui/killifish_genomes/annotation/transcripts/NKH/Trinity.fasta.transdecoder.gff3 > Trinity.fasta.transdecoder.gff3

sed 's/TRINITY_/NFZ_T_/g' /beegfs/group_dv/home/RCui/killifish_genomes/annotation/transcripts/NFZ/Trinity.fasta.transdecoder.mRNA >> Trinity.fasta.transdecoder.mRNA
sed 's/TRINITY_/NFZ_T_/g' /beegfs/group_dv/home/RCui/killifish_genomes/annotation/transcripts/NFZ/Trinity.fasta.transdecoder.gff3 >> Trinity.fasta.transdecoder.gff3

sed 's/TRINITY_/NOR_T_/g' /beegfs/group_dv/home/RCui/killifish_genomes/annotation/transcripts/NOR/Trinity.fasta.transdecoder.mRNA >> Trinity.fasta.transdecoder.mRNA
sed 's/TRINITY_/NOR_T_/g' /beegfs/group_dv/home/RCui/killifish_genomes/annotation/transcripts/NOR/Trinity.fasta.transdecoder.gff3 >> Trinity.fasta.transdecoder.gff3


