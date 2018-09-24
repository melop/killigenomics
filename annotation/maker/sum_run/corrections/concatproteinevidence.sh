cat ../../proteins/uniprot-taxonomyAtherinomorphaeNOTCyprinodontiformes.fasta > all.protein.evidence.fa 

echo  >> all.protein.evidence.fa 

cat ../../proteins/uniprot-taxonomyCyprinodontiformes.fasta >> all.protein.evidence.fa

echo  >> all.protein.evidence.fa 

cat ../../proteins/uniprot-taxonomyTeleosteiNOTAtherinomorphae.fasta >> all.protein.evidence.fa

echo  >> all.protein.evidence.fa 

cat ../../proteins/NCBI_refseq_prot_teleostei_nomt.fasta >> all.protein.evidence.fa

fastahack -i all.protein.evidence.fa
