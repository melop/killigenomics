mkdir -p prot_nogap
mkdir -p hmm

#Link the ensembl hmms and proteins
for i in ../UPhO_Ensembl/mclclusters/UPhO_Seqs/*.fasta ; do
        ln -sf $i  ./prot_nogap/
done

for j in ../UPhO_Ensembl/mclclusters/hmm/*.hmm; do
        ln -sf $j ./hmm/
done


ln -sf ./repeatmasker/NOR_v1.0_pilon/scf.fa.masked ./genome.masked.fa #Link the masked version of the genome, for blast
ln -sf ./nor.pilon.rename.fa ./genome.fa #link the unmasked version of the genome

#index

makeblastdb -in genome.masked.fa -dbtype nucl
fastahack -i genome.fa
