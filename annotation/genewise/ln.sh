ln -sf `realpath ../prot_nogap` ./ #link the probe protein sequences
ln -sf `realpath ../hmm` ./ #link the hmm for each gene (hmmer2 format)



ln -sf ./repeatmasker/NOR_v1.0_pilon/scf.fa.masked ./genome.masked.fa #Link the masked version of the genome, for blast
ln -sf ./nor.pilon.rename.fa ./genome.fa #link the unmasked version of the genome

#index

makeblastdb -in genome.masked.fa -dbtype nucl
fastahack -i genome.fa
