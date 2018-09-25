#!/bin/bash

h2m=/software/source/hal/bin/hal2maf
maffilter=/software/source/MafFilter/maffilter
minlen=15;
arrSpecies=(  Xmac  PLP Alm  Fht ); 
ref=Nfur
for i in "${arrSpecies[@]}" 
do 
echo doing $i ...;
( $h2m progressiveCactusAln_5spp.hal ${i}2${ref}.maf --refGenome ${ref} --targetGenomes ${i} --onlyOrthologs --noAncestors --unique --noDupes ; \
php filtermafbyblocklength.php -i ${i}2${ref}.maf -o ${i}2${ref}.len${minlen}.maf -l $minlen; \
$maffilter \
input.file=${i}2${ref}.len${minlen}.maf \
input.file.compression=none \
output.log=${i}2${ref}.len${minlen}.maffilter.log \
maf.filter=" Subset( species=(${ref} , ${i}) , strict=yes, keep=no, remove_duplicates=yes) , AlnFilter( species=(${ref}, ${i}) ,  window.size=10, window.step=1, max.gap=10, max.ent=0.2, missing_as_gap=yes ) , MaskFilter(species=(${ref}) ,window.size=1, window.step=1, max.masked=0), SequenceStatistics( statistics=( SequenceLength( species=${i}), BlockLength, BlockCounts, PairwiseDivergence(species1=${ref}, species2=${i}) ) , ref_species=${ref}, file=${i}2${ref}.len${minlen}.statistics.csv) , MinBlockLength(min_length=3), VcfOutput( file=${i}2${ref}.len${minlen}.snp.vcf, compression=none, reference=${ref}, genotypes=(${i}) ), Output( file=${i}2${ref}.len${minlen}.filtered.maf, compression=none, mask=no)" > ${i}2${ref}.maffilter.screen ;\
maf-convert sam -f ../notfurScflds.dict ${i}2${ref}.len${minlen}.filtered.maf > ${i}2${ref}.len${minlen}.sam; \
sed "s/$ref\.Gap/Gap/" ${i}2${ref}.len${minlen}.sam > ${i}2${ref}.len${minlen}_fixedname.sam; \
samtools view -b -S ${i}2${ref}.len${minlen}_fixedname.sam > ${i}2${ref}.len${minlen}.bam; \
samtools sort -m 5000000000 ${i}2${ref}.len${minlen}.bam ${i}2${ref}.len${minlen}.sorted ; \
samtools index ${i}2${ref}.len${minlen}.sorted.bam;\
) &
done

wait

