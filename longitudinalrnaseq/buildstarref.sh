STAR=/beegfs/group_dv/software/source/STAR-2.6.0c/bin/Linux_x86_64_static/STAR

mkdir -p NORref

$STAR --runThreadN 20 \
--runMode genomeGenerate \
--genomeDir NORref \
--genomeFastaFiles ref.fa
