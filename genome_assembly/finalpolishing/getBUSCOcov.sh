BUSCO=./release/NOR/1.0/mergegap/busco/v2NOR_v1.0/run_v2NOR_v1.0/singlecopy_coord.bed #this is a bed file containing the coordinates of single copy buscos on the assembly

samtools view -b -L $BUSCO mapped/mergedPE.bam > mapped/mergedPE.sub.bam
bedtools coverage -a "$BUSCO" -b mapped/mergedPE.sub.bam -d > mapped/BUSCO.cov
