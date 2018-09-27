ln -s ../denovo/discovardenovo/NOR/LG/v1.2/Both_Crosses.fasta ./ref.fa
bwa index ref.fa
samtools1.2 faidx ref.fa
PICARDDIR="/beegfs/group_dv/software/source/picard-tools-1.119/"
java -Dsnappy.disable=true -Xmx12g -jar $PICARDDIR/CreateSequenceDictionary.jar R=ref.fa O=ref.dict 
