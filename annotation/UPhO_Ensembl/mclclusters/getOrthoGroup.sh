MINTAXA=3
MINBOOSTRAP=0.50
PROT=allproteins.fa

rm UPhO_orthogroups.csv
rm -R UPhO_branches/
rm -R UPhO_Seqs/

UPhO.py -in trees/*.tre -m $MINTAXA -S $MINBOOSTRAP -ouT -iP -R ../$PROT > ortholog.log 2>&1;

