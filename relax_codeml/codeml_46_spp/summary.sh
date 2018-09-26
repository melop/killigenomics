#clade=Nothobranchius
#clade=Callopanchax
#clade=Aphyosemion
#clade=Scriptaphyosemion
#clade=Archiaphyosemion
#clade=Epiplatys
#clade=Fundulopanchax
#clade=PronothoNothos
#clade=AKN
clade=Nothobranchius_mrca
clade=Aphyosemion_mrca
clade=Callopanchax_mrca
#clade=Scriptaphyosemion_mrca
cat Codeml_${clade}/ret*.txt | grep -P "Success\t" | sort -t $'\t' -k 5n,5n  > sum_${clade}.txt
cat Codeml_${clade}/ret*.txt | grep -P -v  "Success\t" | sort -t $'\t' -k 1,1  > sum_${clade}_err.txt

echo $clade `wc -l sum_${clade}.txt`
echo $clade error `wc -l sum_${clade}_err.txt`
