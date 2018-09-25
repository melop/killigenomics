grep -vP "loc\d+" UPhO_orthogroups.genesymbol.txt | grep -v "unknown" | grep -v "si:" | grep -v "wu:" | wc -l
