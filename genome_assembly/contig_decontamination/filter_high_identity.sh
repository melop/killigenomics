awk '{ if ($0 ~ /^[[:space:]]*#/) {} else { if ($3 >= 99) print $0 } } ' contig.megablast.txt 
