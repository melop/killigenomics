grep -o -P "_AED=([\d\.]+)" genome.all.gff | grep -o -P "[\d\.]+" > AEDscores.txt

