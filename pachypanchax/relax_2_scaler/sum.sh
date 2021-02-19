sp=Pachypanchax
srcdir=Relax_${sp}
out=sum_${sp}.txt
cat $srcdir/ret* | grep "GeneName" | head -n1 > $out
cat $srcdir/ret* | grep -v "GeneName" | grep "Success" >> $out
