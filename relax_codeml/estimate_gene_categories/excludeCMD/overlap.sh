 cat codeml/*_*/* |  cut -f1,1 | sort > pos_all.txt
 cat codeml/*_*/* |  awk '{if (  $5 <0.05 ) print $0; }' | cut -f1,1 | sort > pos.txt

 cat relax/*_*/* |  cut -f1,1 | sort > relax_all.txt
 cat relax/*_*/* |  awk '{if (  $5 <0.05 && $9 < 1 ) print $0; }' | cut -f1,1 | sort > relax.txt

 echo overlap
 comm -12 pos.txt relax.txt | wc -l
 echo outof
 comm -12 pos_all.txt relax_all.txt | wc -l

