arrFolders=(`ls out`)
outfile=final.gtf

echo "" > $outfile
for f in "${arrFolders[@]}"
do 
  echo "Processing $f ...";
 # php pickbestmodel2.php -f out/$f/out.gff -o out/$f/picked.gff -k 1;
  cat out/$f/picked.gff >> $outfile 2>sumret.err.log;
done

