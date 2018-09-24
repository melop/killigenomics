arrFolders=(`ls out`)
outfile=finalrun_out3.gtf

echo "" > $outfile
for f in "${arrFolders[@]}"
do 
  echo "Processing $f ...";
  cat out/$f/out.gtf >> $outfile;
done

