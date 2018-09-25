mkdir -p trees
cp ClusteRs/*_bipartitions.*.out ./trees/
cd ./trees

for i in *.out; do
  FILENAME=`basename $i` ;
  NEWNAME="${FILENAME%_clean*}.tre" ;
  NEWNAME=${NEWNAME#*_bipartitions.}
  mv $i $NEWNAME ;
done
