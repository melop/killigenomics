for j in errmodels/cafe_errormodel_*.txt; do
	sed -i 's/maxcnt:20000/maxcnt:1000/' $j

done
