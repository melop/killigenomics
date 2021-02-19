FILE=../../../fortreemix.txt.gz
treemix=/beegfs/group_dv/software/source/treemix/src/treemix

$treemix -i $FILE -k 1000 -root PML -o mltree > log.txt 2>&1
