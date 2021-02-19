FILE=../../fortreemix.txt.gz
treemix=/beegfs/group_dv/software/source/treemix/src/treemix
Migrate=1

$treemix -i $FILE -k 1000 -root G -m $Migrate -tf ../ML/reroot.txt -se  -o migrate_${Migrate} > log.txt 2>&1

