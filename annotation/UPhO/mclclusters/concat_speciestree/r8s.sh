TREE=rooted.tre
nAlnLen=`head -n 1 alnlen.txt  | grep -oP "\d+"`
python /beegfs/group_dv/software/source/CAFE/cafe_tutorial/python_scripts/cafetutorial_prep_r8s.py -i $TREE -o r8s.conf.txt -s $nAlnLen \
-p 'spotted_gar,zebrafish zebrafish,gadus gadus,takifugu stickleback,tilapia stickleback,takifugu tilapia,medaka medaka,xmac' \
-c '361.2,266.1,139.7,99.75,89.4,86.1,75.38'

#-p 'spotted_gar,zebrafish zebrafish,gadus gadus,takifugu stickleback,tilapia stickleback,takifugu tilapia,medaka medaka,xmac PLP,xmac' \
#-c '361.2,266.1,139.7,99.75,89.4,86.1,75.38,40.5'

r8s -b -f r8s.conf.txt > r8s_tmp.txt
tail -n 1 r8s_tmp.txt | cut -c 16- > ultrametric.tre
