echo " 61 2092416" > one_ind_per_sp.phy
grep -v "PLP[A-Z]" ../../phylipformat/full.phy | grep -v "AKN2" | grep -v  "^ " >> one_ind_per_sp.phy

