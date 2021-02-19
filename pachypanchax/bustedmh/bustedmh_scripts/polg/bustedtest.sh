ln -sf /beegfs/group_dv/home/RCui/killifish_genomes/test_relax_models/polg.nex ./
BUSTEDMH="/beegfs/group_dv/software/source/hyphy-develop_install/bin/hyphy /beegfs/group_dv/home/RCui/killifish_genomes/test_bustedmh/BUSTED-MH.bf"

$BUSTEDMH --alignment polg.nex --branches T > screen.log 2>&1 

