MSMC2=/beegfs/group_dv/software/source/msmc2/build/release/msmc2_linux64bit

arrPops=( PLP_A_PLP-A3_D501_D701 
 PLP_A_PLP-A4_D502_D702  
 PLP_B_PLP-B3_D504_D704 
 PLP_C_PLP-C1_D506_D706
 PLP_C_PLP-C2_D507_D707
 PLP_C_PLP-C3_D508_D708
 PLP_C_PLP-C4_D502_D701
 PLP_D_PLP-D0_ref
 PLP_D_PLP-D1_D503_D702
 PLP_E_PLP-E1_D505_D704
 PLP_E_PLP-E2_D506_D705
 PLP_E_PLP-E3_D507_D706
 PLP_E_PLP-E4_D508_D707
 PLP_E_PLP-E5_D501_D708
 PLP_F_PLP-F1_D503_D701
 PLP_F_PLP-F2_D504_D702
 PLP_F_PLP-F3_D505_D703
 PLP_G_PLP-G1_D506_D704
 PLP_G_PLP-G3_D508_D706
 PLP_G_PLP-G5_D502_D708 )
sIn=formsmc2_in
sOut=msmc2ret
mkdir -p $sOut


for sPop in "${arrPops[@]}"; do
	sFiles="${sIn}/${sPop}*/chr*/formsmc2.multihetsep.txt"
	sCmd="$MSMC2 -o $sOut/$sPop -i 40 -t 2 -r 2.61976354 $sFiles > /dev/null "
	eval $sCmd &
done

wait
