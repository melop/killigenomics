#!/beegfs/group_dv/software/source/CAFE/release/cafe
load -i ../filtered_cafe_input.txt -t 1 -l reports/1rate.log.txt
tree ((((((((xmac:14.621276,poecilia:14.621276):52.754098,(((NOR:22.391744,AAU:22.391744):21.669600,CTO:44.061344):11.674051,PLP:55.735395):11.639979):8.004626,medaka:75.380000):10.720000,tilapia:86.100000):13.650000,((tetraodon:42.699828,takifugu:42.699828):46.700172,stickleback:89.400000):10.350000):39.950000,gadus:139.700000):126.400000,(zebrafish:163.544139,astyanax:163.544139):102.555861):95.100000,spottedgar:361.200000)
errormodel -model errmodels/cafe_errormodel_0.0206884765625.txt -sp gadus
errormodel -model errmodels/cafe_errormodel_0.0724096679687.txt -sp takifugu
errormodel -model errmodels/cafe_errormodel_0.0930981445312.txt -sp PLP
errormodel -model errmodels/cafe_errormodel_0.08275390625.txt -sp tilapia
errormodel -model errmodels/cafe_errormodel_0.124130859375.txt -sp tetraodon
errormodel -model errmodels/cafe_errormodel_0.0.txt -sp astyanax
errormodel -model errmodels/cafe_errormodel_0.103442382812.txt -sp xmac
errormodel -model errmodels/cafe_errormodel_0.155163574219.txt -sp medaka
errormodel -model errmodels/cafe_errormodel_0.134475097656.txt -sp AAU
errormodel -model errmodels/cafe_errormodel_0.124130859375.txt -sp poecilia
errormodel -model errmodels/cafe_errormodel_0.0.txt -sp zebrafish
errormodel -model errmodels/cafe_errormodel_0.103442382812.txt -sp CTO
errormodel -model errmodels/cafe_errormodel_0.0.txt -sp spottedgar
errormodel -model errmodels/cafe_errormodel_0.186196289062.txt -sp NOR
errormodel -model errmodels/cafe_errormodel_0.0620654296875.txt -sp stickleback
lambdamu -l 0.00101637793135 -m 0.00082194099750 -t ((((((((1,1)1,(((1,1)1,1)1,1)1)1,1)1,1)1,((1,1)1,1)1)1,1)1,(1,1)1)1,1)
genfamily  genfamily_sim/rnd -t 1000

