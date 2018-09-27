
sPop1=RACDRY.excludehybrid
sPop2=ORTWET

sOutDIR=${sPop1}_${sPop2}
mkdir $sOutDIR
hhvm computeFST.php $sPop1 $sPop2 $sOutDIR > $sOutDIR/fst.log 2>&1 

