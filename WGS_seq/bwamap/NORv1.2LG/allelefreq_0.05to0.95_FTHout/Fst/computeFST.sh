mkdir ORT
mkdir RAC

hhvm computeFST.php ORTDRY.excludehybrid ORTWET ORT > ORT/fst.log 2>&1 &
hhvm computeFST.php RACDRY.excludehybrid RACWET.excludehybrid RAC > RAC/fst.log 2>&1 &
wait
