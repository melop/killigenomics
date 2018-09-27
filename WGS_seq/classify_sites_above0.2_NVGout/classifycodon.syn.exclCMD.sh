hhvm classifycodon.syn.exclCMD.php -p RACDRY.excludehybrid > RACDRY.excludehybrid.syn.exclCMD.classify.log 2>&1 &
hhvm classifycodon.syn.exclCMD.php -p RACWET.excludehybrid > RACWET.excludehybrid.syn.exclCMD.classify.log 2>&1 &
hhvm classifycodon.syn.exclCMD.php -p ORTDRY.excludehybrid > ORTDRY.excludehybrid.syn.exclCMD.classify.log 2>&1 &
hhvm classifycodon.syn.exclCMD.php -p ORTWET > ORTWET.syn.exclCMD.classify.log 2>&1 &
wait
