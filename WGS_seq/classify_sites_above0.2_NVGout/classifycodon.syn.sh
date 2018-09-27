hhvm classifycodon.syn.php -p RACDRY.excludehybrid > RACDRY.excludehybrid.syn.classify.log 2>&1 &
hhvm classifycodon.syn.php -p RACWET.excludehybrid > RACWET.excludehybrid.syn.classify.log 2>&1 &
hhvm classifycodon.syn.php -p ORTDRY.excludehybrid > ORTDRY.excludehybrid.syn.classify.log 2>&1 &
hhvm classifycodon.syn.php -p ORTWET > ORTWET.syn.classify.log 2>&1 &
wait
