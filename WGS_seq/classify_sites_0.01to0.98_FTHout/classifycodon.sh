hhvm classifycodon.php -p RACDRY.excludehybrid > RACDRY.log 2>&1 &
hhvm classifycodon.php -p RACWET.excludehybrid > RACWET.log 2>&1 &
hhvm classifycodon.php -p ORTDRY.excludehybrid > ORTDRY.log 2>&1 &
hhvm classifycodon.php -p ORTWET > ORTWET.log 2>&1 &
wait
