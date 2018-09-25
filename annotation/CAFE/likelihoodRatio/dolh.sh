mkdir reports
mkdir compute_sims
mkdir slurmlogs

./cafe1ratelambdamu.sh > 1rate.screen.log 2>&1 & #do 1 rate

./cafe2ratelambdamu.sh > 2rate.screen.log 2>&1 & #do 2 rate

wait
