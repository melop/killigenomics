SP=ORT
POPNAME1=${SP}WET
POPNAME2=${SP}DRY
POP1CHR=120
POP2CHR=118
#POP1CHR=2
#POP2CHR=2

#ANCFASTA=anc.fa
ANCFASTA=chr1.NOR.replaced.fa #chr1.small.fa #chr1.NOR.replaced.fa
theta0=0.0056522357187345 # 0.000000006889976117 * 4 * 2.63e-9
rho0=0.01804281451 # 654676.2357 * 4 * 0.000000006889976117
recfile=recmap.txt

CHRSIZE=`grep -v ">" $ANCFASTA | tr -d  "\n" | wc -c`

hhvm simpop.php $POPNAME1 $POPNAME2 $POP1CHR $POP2CHR $theta0 $rho0 $ANCFASTA ' -en  0  2  0.434270805 -ej 0.05393647 2 1' sim${SP}.out.txt gtr4fold.txt $recfile 2>&1 > ${SP}.sim.log

perl simulate_killifishdemo_multi_provideseq.pl sim${SP}.out.txt $POPNAME1 $(( POP1CHR/2 )) 2.2306 0.43664 $CHRSIZE 2>&1 > ${SP}WET.simread.log &

perl simulate_killifishdemo_multi_provideseq.pl sim${SP}.out.txt $POPNAME2 $(( POP2CHR/2 )) 2.2356 0.48717 $CHRSIZE 2>&1 > ${SP}DRY.simread.log &

wait
