sDIR=smcpp_input_rawvcf_varqualmask
BLOCKSIZE=5e6
REPS=30
for i in ../$sDIR/*; do
	hhvm bootstrap.php $i $BLOCKSIZE $REPS $sDIR &
done
wait
