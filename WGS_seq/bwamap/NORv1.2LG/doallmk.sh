sSrc=`pwd`/allelefreq_0.02to0.1_NVGout/MKTEST_syn
WD=`pwd`

for i in `seq 1 8`; do 
    cd allelefreq_0.${i}to0.$(( i + 1 ))_NVGout/MKTEST_syn; 
    cp $sSrc/mk_alpha_per_gene.sh ./;
    cp $sSrc/mk_alpha_per_gene.R ./;
    cp $sSrc/combime_pop_mk.R ./;
    cp $sSrc/correlate_alpha_k.R ./;
    cp $sSrc/correlate_codeml_combined.R ./;
    bash mk_alpha_per_gene.sh  0.${i}to0.$(( i + 1 )) &
    cd $WD ; 

done

wait
