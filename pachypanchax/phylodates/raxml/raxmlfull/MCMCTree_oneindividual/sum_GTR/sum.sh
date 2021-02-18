export PATH=/beegfs/group_dv/software/source/paml4.9i/src/:$PATH

head -n1 ../gtr_approx01/mcmc.txt > mcmc.txt
tail -n400000 ../gtr_approx01/mcmc.txt >> mcmc.txt 
tail -n400000 ../gtr_approx02/mcmc.txt >> mcmc.txt

mcmctree mcmctree.ctl > screenlog.log 2>&1
