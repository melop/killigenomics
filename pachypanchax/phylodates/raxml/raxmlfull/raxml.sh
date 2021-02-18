ln -s ../phylipformat/full.phy ./in.phy
/beegfs/group_dv/software/source/standard-RAxML/raxmlHPC-PTHREADS-SSE3 -s in.phy -T 40 -f a -x 23333 -N 100 -m GTRGAMMA -n allcodons -p 23333

