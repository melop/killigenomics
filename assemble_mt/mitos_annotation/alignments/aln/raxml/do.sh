ln -s ../all.aln.fas ./
raxmlHPC-PTHREADS-SSE3 -s all.aln.fas -q partitions.txt -f a -T 40 -x 23333 -N 100 -m GTRGAMMA -n allcodons -p 23333 
