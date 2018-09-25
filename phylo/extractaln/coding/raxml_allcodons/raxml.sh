#!/bin/bash
#SBATCH -c 40
#SBATCH -p blade
ln -s ../ret_allcodons_table.txt.fasta ./in.fas
raxmlHPC-PTHREADS-SSE3 -s in.fas -T 40 -f a -x 23333 -N 100 -m GTRGAMMA -n allcodons -p 23333
