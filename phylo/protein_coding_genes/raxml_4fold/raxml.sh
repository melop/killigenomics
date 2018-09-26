#!/bin/bash
#SBATCH -c 40
#SBATCH -p blade
ln -s `realpath ../aln/all.4fold.aln.phy` ./in.fas
raxmlHPC-PTHREADS-SSE3 -s in.fas -T 20 -f a -x 23333 -N 100 -m GTRGAMMA -n codon_4fold -p 23333
