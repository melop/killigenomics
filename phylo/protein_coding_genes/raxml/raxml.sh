#!/bin/bash
#SBATCH -c 40
#SBATCH -p blade
ln -s `realpath ../aln/all.aln.fas` ./in.fas
raxmlHPC-PTHREADS-SSE3 -s in.fas -T 40 -f a -x 23333 -N 100 -m GTRGAMMA -n codon_all -p 23333
