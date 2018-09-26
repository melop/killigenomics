ln -s ../all.syn.aln.phy ./
/software/source/standard-RAxML-8.2.9/raxmlHPC-PTHREADS-SSE3 -s all.syn.aln.phy  -f e -T 40  -m GTRGAMMA -n syn -t ../raxml/RAxML_bestTree.allcodons
