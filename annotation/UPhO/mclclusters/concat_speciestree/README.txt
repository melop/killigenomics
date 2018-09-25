1. Run concataln.php to concatenate all amino acid sequences in strict 1-to-1 orthologs represented in all 16 fish species into concat.aln.phy
2. Run raxml with raxmlHPC-PTHREADS-SSE3 -T 40 -f a -p 688 -x 23333 -N 100 -m PROTGAMMAJTTX -s concat.aln.phy -n concattree 
3. Root the tree with FigTrees or a similar tool, save as rooted.tre
4. Edit r8s.conf.txt, and use r8s.sh to obtain an ultrametric tree with dates. This tree was used for CAFE' gene family analyses.
