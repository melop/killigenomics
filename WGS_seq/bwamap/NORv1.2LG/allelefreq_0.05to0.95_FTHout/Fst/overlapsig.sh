nPcutoff=0.01;


nWindowSize=50000;
nStepSize=10000;
nMinSites=10;


bedtools intersect -a RAC.slideFST.win${nWindowSize}.step${nStepSize}.minsite${nMinSites}.sig${nPcutoff}.bed -b ORT.slideFST.win${nWindowSize}.step${nStepSize}.minsite${nMinSites}.sig${nPcutoff}.bed -wa -wb > overlap.win${nWindowSize}.step${nStepSize}.minsite${nMinSites}.sig${nPcutoff}.bed

grep -P "maker\tgene\t" /beegfs/group_dv/home/RCui/killifish_genomes/denovo/discovardenovo/NOR/LG/v1.2/NORv1.2.LG.gff > NOR.gene.coord.gff
bedtools intersect -b overlap.win${nWindowSize}.step${nStepSize}.minsite${nMinSites}.sig${nPcutoff}.bed -a NOR.gene.coord.gff -wa  | sort | uniq > overlap.genes.win${nWindowSize}.step${nStepSize}.minsite${nMinSites}.sig${nPcutoff}.bed
