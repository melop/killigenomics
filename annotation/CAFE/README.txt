1. Run excludeoverlapgenes.php to exclude gene counts not included in the CAFE white list (produced by ../maker/sum_run/corrections/whitelistgenes.php), producing out.seq.mci.I20
2. Prepare an ultrametric tree for the taxa, save as ultrametric.tre.
3. Use cafe.sh to estimate errors. Output in reports/errest_sp/
4. Modify cafefinal.sh to use the estimated error profiles estimated in the previous step for inferring gene family expansion and contractions. Output in reports/errest_sp_allgenes.
5. Check GOAnalysis folder if you want to export gene lists for GO analysis.
6. Check likelihoodRatio folder if you want to perform test of gene birth/death rates
