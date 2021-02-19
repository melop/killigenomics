setwd("/beegfs/group_dv/home/RCui/killifish_genomes/plp_treemix/treemix/optm");
source("/beegfs/group_dv/software/source/treemix/src/plotting_funcs.R")

pdf(file = "test.2.0.pdf", width = 6, height=5);
plot_tree("test.2.0")
dev.off();

pdf(file = "test.2.1.pdf", width = 6, height=5);
plot_tree("test.2.1")
dev.off();

pdf(file = "test.2.2.pdf", width = 6, height=5);
plot_tree("test.2.2")
dev.off();

#0.389779 0.389779 0 <2.22507e-308 D.LaDigue:0.00178441 (C.PraslinWest:0.000668725,B.PraslinPlaineHollandaise:0.000772664):0.0116698
#0.412078 0.411749 0 <2.22507e-308 D.LaDigue:0.0230747 (C.PraslinWest:0.00684906,B.PraslinPlaineHollandaise:0.00724074):0.120865

#plot_resid(stem = "migrate_1", pop_order = "poporder.txt" )

# The tree inferred from the data is in outstem.treeout.gz. The first line of this file is the
# Newick format ML tree, and the remaining lines contain the migration edges. The first column for
# these lines is the weight on the edge, followed (optionally) by the jackknife estimate of the weight,
# the jackknife estimate of the standard error, and the p-values. Then come the subtree below the
# origin of the migration edge, and the subtree below the destination of the migration edge.
