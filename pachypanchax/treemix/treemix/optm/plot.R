setwd("/beegfs/group_dv/home/RCui/killifish_genomes/plp_treemix/treemix/optm");
library(OptM);

folder <- '.';


test.optM <- optM(folder, tsv="evano.tsv")
plot_optM(test.optM, method = "Evanno", pdf = "evano.pdf")


test.Sizer <- optM(folder, tsv="sizer.tsv", method = "SiZer")
plot_optM(test.Sizer, pdf = "sizer.pdf", method = "SiZer")

test.Linear <- optM(folder, tsv="linear.tsv", method = "linear")
plot_optM(test.Linear, pdf = "linear.pdf", method = "linear")


