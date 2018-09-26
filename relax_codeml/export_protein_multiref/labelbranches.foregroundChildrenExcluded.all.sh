
arrOutgroupSpp='APL,PLP';
arrUnusedClades='Aplocheilidae,root'; #c('Aplocheilidae', 'root'); #mark as unused clades "only for hyphy relax"
arrMarkUnusedCladeChildren='T,F'; #c(T , F);

#######################################################################
arrForegroundClades='Callopanchax'; #mark as foreground
arrMarkForegroundChildren='T';

sTaxonRequirements="min_taxon_requirements_Callopanchax.txt";

nMarkStyle='relax'; #either "codeml" or "relax"
sOutDIR="Relax_Callopanchax_childrenexcl";

mkdir -p $sOutDIR
Rscript labelbranches.foregroundChildrenExcluded.R $sOutDIR  $nMarkStyle  $sTaxonRequirements  $arrOutgroupSpp  $arrForegroundClades  $arrMarkForegroundChildren  $arrUnusedClades  $arrMarkUnusedCladeChildren > $sOutDIR/log.txt &


########################################################################

#######################################################################
arrForegroundClades='Nothobranchius'; #mark as foreground
arrMarkForegroundChildren='T';

sTaxonRequirements="min_taxon_requirements_Nothos.txt";

nMarkStyle='relax'; #either "codeml" or "relax"
sOutDIR="Relax_Nothobranchius_childrenexcl";

mkdir -p $sOutDIR
Rscript labelbranches.foregroundChildrenExcluded.R $sOutDIR  $nMarkStyle  $sTaxonRequirements  $arrOutgroupSpp  $arrForegroundClades  $arrMarkForegroundChildren  $arrUnusedClades  $arrMarkUnusedCladeChildren > $sOutDIR/log.txt &

########################################################################


#######################################################################
arrForegroundClades='Nothobranchius,Callopanchax'; #mark as foreground
arrMarkForegroundChildren='T,T';

sTaxonRequirements="min_taxon_requirements_CallNothos.txt";

nMarkStyle='relax'; #either "codeml" or "relax"
sOutDIR="Relax_CallNothos_childrenexcl";

mkdir -p $sOutDIR
Rscript labelbranches.foregroundChildrenExcluded.R $sOutDIR  $nMarkStyle  $sTaxonRequirements  $arrOutgroupSpp  $arrForegroundClades  $arrMarkForegroundChildren  $arrUnusedClades  $arrMarkUnusedCladeChildren > $sOutDIR/log.txt &

########################################################################


#######################################################################
arrForegroundClades='PronothoNothos'; #mark as foreground
arrMarkForegroundChildren='T';

sTaxonRequirements="min_taxon_requirements_ProNothos.txt";

nMarkStyle='relax'; #either "codeml" or "relax"
sOutDIR="Relax_PronothoNothos_childrenexcl";

mkdir -p $sOutDIR
Rscript labelbranches.foregroundChildrenExcluded.R $sOutDIR  $nMarkStyle  $sTaxonRequirements  $arrOutgroupSpp  $arrForegroundClades  $arrMarkForegroundChildren  $arrUnusedClades  $arrMarkUnusedCladeChildren > $sOutDIR/log.txt &

########################################################################


#######################################################################
arrForegroundClades='AKN'; #mark as foreground
arrMarkForegroundChildren='T';

sTaxonRequirements="min_taxon_requirements_AKN.txt";

nMarkStyle='relax'; #either "codeml" or "relax"
sOutDIR="Relax_AKN_childrenexcl";

mkdir -p $sOutDIR
Rscript labelbranches.foregroundChildrenExcluded.R $sOutDIR  $nMarkStyle  $sTaxonRequirements  $arrOutgroupSpp  $arrForegroundClades  $arrMarkForegroundChildren  $arrUnusedClades  $arrMarkUnusedCladeChildren > $sOutDIR/log.txt &

########################################################################


#######################################################################
arrForegroundClades='Aphyosemion'; #mark as foreground
arrMarkForegroundChildren='T';

sTaxonRequirements="min_taxon_requirements_Aphyosemion.txt";

nMarkStyle='relax'; #either "codeml" or "relax"
sOutDIR="Relax_Aphyosemion_childrenexcl";

mkdir -p $sOutDIR
Rscript labelbranches.foregroundChildrenExcluded.R $sOutDIR  $nMarkStyle  $sTaxonRequirements  $arrOutgroupSpp  $arrForegroundClades  $arrMarkForegroundChildren  $arrUnusedClades  $arrMarkUnusedCladeChildren > $sOutDIR/log.txt &

########################################################################


#######################################################################
arrForegroundClades='Epiplatys'; #mark as foreground
arrMarkForegroundChildren='T';

sTaxonRequirements="min_taxon_requirements_Epiplatys.txt";

nMarkStyle='relax'; #either "codeml" or "relax"
sOutDIR="Relax_Epiplatys_childrenexcl";

mkdir -p $sOutDIR
Rscript labelbranches.foregroundChildrenExcluded.R $sOutDIR  $nMarkStyle  $sTaxonRequirements  $arrOutgroupSpp  $arrForegroundClades  $arrMarkForegroundChildren  $arrUnusedClades  $arrMarkUnusedCladeChildren > $sOutDIR/log.txt &

########################################################################


#######################################################################
arrForegroundClades='Archiaphyosemion'; #mark as foreground
arrMarkForegroundChildren='T';

sTaxonRequirements="min_taxon_requirements_Archiaphyosemion.txt";

nMarkStyle='relax'; #either "codeml" or "relax"
sOutDIR="Relax_Archiaphyosemion_childrenexcl";

mkdir -p $sOutDIR
Rscript labelbranches.foregroundChildrenExcluded.R $sOutDIR  $nMarkStyle  $sTaxonRequirements  $arrOutgroupSpp  $arrForegroundClades  $arrMarkForegroundChildren  $arrUnusedClades  $arrMarkUnusedCladeChildren > $sOutDIR/log.txt &

########################################################################


#######################################################################
arrForegroundClades='Scriptaphyosemion'; #mark as foreground
arrMarkForegroundChildren='T';

sTaxonRequirements="min_taxon_requirements_Scriptaphyosemion.txt";

nMarkStyle='relax'; #either "codeml" or "relax"
sOutDIR="Relax_Scriptaphyosemion_childrenexcl";

mkdir -p $sOutDIR
Rscript labelbranches.foregroundChildrenExcluded.R $sOutDIR  $nMarkStyle  $sTaxonRequirements  $arrOutgroupSpp  $arrForegroundClades  $arrMarkForegroundChildren  $arrUnusedClades  $arrMarkUnusedCladeChildren > $sOutDIR/log.txt &
########################################################################


#######################################################################
arrForegroundClades='Fundulopanchax'; #mark as foreground
arrMarkForegroundChildren='T';

sTaxonRequirements="min_taxon_requirements_Fundulopanchax.txt";

nMarkStyle='relax'; #either "codeml" or "relax"
sOutDIR="Relax_Fundulopanchax_childrenexcl";

mkdir -p $sOutDIR
Rscript labelbranches.foregroundChildrenExcluded.R $sOutDIR  $nMarkStyle  $sTaxonRequirements  $arrOutgroupSpp  $arrForegroundClades  $arrMarkForegroundChildren  $arrUnusedClades  $arrMarkUnusedCladeChildren > $sOutDIR/log.txt &

########################################################################


wait


