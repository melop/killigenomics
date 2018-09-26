library(ape)
library(phytools)

#setwd("/beegfs/group_dv/home/RCui/killifish_genomes/hyphyrelax/export_protein_multiref_23_july_2017");

arrAUResults <- Sys.glob("genetrees_improved/ret_part*.txt"); #c("genetrees/ret_part0.of.200.txt");#Sys.glob("genetrees/ret_part*.txt")

############ Change this ##############################
sOutDIR <- "Relax_Nothos";
arrForegroundClades <- c('Nothobranchius'); #mark as foreground
arrMarkForegroundChildren <- c(T );

arrUnusedClades <- c(); #c('Aplocheilidae', 'root'); #mark as unused clades "only for hyphy relax"
arrMarkUnusedCladeChildren <- c();#c(T , F);

nMarkStyle <- 'relax';#either "codeml" or "relax"
sTaxonRequirements <- "min_taxon_requirements_Nothos.txt";
arrOutgroupSpp <- c('APL' , 'PLP');
############ Change this ##############################

args = commandArgs(trailingOnly=TRUE);

nArgCount <- 0;

nArgCount <- nArgCount + 1;
if (length(args) >= nArgCount) {
  sOutDIR <- args[nArgCount];
}

nArgCount <- nArgCount + 1;
if (length(args) >= nArgCount) {
  nMarkStyle <- args[nArgCount];
}


nArgCount <- nArgCount + 1;
if (length(args) >= nArgCount) {
  sTaxonRequirements <- args[nArgCount];
}

nArgCount <- nArgCount + 1;
if (length(args) >= nArgCount) {
  arrOutgroupSpp <- unlist(strsplit(args[nArgCount] , split = ","));
}

nArgCount <- nArgCount + 1;
if (length(args) >= nArgCount) {
  arrForegroundClades <- unlist(strsplit(args[nArgCount] , split = ","));
}


nArgCount <- nArgCount + 1;
if (length(args) >= nArgCount) {
  arrMarkForegroundChildren <- as.logical( unlist(strsplit(args[nArgCount] , split = ",")));
}


nArgCount <- nArgCount + 1;
if (length(args) >= nArgCount) {
  arrUnusedClades <- unlist(strsplit(args[nArgCount] , split = ","));
}


nArgCount <- nArgCount + 1;
if (length(args) >= nArgCount) {
  arrMarkUnusedCladeChildren <- as.logical( unlist(strsplit(args[nArgCount] , split = ",")));
}



dir.create(sOutDIR, recursive = T , mode = "0777")

sMonophylyDef <- "monophyly_desc.txt";

datTaxonReq <- read.table(sTaxonRequirements , header=T, stringsAsFactors = F);

fnBuildTree <- function(sNodeName, datMap, arrIncludeTaxon) {
  
  arrTree <- list();
  
  if ( !( sNodeName %in% datMap$V1) ) {
    return(paste("Node name", sNodeName," undefined in map\n"));
  }
  
  arrChilds <- unlist( strsplit( datMap[datMap$V1==sNodeName , 'V2'], ','));
  
  
  for( sChild in arrChilds ) {
    sChild <- trimws(sChild);
    if ( substr(sChild , 1,1) == '*' ) { 
      sChildNodeName <- substr(sChild, 2, nchar(sChild) ); #delete * from taxon name
      #echo(sChildNodeName."\n");
      oChildNode <- fnBuildTree(sChildNodeName , datMap , arrIncludeTaxon);
      if ( length(oChildNode) > 0 ) {
        arrTree[[sChildNodeName]] <- oChildNode;
      } 
    } else {
      if ( isTRUE(arrIncludeTaxon) || (sChild %in% arrIncludeTaxon) ) {
        arrTree[sChild] <- T; #terminal taxon
      }
    }
  }
  
  return(arrTree);
}


fnParseMonoDesc <- function(sIn, arrIncludeTaxon = T) {
  datMap <- read.table(sMonophylyDef , header=F, stringsAsFactors = F);
  
  arrTree <-fnBuildTree('root', datMap , arrIncludeTaxon);
    
  return(arrTree);
}

fnGetCountForTaxon <- function(sTaxon, lsGroup) {
  nCounts <- 0;
  
  for(sName in names(lsGroup) ) {
    if (sName == sTaxon) {
      #cat(sName , " ", sTaxon, "\n");
      if ( isTRUE(lsGroup[[sName]]) ) {
        return(1);
      } else {
        #get all children
        for(sChildName in names(lsGroup[[sName]]) ) {
          nCounts <- nCounts + fnGetCountForTaxon(sChildName , lsGroup[[sName]]);
        }
        return(nCounts);
      }
    } else { #search in lower level
      #get all children
        nCounts <- nCounts + fnGetCountForTaxon(sTaxon , lsGroup[[sName]]);
    }
  }
  
  return(nCounts);
}

fnGetAllDescendantTips <- function(sTaxon, lsGroup)  {
  arrTips <- c();
  
  for(sName in names(lsGroup) ) {
    if (sName == sTaxon) {
      #cat(sName , " ", sTaxon, "\n");
      if ( isTRUE(lsGroup[[sName]]) ) {
        return(sName);
      } else {
        #get all children
        for(sChildName in names(lsGroup[[sName]]) ) {
          arrTips <- c(arrTips,  fnGetAllDescendantTips(sChildName , lsGroup[[sName]] ));
        }
        return(arrTips);
      }
    } else { #search in lower level
      #get all children
      arrTips <- c( arrTips , fnGetAllDescendantTips(sTaxon , lsGroup[[sName]]) );
    }
  }
  
  return(arrTips);
}


datAU <- NULL;
for( sAUFile in arrAUResults) {
  datAU <- rbind(datAU, read.table(sAUFile,header=F, fill = T, strip.white = F, stringsAsFactors = F, na.strings = "") );
}

colnames(datAU) <- c('GeneId' , 'nTaxa' , 'status' , 'obs', 'au' , 'tree_free', 'tree_constraint');

datAUFilter <- datAU[datAU$status=='accepted_monophyly', ];
#View(datAUFilter);

for(nR in 1:nrow(datAUFilter)) {
#for(nR in 1) {
    
  sGeneId <- datAUFilter[nR, 'GeneId'];
  sTree <- datAUFilter[nR, 'tree_constraint'];
  oTree <- read.tree(text=sTree);
  file.remove( paste(sOutDIR ,'/', sGeneId , '.txt' , sep=""), showWarnings=F );

  if (is.null(oTree)) {
	cat("Error reading tree for ", sGeneId, "\n", "Tree:", sTree , "\n");
        next;
  }
  oTree <- unroot( oTree );
  arrOutgroupSppFound <- intersect(arrOutgroupSpp, oTree$tip.label);
  cat("Doing " , sGeneId,"\n");
  if (length(arrOutgroupSppFound) ==0) {
    cat("ERROR: All outgroup species are missing, cannot root");
    next;
  }
  
  if (length(arrOutgroupSppFound ) == 1) {
    cat( sGeneId," only one outgroup \n");
    
    oTree <- reroot(tree = oTree, node.number = which(oTree$tip.label == arrOutgroupSppFound[1]));
  } else {
    oTree <- reroot(tree = oTree, node.number = findMRCA(oTree, tips=arrOutgroupSppFound,type = 'node' ));
  }
  oTree$edge.length <- NULL;

  #now check the number of taxa for each clade for this analysis:
  oClades <- list();
  oClades[['root']] <- fnParseMonoDesc(sMonophylyDef , oTree$tip.label);
  
  bValidated <- T;
  for(nR2 in 1:nrow(datTaxonReq) ){
    sTaxonToCheck <- datTaxonReq[nR2, 1];
    nMinTaxon <- datTaxonReq[nR2, 2];
    nCurrentTaxaCount <- fnGetCountForTaxon(sTaxonToCheck , oClades) ;
    if ( nCurrentTaxaCount < nMinTaxon) {
      cat("Not enough taxa for gene " , sGeneId , " in clade " , sTaxonToCheck ," min : ", nMinTaxon, " only ", nCurrentTaxaCount ,  "\n");
      bValidated <- F;
    }
  }
  
  if (!bValidated) next;
  if (nMarkStyle == 'relax') {
    #this is relax style marking. first need to mark every internal node as reference {R}
    nInternalNodeStart <- length(oTree$tip.label) ;
    nInternalNodeEnd <- max(oTree$edge);
    oTree$node.label <- rep('{R}' , nInternalNodeEnd - nInternalNodeStart);
    
    bErr <- F;
    
    
    arrForSeq <- c();
    if (length(arrForegroundClades)>0) {
      arrForSeq <- 1:length(arrForegroundClades);
    }
    
  
    for(nForeground in arrForSeq)  {
      sForegroundClade <- arrForegroundClades[nForeground];
      bMarkChildren <- arrMarkForegroundChildren[nForeground];
      arrDescTips <- fnGetAllDescendantTips(sForegroundClade , oClades);
      if (is.null(arrDescTips) ) {
        bErr <-T;
        break;
      }
      if (length(arrDescTips) ==1) { #if this internal node only contains a single descendant
        sTip1 <- arrDescTips[1];
        oTree$tip.label[oTree$tip.label == sTip1] <- paste(sTip1, '{T}', sep=''); 
        next;
      }

      if (bMarkChildren) {
        for(sTip1 in arrDescTips) {
          oTree$tip.label[oTree$tip.label == sTip1] <- paste(sTip1, '{U}', sep=''); 
          for(sTip2 in arrDescTips) {
            if (sTip1 == sTip2) next;
            #cat(sTip1, " " , sTip2, "\n");
            oTree <- makeNodeLabel(oTree, method='u', nodeList = list('{U}' =  c(sTip1, sTip2) ) );
          }
        }
      }
      oTree <- makeNodeLabel(oTree, method='u', nodeList = list('{T}' =  arrDescTips ) );
    }
    
    if (bErr) {
      cat("Some foreground clades do not exist in the current phylogeny, may be due to missing data. Skip \n"  );
      next;
    }
    
    arrForSeq <- c();
    if (length(arrUnusedClades)>0) {
      arrForSeq <- 1:length(arrUnusedClades);
    }
    
    bErr <- F;
    for(nUnused in arrForSeq ) {
      sUnusedClade <- arrUnusedClades[nUnused];
      bMarkChildren <- arrMarkUnusedCladeChildren[nUnused];
      arrDescTips <- fnGetAllDescendantTips(sUnusedClade , oClades);
      if (is.null(arrDescTips) ) {
        bErr <-T;
        break;
      }
      if (length(arrDescTips) ==1) { #if this internal node only contains a single descendant
        sTip1 <- arrDescTips[1];
        oTree$tip.label[oTree$tip.label == sTip1] <- paste(sTip1, '{U}', sep=''); 
        next;
      }
      oTree <- makeNodeLabel(oTree, method='u', nodeList = list('{U}' =  arrDescTips ) );
      if (bMarkChildren) {
        for(sTip1 in arrDescTips) {
          oTree$tip.label[oTree$tip.label == sTip1] <- paste(sTip1, '{U}', sep=''); 
          for(sTip2 in arrDescTips) {
            if (sTip1 == sTip2) next;
            cat(sTip1, " " , sTip2, "\n");
            oTree <- makeNodeLabel(oTree, method='u', nodeList = list('{U}' =  c(sTip1, sTip2) ) );
          }
        }
      }
    }
    
    if (bErr) {
      cat("Some unconsidered {U} clades do not exist in the current phylogeny, may be due to missing data. Skip \n"  );
      next;
    }
    
    #at the end, mark any unmarked tips as {R}
    
    for(sTip in oTree$tip.label) {
        if ( length( grep('{' , sTip, perl=T)) == 1) {
          next; #already marked
        }
      oTree$tip.label[oTree$tip.label == sTip] <- paste(sTip, '{R}', sep=''); 
      
    }
    
    
    
  }
  if (nMarkStyle == 'codeml') {
    #this is codeml style marking. 

    bErr <- F;
    for(nForeground in 1:length(arrForegroundClades) ) {
      sForegroundClade <- arrForegroundClades[nForeground];
      bMarkChildren <- arrMarkForegroundChildren[nForeground];
      arrDescTips <- fnGetAllDescendantTips(sForegroundClade , oClades);
      if (is.null(arrDescTips) || length(arrDescTips) ==1 ) {
        if (is.null(arrDescTips)) {
          bErr <-T;
        }
        break;
      }
      sCladeMark <- "#1";
      if (bMarkChildren) {
        sCladeMark <- "$1";
      }
      
      if (length(arrDescTips) == 1) {
        sTip1 <- arrDescTips[0];
        oTree$tip.label[oTree$tip.label == sTip1] <- paste(sTip1, '#1', sep=''); 
      } else {
        lsCm <- list();
        lsCm[[sCladeMark]] <- arrDescTips;
        oTree <- makeNodeLabel(oTree, method='u', nodeList =  lsCm  );
      }
 
    }
    
    if (bErr) {
      cat("Some foreground clades do not exist in the current phylogeny, may be due to missing data. Skip \n"  );
      next;
    }
    
   
  }
  
  #plot(oTree , show.node.label  = T);
  write.tree(oTree , file=paste(sOutDIR ,'/', sGeneId , '.txt' , sep="") );
}




