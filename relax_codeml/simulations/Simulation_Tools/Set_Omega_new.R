### Script reads Tree, Forground Clade,Omega value for Forground Cladeand
#Omega value to all others as Background Clade from standart input###
## It expects only one Tree !!!!
#reroots tree
#sets Forground Clade
#sets Omega value for Forground Clade
#sets Omega value to all others as Background Clade

library(ape)
library(phytools)

# read a Tree from standart input

# transformes the phylogenie input
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

# calls fnBuildTree to set the root of the tree
fnParseMonoDesc <- function(sIn, arrIncludeTaxon = T) {
  datMap <- read.table(sMonophylyDef , header=F, stringsAsFactors = F);
  
  arrTree <-fnBuildTree('root', datMap , arrIncludeTaxon);
  
  return(arrTree);
}


#print("Before loop")
fStdin <- file("stdin")
open(fStdin)

### change this according to directory on amalia ###
sMonophylyDef <- "/beegfs/group_dv/home/LIasi/Simulation_Tools/monophyly_desc.txt";

#sMonophylyDef <- "/Users/liasi/Desktop/PHP_master_script/Set_Omega/monophyly_desc.txt";

fnProcessFromPHPInput <- function (Input1) {
  #Input1 <- readLines(x,n=1, warn = FALSE)
  #print("Starting loop")
  if (is.null(Input1)) {
    cat(" STDIN Input is not there \n")
  }
  #print(as.character(unlist(Input1)))
  Input <- read.table(text = Input1, header=F, fill = T, strip.white = F, stringsAsFactors = F, na.strings = "", sep = "\t")
  if (is.null(Input)) {
    cat("ERROE: No Input")
  } 

  TreeIn <-  Input[,1]
  if (is.null(TreeIn)) {
    cat("ERROE: Input not assigned to variable = TreeIn \n")
  } 
  #cat(as.character(unlist(TreeIn)))
  
  # variable defining the name of the forground clade
  IsForgroundClade <- as.character(tail(Input$V2,n=1))
  if (length(IsForgroundClade)==0) {
    print("ERROE: No Foreground Clade specified")
    
  }
  # variable set to true for including al children of the forground clade
  arrMarkForegroundChildren <- c(T );
  
  # variable for defining the Omega value for the Forground clade
  OmegaForground <- as.numeric(tail(Input$V3,n=1))
  if (length(OmegaForground)==0) {
    print("ERROE: Foreground Omega not given")
    
  }
  # variable for defining the Omega value for the Background clade
  OmegaBackground <- as.numeric(tail(Input$V4,n=1))
  if (length(OmegaBackground)==0) {
    print("ERROR: Background Omega not given")
    
  }
  
  # variable setting the Outgroups
  areOutgroups <- c('APL' , 'PLP');
  
  
  
  ## read in the tree given by standart input
  #open(TreeIn)
  datAU <- NULL;
  #datAU <- rbind(datAU, read.table(file = TreeIn,header=F, fill = T, strip.white = F, stringsAsFactors = F, na.strings = "") );
  datAU <- rbind(datAU,TreeIn)
  colnames(datAU) <- c('tree_constraint');
  

  sTree <- datAU[, 'tree_constraint'];
  if (is.null(sTree)) {
    cat("sTree unsuccesfull \n");
  }
  oTree <- read.tree(text=sTree);
  if (is.null(oTree)) {
    cat("oTree unsuccesfull \n");
  } 
  # reroots tree woth Outgroups set in variable areOutgroups
  oTree <- unroot( oTree );
  arrOutgroupSppFound <- intersect(areOutgroups, oTree$tip.label);
  
  if (length(arrOutgroupSppFound ) == 1) {
    #cat("* Warning: Tree has only one outgroup \n");
    oTree <- reroot(tree = oTree, node.number = which(oTree$tip.label == arrOutgroupSppFound[1]));
  } else {
    oTree <- reroot(tree = oTree, node.number = findMRCA(oTree, tips=arrOutgroupSppFound,type = 'node' ));
  }
  
  oTree$edge.length <- NULL;
  
  # setting the node lable
  nInternalNodeStart <- length(oTree$tip.label) ;
  nInternalNodeEnd <- max(oTree$edge);
  oTree$node.label <- rep('#BcGr' , nInternalNodeEnd - nInternalNodeStart);
  
  # reads in the phylogenie of the Aplocheiidae
  #sMonophylyDef <- "monophyly_desc.txt";
  if (length(sMonophylyDef)==0) {
    print("ERROE: n Monopholy def incorrect \n")
    
  }
  

  

  
  # sets the Tips for the Tree

  #cat("middle of loop")
  oClades <- list();
  oClades[['root']] <- fnParseMonoDesc(sMonophylyDef , oTree$tip.label);
  
  arrForSeq <- c();
  if (length(IsForgroundClade)>0) {
    arrForSeq <- 1:length(IsForgroundClade);
  }
  
  # sets the forground clades and marks them
  for(nForeground in arrForSeq)  {
    sForegroundClade <- IsForgroundClade[nForeground];
    bMarkChildren <- arrMarkForegroundChildren[nForeground];
    arrDescTips <- fnGetAllDescendantTips(sForegroundClade , oClades);
    if (is.null(arrDescTips) ) {
      bErr <-T;
      break;
    }
    if (length(arrDescTips) ==1) { #if this internal node only contains a single descendant
      sTip1 <- arrDescTips[1];
      oTree$tip.label[oTree$tip.label == sTip1] <- paste(sTip1, '#ForGr', sep=''); 
      next;
    }
    oTree <- makeNodeLabel(oTree, method='u', nodeList = list('#ForGr' =  arrDescTips ) );
    if (bMarkChildren) {
      for(sTip1 in arrDescTips) {
        oTree$tip.label[oTree$tip.label == sTip1] <- paste(sTip1, '#ForGr', sep=''); 
        for(sTip2 in arrDescTips) {
          if (sTip1 == sTip2) next;
          #cat(sTip1, " " , sTip2, "\n");
          oTree <- makeNodeLabel(oTree, method='u', nodeList = list('#ForGr' =  c(sTip1, sTip2) ) );
        }
      }
    }
  }
  # marks every unmarked tip as background
  for(sTip in oTree$tip.label) {
    if ( length( grep('#' , sTip, perl=T)) == 1) {
      next; #already marked
    }
    oTree$tip.label[oTree$tip.label == sTip] <- paste(sTip, '#BcGr', sep=''); 
    
  }
  
  oTree$tip.label=sub('BcGr', paste(OmegaBackground), oTree$tip.label)
  oTree$tip.label=sub('ForGr', paste(OmegaForground), oTree$tip.label)
  oTree$node.label=sub('BcGr', paste(OmegaBackground), oTree$node.label)
  oTree$node.label=sub('ForGr', paste(OmegaForground), oTree$node.label)
  
  #print("Still in loop")
  
    if (is.null(oTree)){
      cat("no tree written \n");
    }

  #write(file="",x = paste( (write.tree(oTree )), "\n", sep=""))

  #write(file = "",x=write.tree(oTree ))
  #print("Allmost finished with loop")
  #cat(paste("HERE IS THE TREE: ",write.tree(oTree )," END OF THE TREE"))
  cat(write.tree( oTree ),"\n");
  #cat("aaaa\n", sep = "");

}

while(length(sLn <- readLines(fStdin,n=1, warn = T)) > 0) {
  #cat("doing line ", sLn );
  fnProcessFromPHPInput(sLn);
  }
 
#fnReadFromPHP(readsource=fStdin)
#print("Out of Loop")
