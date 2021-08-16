# ignore comment below we are running this script 1000 times ones 

# housekeeping
rm(list=ls())
graphics.off()

# packages
require(caper)
require(phytools)
require(phylobase)
require(data.table)
require(geiger)
require(pez)

#### 

setwd("../../../oenon/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/oslin/EDGE2/")
load("NaResults.Rdata")

iters <- na.results$Am

### 

# which pext value to pick

load("Am_newpext.Rdata")

# and make sure it's in the right format

pext.1000$X1 <- as.character(pext.1000$X1)

# decide which tree to load

# and put it in the correct formet

# also load in species data
load("AmSpecies.Rdata")


# select the pext vals to use

pext <- data.frame(Species = as.character(pext.1000$X1), pext = NA, stringsAsFactors = F)

# Now we calculate EDGE 2.0 for each phylogeny and combine the results

# EDGE 2.0 functionhange 
# tree is a single phylo object and pext must be a dataframe with columns: 
# Species - the names of all species in the tree
# pext - the pext of each species
# iteration - the numeric value to be assigned to the run - e.g. the i in a for loop of 'i in 1:n'

EDGE.2.calc <- function(tree, pext){
  require(phylobase)
  require(data.table)
  require(caper)
  # create df for data
  if(!class(tree) == "phylo"){
    tree <- as(tree, "phylo")
  }
  tree_dat <- data.frame(Species = as.character(unique(tree$tip.label)),
                         TBL = tree$edge.length[sapply(c(1:length(tree$tip.label)),
                                                       function(x,y) which(y==x),y=tree$edge[,2])], 
                         pext = NA, ED = NA, EDGE = NA)
  
  # converts tree to phylo4 object
  tree <- as(tree, "phylo4")
  names(pext) <- c("species","pext")
  for(i in 1:length(tree_dat$Species)){
    tree_dat$pext[i] <- pext$pext[pext$species == tree_dat$Species[i]]
  }
  # calculate IWE and TBL for each species
  nodes <- descendants(tree, rootNode(tree), "all")
  for(i in 1:length(nodes)){
    tips <- descendants(tree, nodes[i], "tips")
    tips <- names(tips)
    tipscores <- which(pext$species %in% tips)
    tree@edge.length[which(tree@edge[,2] == nodes[i])] <- edgeLength(tree, nodes[i])*prod(pext$pext[tipscores])
    #print(paste("Node",i,"of",length(nodes),"transformed!", sep = " "))
  }
  #plot(tree)
  print("Tree transformed")
  for(i in 1:length(tree_dat$Species)){
    tree_dat$EDGE[i] <- sum(tree@edge.length[which(tree@edge[,2] %in% ancestors(tree, 
                                                                                which(tipLabels(tree) == tree_dat$Species[i]), "ALL"))], na.rm=T)
    tree_dat$ED[i] <- tree_dat$EDGE[i] / tree_dat$pext[i] 
    #print(paste("EDGE 2.0 calculated for species",i,"of",length(tipLabels(tree)),"!",sep=" ")) 
  }
  print("Edge scores calculated")
  tree <- as(tree, "phylo")
  edge.res <- list(tree_dat,tree)
  return(edge.res)
}


# now calculate EDGE2

iters <- iters[-c(1:8)] # want to start at 58, as that hasn't been completed yet
iters <- iters[-c(1:21)]
iters <- iters[-1] # start with 82

  # calculate EDGE 2

for (i in iters){
    print(i)
     # create lists
  # List to store EDGE scores across distribution of trees
    EDGE.2.list <- list()
    exp.PD.list <- list()
    exp.PD.trees <- list()
    # load tree
    tree.to.load <- paste("D:/Documents/ResearchProject/EDGE2/Amphibians/AmTrees/AmTree_",
          i, ".Rdata", sep = "")
    load(tree.to.load)
    
    # select pext column
    pext$pext <- pext.1000[,(i+1)]
    
    # run EDGE code
    res  <- EDGE.2.calc(tree, pext)
    res2 <- data.frame(res[[1]],above.median = 0)
    res2$above.median[res2$EDGE > median(res2$EDGE)] <- 1
    EDGE.2.list[[length(EDGE.2.list)+1]] <- res2[order(res2$EDGE,decreasing = T),]
    exp.PD.list[[length(exp.PD.list)+1]] <- data.frame(PD = sum(tree$edge.length),ePDloss = sum(res[[2]]$edge.length),
                                                       ePD = (sum(tree$edge.length) - sum(res[[2]]$edge.length)))
    exp.PD.trees[[length(exp.PD.trees)+1]] <- res[[2]]
   
    # save down results
    name_to_save <- paste("AmResults/Am_New_EDGE2",i, ".Rdata", sep = "_")
  
    save(exp.PD.trees, exp.PD.list, EDGE.2.list, file = name_to_save)

}
