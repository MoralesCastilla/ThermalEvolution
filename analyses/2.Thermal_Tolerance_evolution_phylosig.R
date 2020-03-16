#'##############################################################################
#' Script to run:
#' * Phylogenetic comparative analyses of thermal tolerance data
#'   
#'  by Ignacio Morales-Castilla
#'  started June 2016
#'##############################################################################


#### to start ####
rm(list=ls())
options(stringsAsFactors=FALSE)



#### source data ####
source("1.Thermal_Tolerance_evolution_data.R")


#### load extra packages ####
library(phangorn)



#'###############################
#### some more data cleaning ####


## define parameters to subset data
parslist <- list(c("ecto","No.plants"),   # all ectothermic organisms
               c("endo","No.plants"),   # all endotherms (mammals & birds)
               c("ecto","Plants.only"), # only plants (excluding algae)
               c("ecto","Plants"))      # all plants (including algae)


## adding species lacking from the phylogeny using congeneric.merge (package PEZ)
names.all.data <- paste(sub("^\\s+|\\s+$", "", GTherm.data$Genus, fixed=F),
                      sub("^\\s+|\\s+$", "", GTherm.data$Species, fixed=F), sep="_")
names.all.data <- sub("[.]", "", names.all.data)
names.clean <- unlist(lapply(strsplit(names.all.data, "_"),
                             function(x){paste(x[1], x[2], sep="_")}))


names.to.add <- names.clean[which(!names.clean %in% ttol.GTherm.phy$tip.label)]
ttol.GTherm.phy <- congeneric.merge(ttol.GTherm.phy, names.to.add, split = "_")



#'############################################
#### calculate species level phylo-signal ####

## generate array to store results
storing.phylosig <- array(NA, dim = c(10, 6, nrep))
colnames(storing.phylosig) <- c("n","sigma2","OU.alpha","lnl.OU","lnl.BM","lnl.WN")
rownames(storing.phylosig) <- c("all.up","all.low",
                              "ecto.up","ecto.low",
                              "endo.up","endo.low",
                              "plants.only.up","plants.only.low",
                              "plants.up","plants.low")

## check for which metrics there are data
table(GTherm.data$max_metric)
table(GTherm.data$min_metric)


## loop to compute signal across subsets of data
impar = seq(1, 100, 2)
par = seq(2, 100, 2)

for(i in 1:5){
  print(i)
  
  ## subset species by groups 
  if(i == 1){GTherm.final <- GTherm.data}
  
  if(i %in% 2:3){
    GTherm.final <- subset(GTherm.data,GTherm.data$thermy == parslist[[i-1]][1]
                        & GTherm.data$Plant == parslist[[i-1]][2])
  } 
  if (i == 4){
    GTherm.final <- subset(GTherm.data,GTherm.data$thermy == parslist[[i-1]][1]
                        & GTherm.data$Plant.only == parslist[[i-1]][2])
  }
  if (i == 5){
    GTherm.final <- subset(GTherm.data,GTherm.data$thermy == parslist[[i-1]][1]
                        & GTherm.data$Plant == parslist[[i-1]][2])
  }
  
  
  GTherm.final$species.phylo <- paste(sub("^\\s+|\\s+$", "",
                                          GTherm.final$Genus,fixed=F),
                                   sub("^\\s+|\\s+$", "",
                                       GTherm.final$Species,fixed=F), sep="_")
  
  GTherm.final$species.phylo <- sub("[.]", "", GTherm.final$species.phylo)
  names.clean <- unlist(lapply(strsplit(GTherm.final$species.phylo, "_"),
                               function(x){paste(x[1], x[2], sep = "_")}))
  GTherm.final$species.phylo <- names.clean
  
  ## subset by tolerance metrics (for each group)
  if(i == 1){
    x.up <- GTherm.final$Tmax[which(GTherm.final$max_metric == "ctmax" | 
                                      GTherm.final$max_metric == "UTNZ" | 
                                      GTherm.final$max_metric == "LT0")]
    x.low <- GTherm.final$tmin[which(GTherm.final$min_metric == "ctmin" | 
                                       GTherm.final$min_metric == "LTNZ"| 
                                       GTherm.final$min_metric == "LT50")]
    
    names(x.up) <- GTherm.final$species.phylo[which(GTherm.final$max_metric == "ctmax" | 
                                                      GTherm.final$max_metric == "UTNZ"| 
                                                      GTherm.final$max_metric == "LT0")]
    names(x.low) <- GTherm.final$species.phylo[which(GTherm.final$min_metric == "ctmin" | 
                                                       GTherm.final$min_metric == "LTNZ"| 
                                                       GTherm.final$min_metric == "LT50")]
    
  }
  
  if(i %in% 2:3){
    
    x.up <- GTherm.final$Tmax[which(GTherm.final$max_metric == "ctmax" | 
                                      GTherm.final$max_metric == "UTNZ")]
    x.low <- GTherm.final$tmin[which(GTherm.final$min_metric == "ctmin" | 
                                       GTherm.final$min_metric == "LTNZ")]
    
    names(x.up) <- GTherm.final$species.phylo[which(GTherm.final$max_metric == "ctmax" | 
                                                      GTherm.final$max_metric == "UTNZ")]
    names(x.low) <- GTherm.final$species.phylo[which(GTherm.final$min_metric == "ctmin" | 
                                                       GTherm.final$min_metric == "LTNZ")]
    
  } 
  
  if (i == 4){
    x.up <- GTherm.final$Tmax[which(GTherm.final$max_metric == "ctmax")]
    x.low <- GTherm.final$tmin[which(GTherm.final$min_metric == "LT50")]
    
    names(x.up) <- GTherm.final$species.phylo[which(GTherm.final$max_metric == "ctmax")]
    names(x.low) <- GTherm.final$species.phylo[which(GTherm.final$min_metric == "LT50")]
  
    }
  
  if (i == 5){
    x.up <- GTherm.final$Tmax[which(GTherm.final$max_metric == "LT0")]
    x.low <- GTherm.final$tmin[which(GTherm.final$min_metric == "LT50")]
    
    names(x.up) <- GTherm.final$species.phylo[which(GTherm.final$max_metric == "LT0")]
    names(x.low) <- GTherm.final$species.phylo[which(GTherm.final$min_metric == "LT50")]
    
  }
  
  ## manipulate phylogeny for data subset (prune, randomize polytomies)
  phy.t <- drop.tip(ttol.GTherm.phy,which(!ttol.GTherm.phy$tip.label %in% 
                                             GTherm.final$species.phylo))
  if(Ntip(phy.t)-Nnode(phy.t)>1){
    foo  <-  function(node, tree) length(Children(tree, node))
    nodes  <-  1:phy.t$Nnode + Ntip(phy.t)
    nodes  <-  nodes[sapply(1:phy.t$Nnode + Ntip(phy.t), 
                            foo, tree = phy.t) >2]
    nodes.2.solve  <-  unlist(lapply(lapply(nodes,
                                            function(x){Children(phy.t, x)}),length))
    nnodes.solve <- sum(nodes.2.solve<6)
    list.objs <- list()
    counter <- 1
    for(ob in 1:nnodes.solve){
      #counter <- counter+1
      treessolved <- resolveNode(phy.t,nodes[which(nodes.2.solve<6)[ob]])
      if(class(treessolved)=='phylo'){
        
        list.objs[[counter]] <- treessolved
        counter <- counter+1
      } else {
        ntres <- length(treessolved)
        for(nt in 1:ntres){
          #print(paste(ob,nt,counter))
          list.objs[[counter]] <- treessolved[[nt]]
          counter <- counter+1
        }
      }
      
    }
    class(list.objs) <- "multiPhylo"
  } else{
    list.objs <- list(phy.t)
    class(list.objs) <- "multiPhylo"
  }
  
  # set interval (how many trees to use to compute phylogenetic signal)
  interval=1:length(list.objs) # more options below if 
  #if(length(list.objs) == 1){interval = 1}
  #if(length(list.objs)>100){interval=1:100}
  #if(length(list.objs)>10){interval = 1:10}
  #if(length(list.objs)>1 & length(list.objs)<=10){ interval=1:length(list.objs)}
  
  for(j in interval){
    
    print(paste(i,j))
    phy.t <- multi2di(list.objs[[j]])
    phy.t <- nnls.tree(cophenetic(phy.t),phy.t,rooted=TRUE)
    
    
    ## create comparative data
    comp.sps.tmax <- comparative.data(phy.t, 
                                      data.frame(x.up, namestmax = names(x.up)),
                                      names.col="namestmax")
    comp.sps.tmin <- comparative.data(phy.t,
                                      data.frame(x.low,namestlow = names(x.low)),
                                      names.col="namestlow")
    
    
    ## declare phylogenies
    phy.low <- comp.sps.tmin$phy
    phy.up <- comp.sps.tmax$phy
    
    ## declare and name trait vectors
    vec.up <- as.numeric(as.vector(comp.sps.tmax$data[,1]))
    vec.low <- as.numeric(as.vector(comp.sps.tmin$data[,1]))
    names(vec.up) <- rownames(comp.sps.tmax$data)
    names(vec.low) <- rownames(comp.sps.tmin$data)
    
    ## compute phylogenetic signal
    phy.upOU <- rescaleTree(comp.sps.tmax$phy,1)
    phy.lowOU <- rescaleTree(comp.sps.tmin$phy,1)
    up.OU <- fitContinuous(phy.upOU,vec.up,model = "OU", bounds = list(alpha = c(0.0005,20)))
    low.OU <- fitContinuous(phy.lowOU,vec.low,model = "OU",bounds = list(alpha = c(0.0005,20)))
    up.BM <- fitContinuous(comp.sps.tmax$phy,vec.up,model = "BM")
    low.BM <- fitContinuous(comp.sps.tmin$phy,vec.low,model = "BM")
    up.lambda <- fitContinuous(comp.sps.tmax$phy,vec.up,model = "white")
    low.lambda <- fitContinuous(comp.sps.tmin$phy,vec.low,model = "white")
    
    ## store phylogenetic signal results
    storing.phylosig[impar[i],1,j] = length(comp.sps.tmax$phy$tip.label)
    storing.phylosig[impar[i],2,j] = up.BM$opt$sigsq
    storing.phylosig[impar[i],3,j] = log(up.OU$opt$alpha,10)*(-1)
    storing.phylosig[impar[i],4,j] = up.OU$opt$lnL
    storing.phylosig[impar[i],5,j] = up.BM$opt$lnL
    storing.phylosig[impar[i],6,j] = up.lambda$opt$lnL
    storing.phylosig[par[i],1,j] = length(comp.sps.tmin$phy$tip.label)
    storing.phylosig[par[i],2,j] = low.BM$opt$sigsq
    storing.phylosig[par[i],3,j] = log(low.OU$opt$alpha,10)*(-1)
    storing.phylosig[par[i],4,j] = low.OU$opt$lnL
    storing.phylosig[par[i],5,j] = low.BM$opt$lnL
    storing.phylosig[par[i],6,j] = low.lambda$opt$lnL
    
  }
}

## inspect results
storing.phylosig


#### end ####







