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
library(party)



#'############################### 
#### Random Forest modelling ####

## The models analyze 3 hypotheses linked to 3 predictor variables:
## Age at the order level; Origination under cold-warm palaeocliamte; 
## Temperature where sampled


## Generate multi-dimensional array to store results
store.random.forest<-array(NA,dim=c(10,8,100))


## Generate lists of parameters to subset data
parslist <- list(c("Annelida", "Arthropoda", "Brachiopoda", "Bryozoa", "Chordata", "Echinodermata", "Mollusca"), 
               c("Chordata"), 
               c("Streptophyta"), 
               c("Chlorophyta", "Phaeophyceae", "Rhodophyta"), 
               c("Streptophyta", "Chlorophyta", "Phaeophyceae", "Rhodophyta"), 
               c("Ascomycota", "Basidiomycota"))
parslist2 <- c("ecto", "endo", "ecto", "ecto", "ecto")
metrics.max <- c("ctmax", "UTNZ", "ctmax", "LT0", "LT0")
metrics.min <- c("ctmin", "LTNZ", "LT50", "LT50", "LT50")


## Start loop to run models data
par = seq(2, 100, 2)
impar = seq(1, 100, 2)
set.seed(0)
for(i in 1:5){
  print(i)
  subs.data <- subset(GTherm.data, GTherm.data$thermy == parslist2[i] & 
                        GTherm.data$Phylum %in% parslist[[i]])
  
  ## make sure variables are numeric
  subs.data$order.age <- as.numeric(subs.data$order.age)
  subs.data$temp_max <- as.numeric(subs.data$temp.max)
  subs.data$temp_min <- as.numeric(subs.data$temp.min)
  subs.data$tmin <- as.numeric(subs.data$tmin)
  subs.data$Tmax <- as.numeric(subs.data$Tmax)
  #subs.data$max_metric
  
  
  subs.data.max = subset(subs.data, !is.na(Tmax) & max_metric  ==  metrics.max[i])
  subs.data.min = subset(subs.data, !is.na(tmin) & min_metric  ==  metrics.min[i])
  
  
  ## iterate models
  for(j in 1:100){
    set.seed(j)
    print(paste(i, j))
    mod.tmax.cf  <-  cforest(Tmax ~ order.age + warm.cold + temp_max,  data  =  subs.data.max,  
                           control  =  cforest_unbiased(mtry  =  2, ntree = 500))
    mod.tmin.cf  <-  cforest(tmin ~ order.age + warm.cold + temp_min,  data  =  subs.data.min,  
                           control  =  cforest_unbiased(mtry  =  2, ntree = 500))
    mod.tmax <- randomForest(Tmax ~ order.age+ warm.cold + temp_max, data = subs.data, 
                           ntree = 500, na.action = na.omit, importance = T, 
                           nPerm = 50, mtry = 2, corr.bias = T, localImp = T)
    mod.tmin <- randomForest(tmin ~ order.age + warm.cold + temp_min, data = subs.data, 
                           ntree = 500, na.action = na.omit, importance = T, 
                           nPerm = 50, mtry = 2, corr.bias = T, localImp = T)
    
    
    store.random.forest[impar[i], 1, j] <- mod.tmax.cf@responses@nobs
    store.random.forest[par[i], 1, j] <- mod.tmin.cf@responses@nobs

    store.random.forest[impar[i], 2, j] <- mod.tmax$rsq[500]
    store.random.forest[par[i], 2, j] <- mod.tmin$rsq[500]
    
    store.random.forest[impar[i], 3:5, j] <- mod.tmax$importance[, 1]/sum(abs(mod.tmax$importance[, 1]))
    store.random.forest[par[i], 3:5, j] <- mod.tmin$importance[, 1]/sum(abs(mod.tmin$importance[, 1]))

    
  }
}




#### saving and plotting results ####
summarizedresults <- array(NA,  dim = c(10, 9))
rownames(summarizedresults) <- c("Ectotherms CTmax", "Ectotherms CTmin", 
                               "Endotherms CTmax", "Endotherms CTmin", 
                               "Terr.Plants CTmax", "Terr.Plants CTmin", 
                               "Algae CTmax", "Algae CTmin", 
                               "All.Plants CTmax", "All.Plants CTmin")
colnames(summarizedresults) <- c("n", "Rsq", "Rsq.sd", "Age", "Age.sd", "Paleo", 
                                 "Paleo.sd", "Amb.Temp", "Amb.Temp.sd")


summarizedresults[, 1] <- store.random.forest[, 1, 1]
evens <- seq(2,9,2)
alls <- seq(2,5,1)

for(i in 1:4){
  
summarizedresults[, evens[i]] <- apply(store.random.forest[, alls[i], ], 1, 
                                       mean, na.rm = T)
summarizedresults[, evens[i]+1] <- apply(store.random.forest[, alls[i], ], 1, 
                                      sd, na.rm = T)

}



## plotting
dev.off()
par(mfrow = c(3, 1), mar = c(3, 3, 1, 1))
for(j in c(1, 2, 5)){#j = 1
  plot(x = NULL, y = NULL, ylim = c(0, 1), 
       xlim = c(0.5, 3.5), xaxt = "n", yaxt = "n", 
       ylab = "Variable Importance", 
       xlab = "", cex.lab = 1.2)
  axis(1, 1:3, c("Current\ temperature", "Age", "Paleo-\ temperature"))
  axis(2, seq(0, 1, 0.2), seq(0, 1, 0.2), las = 2)
  
  text(0.5, 0.95, c("a", "b", "", "c", "c")[j], cex = 2)
  
  points(c(0.9, 1.9, 2.9), summarizedresults[impar[j], c(8, 4, 6)]/
           sum(summarizedresults[impar[j], c(8, 4, 6)]), pch = 19, 
         col = "red", cex = 1.2)
  points(c(1.1, 2.1, 3.1), summarizedresults[par[j], c(8, 4, 6)]/
           sum(summarizedresults[par[j], c(8, 4, 6)]), pch = 19, 
         col = "blue", cex = 1.2)
  text(3.2, 0.95, paste("R²  =  ", 
                      round(summarizedresults[impar[j], 1], 2)), cex = 1.2, 
       col = "red")
  text(3.2, 0.85, paste("R²  =  ", 
                      round(summarizedresults[par[j], 1], 2)), cex = 1.2, 
       col = "blue")
  
  
  for(i in 1:3){#i = 1
    
    lines(c(rep(i-0.1, 2)), c(summarizedresults[impar[j], c(8, 4, 6)[i]] - 
                                (summarizedresults[impar[j], c(9, 5, 7)[i]])*2, 
                            summarizedresults[impar[j], c(8, 4, 6)[i]] + 
                              (summarizedresults[impar[j], c(9, 5, 7)[i]])*2), 
          col = "red", lwd = 2,  las = 1)
    lines(c(rep(i+0.1, 2)), c(summarizedresults[par[j], c(8, 4, 6)[i]] - 
                                (summarizedresults[par[j], c(9, 5, 7)[i]])*2, 
                            summarizedresults[par[j], c(8, 4, 6)[i]] + 
                              (summarizedresults[par[j], c(9, 5, 7)[i]])*2), 
          col = "blue", lwd = 2,  las = 1)
  }
}






#### end ####







