#'##############################################################################
#' Script and functions to run:
#' * Clean and format GlobTherm data, climate data, phylogenetic data
#'   
#'
#'  by Ignacio Morales-Castilla
#'  started June 2016
#'##############################################################################


#### to start ####
rm(list=ls())
options(stringsAsFactors=FALSE)


#### load packages ####
packs.to.extract<-list('phytools','geiger','pez','caper','randomForest')
lapply(packs.to.extract,require, character.only=T)
rm(packs.to.extract)



#setwd("~/MEGA/Work_Montreal_postdoc/GTherm/Joanne Project/TTOL data")
setwd("~/GitHub/ThermalEvolution/analyses/")


#'#######################################
#### load tolerance and climate data ####


## Load the GlobTherm dataset (https://www.nature.com/articles/sdata201822). 
## Please cite Bennett et al. (2018) Scientific Data, 5, Article number: 180022  

GTherm.data <- read.csv("../data/GlobTherm.csv")
#dim(GTherm.data)
#str(GTherm.data)


## Load climatic data from: 
## Bio-ORACLE v2.0 (https://www.bio-oracle.org/) 
## Assis, J., et al. (2017) Global Ecology and Biogeography.
## and Worldclim 2 (https://www.worldclim.org/)
## Fick, S.E. and R.J. Hijmans, 2017 International Journal of Climatology 37 (12): 4302-4315. 
## Please note that the following data has been extracted for the GlobTherm
## dataset from the above sources. Make sure they are cited appropriately.

biooracle <- read.table("../data/GlobTherm_coords_biooracle.txt",header = T)

worldclim <- read.table("../data/GlobTherm_coords_worldclim.txt",header = T)

coords <- read.table("../data/GlobTherm_coords_excel.txt",header = T)


#'##########################################
#### merge climate data with GTherm data ###

index.clim <- which(paste(coords$Genus, coords$Species) %in%
                      paste(GTherm.data$Genus, GTherm.data$Species))

matched.clim <- match(paste(coords$Genus[index.clim], coords$Species[index.clim]),
           paste(GTherm.data$Genus, GTherm.data$Species))


GTherm.data$temp.min.marine[matched.clim] <- biooracle$Present.Surface.Temperature.Min[index.clim]
GTherm.data$temp.max.marine[matched.clim] <- biooracle$Present.Surface.Temperature.Max[index.clim]

GTherm.data$temp.min.terrest[matched.clim] <- worldclim$wc2.0_bio_2.5m_11[index.clim]
GTherm.data$temp.max.terrest[matched.clim] <- worldclim$wc2.0_bio_2.5m_10[index.clim]

GTherm.data$temp.min <- ifelse(!is.na(GTherm.data$temp.min.marine),
                               GTherm.data$temp.min.marine, GTherm.data$temp.min.terrest)
GTherm.data$temp.max <- ifelse(!is.na(GTherm.data$temp.max.marine),
                               GTherm.data$temp.max.marine, GTherm.data$temp.max.terrest)

rm(biooracle,worldclim,coords)





#'####################################################
#### clean and complete data with thermy, aquatic ####

## generate classifying variables to subset groupings
GTherm.data$thermy <- ifelse(GTherm.data$Class %in% c("Mammalia", "Aves"), 
                             "endo", "ecto")

GTherm.data$is.aquatic <- ifelse(GTherm.data$Realm %in% c("Marine", "Freshwater", "Marine.Freshwater"), 
                                 "Aquatic",ifelse(GTherm.data$Realm %in% "Intertidal" , 
                                                  NA, "Terrestrial"))

GTherm.data$Palaeotemp <- ifelse(GTherm.data$order.temp.origin == 5, 3, 
                              ifelse(GTherm.data$order.temp.origin == 4, 2, 1))

GTherm.data$warm.cold <- ifelse(GTherm.data$Palaeotemp == 1, 1, 2) 

GTherm.data$Plant <- ifelse(GTherm.data$Phylum %in% c("Rhodophyta","Streptophyta","Chlorophyta","Phaeophyceae"),
                            "Plants", "No.plants")

GTherm.data$Plant.only <- ifelse(GTherm.data$Phylum == "Streptophyta", 
                                "Plants.only", "No.plants")

ectolandwarmo <- subset(GTherm.data, GTherm.data$thermy=="ecto" & 
                        GTherm.data$Realm=="Terrestrial" & 
                        GTherm.data$order.temp.origin==5)



#'######################
#### load phylogeny ####

## load phylogeny from: https://academic.oup.com/mbe/article/32/4/835/1078218
## Hedges, S. B., et al. (2015) Molecular biology and evolution, 32(4), 835-845.
## Please note that the following phylogeny has been pruned for the GlobTherm
## dataset from the above source. Make sure it is cited appropriately.

ttol <- read.tree("../data/Phylo_GlobTherm.tre")

names.inTTOL <- ttol$tip.label


## Generate vector of species names
GTherm.data$species.phylo <- paste(GTherm.data$Genus,GTherm.data$Species,sep="_")


## prunning phylogeny by species names (binomial)
sps.to.prune <- which(!names.inTTOL %in% GTherm.data$species.phylo)

ttol.GTherm.phy <- drop.tip(ttol,sps.to.prune)
rm(ttol,index.clim,matched.clim,sps.to.prune)



## end
