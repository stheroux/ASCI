# modeling metrics for reference sites 

setwd(paste0("/Users/susannat/Documents/R/ASCI/pMMI_newnew/",assemblage))


# Import metric observed scores 

if (diatoms == T) { 
obs<-read.csv("diatoms.all.metrics.csv", row.names=1)
obs<-obs[order(row.names(obs)),]
dir.create("~/Documents/R/ASCI/pMMI_newnew/diatoms/diatoms.rf.models/")
  }

if (sba == T) { 
  obs<-read.csv("sba.all.metrics.csv", row.names=1)
  obs<-obs[order(row.names(obs)),]
  dir.create("~/Documents/R/ASCI/pMMI_newnew/sba/sba.rf.models/")
  }

if (hybrid == T) { 
  obs<-read.csv("hybrid.all.metrics.csv", row.names=1)
  obs<-obs[order(row.names(obs)),]
  dir.create("~/Documents/R/ASCI/pMMI_newnew/hybrid/hybrid.rf.models/")
  }

# import stations data 
stations<-read.csv("~/Documents/R/ASCI/DATA/algae.site.data.04052019.csv") # most recent + not most recent
row.names(stations)<-stations$SampleID_old
stations<-stations[order(row.names(stations)),]

source("~/Documents/R/ASCI/DATA/cand.var.bu.R")
cand.var
cand.var %in% names(stations) #do you have them all?

cal<-subset(stations, stations$RefCalVal=="Cal")
val<-subset(stations, stations$RefCalVal=="Val")

stations.rc<-subset(stations, row.names(stations) %in% row.names(cal))
stations.rv<-subset(stations, row.names(stations) %in% row.names(val))
obs.rc<-subset(obs, row.names(obs) %in% row.names(cal))
obs.rv<-subset(obs, row.names(obs) %in% row.names(val))

stations.rc<-droplevels(subset(stations.rc, row.names(stations.rc) %in% row.names(obs.rc)))

row.names(stations.rc) == row.names(obs.rc)

#foo<-which(is.na(colSums(obs.rc)))
#if (length(foo) > 0) { obs.rc<-obs.rc[,-foo] }
#foo<-which(colSums(obs.rc)==0)
#if (length(foo) > 0) { obs.rc<-obs.rc[,-foo] }


metrics.list<-sort(colnames(obs.rc))

##################################################
##### BEGIN LOOP #######
##### this loop will build models to predict metric distributions #######


source("~/Documents/R/ASCI/pMMI_newnew/model.loop.R")
source("~/Documents/R/ASCI/pMMI_newnew/screen.metrics.R")
library(parallel)

# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores)
base<-model.loop
wd<-"~/Documents/R/ASCI/pMMI_newnew/"
clusterExport(cl, c("base",
                    "metrics.list", 
                    "obs.rc", "stations.rc", 
                    "ctrl", "assemblage", "cand.var", "wd" ))

parSapply(cl, metrics.list,
          function(x)
            fun=base(x))

stopCluster(cl)


##### END LOOP #######
######################### 




