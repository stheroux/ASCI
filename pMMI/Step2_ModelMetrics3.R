# modeling metrics for reference sites 

setwd("~/Documents/R/ASCI/pMMI/")
library(pscl)
library(caret)
library(randomForest)
source("~/Documents/R/bin/OE.load.and.source.R")
source("~/Documents/R/bin/OE.caret.load.and.source.R")
source("~/Documents/R/bin/cbind.fill.R")
source("~/Documents/R/bin/cbind.na.R")

diatoms = F #pick one 
sba = F
hybrid = T

if (diatoms == T) { assemblage <- "diatoms" }
if (sba == T) { assemblage <- "sba" }
if (hybrid == T) { assemblage <- "hybrid" }

# Import metric observed scores 

if (diatoms == T) { 
obs<-read.csv("diatoms/diatoms.all.metrics.csv", row.names=1)
obs<-obs[order(row.names(obs)),]

}


if (sba == T) { 
  obs<-read.csv("sba/sba.all.metrics.csv", row.names=1)
  obs<-obs[order(row.names(obs)),]
  
  }

if (hybrid == T) { 
  obs<-read.csv("hybrid/hybrid.all.metrics.csv", row.names=1)
  obs<-obs[order(row.names(obs)),]
  
  }

# import stations data 
stations<-read.csv("~/Documents/R/ASCI/DATA/algae.site.data.07202019.csv") # most recent + not most recent
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

foo<-which(is.na(colSums(obs.rc)))
if (length(foo) > 0) { obs.rc<-obs.rc[,-foo] }
foo<-which(colSums(obs.rc)==0)
if (length(foo) > 0) { obs.rc<-obs.rc[,-foo] }

##################################################
##### BEGIN LOOP #######
##### this loop will build models to predict metric distributions #######

metrics.list<-colnames(obs.rc)

source("~/Documents/R/ASCI/pMMI/model.loop.R")
source("~/Documents/R/ASCI/pMMI/screen.metrics.R")
library(parallel)

# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores)
base<-model.loop
wd<-"~/Documents/R/ASCI/pMMI/"
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

# put all rf model output in folder .rf.models 
# edit txt file to remove quotes and replace space with tab 

results<-read.table("refcal.lm.out.txt", sep='\t')
colnames(results)<-c("top.metrics", "foo", "rsq")

if (diatoms ==T) { 
write.csv(results, 'diatoms/diatoms.refcal.rsq.csv') }

if (sba ==T) { 
  write.csv(results, 'sba/sba.refcal.rsq.csv') }

if (hybrid ==T) { 
  write.csv(results, 'hybrid/hybrid.refcal.rsq.csv') }

metrics<-results$top.metrics
#results$rsq<-as.numeric(as.character(results$rsq))
tops <-droplevels(subset(results, rsq>0.2)) # or pr2.1>0.1
top.metrics<-tops$top.metrics


if (diatoms == T ) {write.csv(top.metrics, "diatoms/diatoms.top.metrics.csv") }
if (sba == T ) {write.csv(top.metrics, "sba/sba.top.metrics.csv") }
if (hybrid == T ) {write.csv(top.metrics, "hybrid/hybrid.top.metrics.csv") }

##########################################################
#using RF models, predict values 
##########################################################

predicted<-data.frame(row.names(stations))
row.names(predicted)<-predicted[,1]

for (i in metrics) {        # i<-"prop.spp.BCG45"
 
  # import RF model 
  filename<-paste0("~/Documents/R/ASCI/pMMI/", assemblage,".rf.models/",assemblage,".",i,"RF.model.Rdata" )
  load(filename)
  pred.var<-row.names(rf.out$importance)
  
  ### apply new RF model to ALL data 
  stations1<-stations[,pred.var]
  foo<-which(is.na(rowSums(stations1)))
  if (length(foo) > 0) {stations1<-stations1[-foo,]}
  set.seed(10)
  pred<-as.data.frame(predict(rf.out, stations1))
  names(pred)[1]<-i
  predicted<-cbind.na(predicted, pred) } 

predicted<-predicted[,-1]
colnames(predicted)<-metrics

if (diatoms == T ) {write.csv(predicted, "diatoms/diatoms.predicted.metrics.csv") }
if (sba == T ) {write.csv(predicted, "sba/sba.predicted.metrics.csv") }
if (hybrid == T ) {write.csv(predicted, "hybrid/hybrid.predicted.metrics.csv") }

predicted.top<-predicted[,top.metrics]

if (diatoms == T ) {write.csv(predicted.top, "diatoms/diatoms.predicted.top.metrics.csv") }
if (sba == T ) {write.csv(predicted.top, "sba/sba.predicted.top.metrics.csv") }
if (hybrid == T ) {write.csv(predicted.top, "hybrid/hybrid.predicted.top.metrics.csv") }





