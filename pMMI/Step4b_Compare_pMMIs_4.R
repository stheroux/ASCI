# compare MMI

setwd("~/Documents/R/ASCI/pMMI/")
source("~/Documents/R/bin/cbind.na.R")
library(plyr); library(tidyverse)

diatoms = F # pick one 
sba = F
hybrid = T

if (diatoms == T) {assemblage <- "diatoms"}
if (sba == T) {assemblage <- "sba"}
if (hybrid == T) {assemblage <- "hybrid"}

#################################
# Import data ------------------
#################################
if (diatoms == T) {setwd("~/Documents/R/ASCI/pMMI/diatoms/")}
if (sba == T) {setwd("~/Documents/R/ASCI/pMMI/sba/")}
if (hybrid == T) {setwd("~/Documents/R/ASCI/pMMI/hybrid/")}

# import scaled metrics 
filename<-paste0(assemblage, ".combined.scaled.csv")
combined.sc<-read.csv(filename, row.names=1)

# import stations 
stations<-read.csv('~/Documents/R/ASCI/DATA/algae.site.data.07202019.csv')
row.names(stations)<-stations$SampleID_old
stations<-stations[order(row.names(stations)),]
stations<-subset(stations, row.names(stations) %in% row.names(combined.sc), drop=T)
foo<-sapply(stations, is.numeric)
stations.num<-stations[,foo] # only numeric columns 

row.names(combined.sc) == row.names(stations)
row.names(combined.sc) == row.names(stations.num)

calval<-subset(stations, stations$RefCalVal!="")

# subset 
combined.sc.ref<-subset(combined.sc, row.names(combined.sc) %in% row.names(calval))
stations.ref<-subset(stations, row.names(stations) %in% row.names(calval))

stations.str<-subset(stations, stations$SiteSetSample2=="Stressed")
combined.sc.str<-subset(combined.sc, row.names(combined.sc) %in% row.names(stations.str))
row.names(stations.str) == row.names(combined.sc.str)

metrics.list<-colnames(combined.sc)

if (assemblage == "diatoms") {results3<-read.csv("diatoms.results3.csv", header = T, row.names=1) }
if (assemblage == "sba") {results3<-read.csv("sba.results3.csv", header = T, row.names=1) }
if (assemblage == "hybrid") {results3<-read.csv("hybrid.results3.csv", header = T, row.names=1) }

rc<-subset(calval, calval$RefCalVal=="Cal")
rv<-subset(calval, calval$RefCalVal=="Val")



