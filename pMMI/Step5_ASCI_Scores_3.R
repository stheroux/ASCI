# Step5_ASCI_SCORES

setwd("~/Documents/R/ASCI/pMMI/")
source("~/Documents/R/bin/cbind.na.R")
library(ggplot2)

stations<-read.csv("~/Documents/R/ASCI/DATA/algae.site.data.07202019.csv", row.names=1)

# O/E 
OE.d<-read.csv("~/Documents/R/ASCI/OE/diatoms.best/Scores.all.csv", header=T, row.names=1)
OE.d<-OE.d[,c(1:3,6)]
names(OE.d) <- gsub(x = names(OE.d), pattern = "\\OE.scores.", replacement = "")  

OE.sba<-read.csv("~/Documents/R/ASCI/OE/sba.best/Scores.all.csv", header=T, row.names=1)
OE.sba<-OE.sba[,c(1:3,6)]
names(OE.sba) <- gsub(x = names(OE.sba), pattern = "\\OE.scores.", replacement = "")  

OE.hybrid<-read.csv("~/Documents/R/ASCI/OE/hybrid.best/Scores.all.csv", header=T, row.names=1)
OE.hybrid<-OE.hybrid[,c(1:3,6)]
names(OE.hybrid) <- gsub(x = names(OE.hybrid), pattern = "\\OE.scores.", replacement = "")  


# MMI
pMMI.d<-read.csv("~/Documents/R/ASCI/pMMI/diatoms/diatoms.combined.win.metrics.scores.csv", header=T, row.names=1)
pMMI.sba<-read.csv("~/Documents/R/ASCI/pMMI/sba/sba.combined.win.metrics.scores.csv", header=T, row.names=1)
pMMI.hybrid<-read.csv("~/Documents/R/ASCI/pMMI/hybrid/hybrid.combined.win.metrics.scores.csv", header=T, row.names=1)

pMMI.d.null<-read.csv("~/Documents/R/ASCI/pMMI/diatoms/diatoms.combined.win.metrics.scores.null.csv", header=T, row.names=1)
pMMI.sba.null<-read.csv("~/Documents/R/ASCI/pMMI/sba/sba.combined.win.metrics.scores.null.csv", header=T, row.names=1)
pMMI.hybrid.null<-read.csv("~/Documents/R/ASCI/pMMI/hybrid/hybrid.combined.win.metrics.scores.null.csv", header=T, row.names=1)

colnames(pMMI.d.null)<-paste0(colnames(pMMI.d.null), ".null")
colnames(pMMI.sba.null)<-paste0(colnames(pMMI.sba.null), ".null")
colnames(pMMI.hybrid.null)<-paste0(colnames(pMMI.hybrid.null), ".null")

pMMI.d<-as.data.frame(cbind.na(pMMI.d, pMMI.d.null))
pMMI.sba<-as.data.frame(cbind.na(pMMI.sba, pMMI.sba.null))
pMMI.hybrid<-as.data.frame(cbind.na(pMMI.hybrid, pMMI.hybrid.null))

# both 

d<-cbind.na(OE.d[,3:4], pMMI.d[,c(2,4)])
sba<-cbind.na(OE.sba[,3:4], pMMI.sba[,c(2,4)])
hybrid<-cbind.na(OE.hybrid[,3:4], pMMI.hybrid[,c(2,4)])

##### 
colnames(d)<-paste0(colnames(d), ".d")
colnames(sba)<-paste0(colnames(sba), ".sba")
colnames(hybrid)<-paste0(colnames(hybrid), ".hybrid")

site.list<-c(row.names(d), row.names(hybrid), row.names(sba))
site.list<-unique(site.list)

asci.scores<-as.data.frame(site.list)
row.names(asci.scores)<-asci.scores$site.list
  
asci.scores<-cbind.na(asci.scores,d)
asci.scores<-cbind.na(asci.scores, sba)
asci.scores<-cbind.na(asci.scores, hybrid)


colnames(asci.scores)

write.csv(asci.scores, "asci.scores.csv")


#################

setwd("~/Documents/R/ASCI/pMMI/")

source("~/Documents/R/bin/cbind.na.R")

asci.scores<-read.csv("asci.scores.csv", header=T, row.names=1)
asci.scores<-asci.scores[order(row.names(asci.scores)),]
asci.scores$zed<-row.names(asci.scores)


stations<-read.csv("~/Documents/R/ASCI/DATA/algae.site.data.07202019.csv", header=T, row.names=1)
stations<-stations[order(row.names(stations)),]
stations$zed2<-row.names(stations)
stations<-subset(stations, row.names(stations) %in% row.names(asci.scores))
row.names(asci.scores) == row.names(stations)

asci.scores.forRafi<-cbind.na(stations, asci.scores)
write.csv(asci.scores.forRafi, "asci.scores.forRafi.csv")
names(asci.scores.forRafi) <- gsub(x = names(asci.scores.forRafi), pattern = "\\Means", replacement = "MMI")  
write.csv(asci.scores.forRafi, "asci.scores.forRafi2.csv")

