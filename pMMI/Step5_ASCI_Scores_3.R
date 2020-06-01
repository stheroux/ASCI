# Step5_ASCI_SCORES

setwd("~/Documents/R/ASCI/pMMI_newnew/")

stations<-read.csv("~/Documents/R/ASCI/DATA/algae.site.data.04052019.csv", row.names=1)
source("~/Documents/R/bin/cbind.na.R")

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
pMMI.d<-read.csv("~/Documents/R/ASCI/pMMI_newnew/diatoms/diatoms.combined.win.metrics.csv", header=T, row.names=1)
pMMI.sba<-read.csv("~/Documents/R/ASCI/pMMI_newnew/sba/sba.combined.win.metrics.csv", header=T, row.names=1)
pMMI.hybrid<-read.csv("~/Documents/R/ASCI/pMMI_newnew/hybrid/hybrid.combined.win.metrics.csv", header=T, row.names=1)

pMMI.d.null<-read.csv("~/Documents/R/ASCI/pMMI_newnew/diatoms/diatoms.combined.win.metrics.null.csv", header=T, row.names=1)
pMMI.sba.null<-read.csv("~/Documents/R/ASCI/pMMI_newnew/sba/sba.combined.win.metrics.null.csv", header=T, row.names=1)
pMMI.hybrid.null<-read.csv("~/Documents/R/ASCI/pMMI_newnew/hybrid/hybrid.combined.win.metrics.null.csv", header=T, row.names=1)

colnames(pMMI.d.null)<-paste0(colnames(pMMI.d.null), ".null")
colnames(pMMI.sba.null)<-paste0(colnames(pMMI.sba.null), ".null")
colnames(pMMI.hybrid.null)<-paste0(colnames(pMMI.hybrid.null), ".null")

#pMMI.d<-as.data.frame(cbind.na(pMMI.d, pMMI.d.null))
#pMMI.sba<-as.data.frame(cbind.na(pMMI.sba, pMMI.sba.null))
#pMMI.hybrid<-as.data.frame(cbind.na(pMMI.hybrid, pMMI.hybrid.null))

# both 

d1<-data.frame(row.names(OE.d),OE.d$OoverE, OE.d$OoverE.null)
d2<-data.frame(row.names(pMMI.d),pMMI.d$MMI_scaled, pMMI.d.null$MMI_scaled.null)
names(d1)<-c("SampleID", "OoverE", "OoverE.null")
names(d2)<-c("SampleID", "MMI_scaled", "MMI_scaled.null")
d<-merge(d1,d2, by="SampleID")

sba1<-data.frame(row.names(OE.sba),OE.sba$OoverE, OE.sba$OoverE.null)
sba2<-data.frame(row.names(pMMI.sba),pMMI.sba$MMI_scaled, pMMI.sba.null$MMI_scaled.null)
names(sba1)<-c("SampleID", "OoverE", "OoverE.null")
names(sba2)<-c("SampleID", "MMI_scaled", "MMI_scaled.null")
sba<-merge(sba1,sba2, by="SampleID")

hybrid1<-data.frame(row.names(OE.hybrid),OE.hybrid$OoverE, OE.hybrid$OoverE.null)
hybrid2<-data.frame(row.names(pMMI.hybrid),pMMI.hybrid$MMI_scaled, pMMI.hybrid.null$MMI_scaled.null)
names(hybrid1)<-c("SampleID", "OoverE", "OoverE.null")
names(hybrid2)<-c("SampleID", "MMI_scaled", "MMI_scaled.null")
hybrid<-merge(hybrid1,hybrid2, by="SampleID")

##### 
colnames(d)<-paste0(colnames(d), ".d")
colnames(sba)<-paste0(colnames(sba), ".sba")
colnames(hybrid)<-paste0(colnames(hybrid), ".hybrid")

colnames(d)[1]<-"SampleID"
colnames(sba)[1]<-"SampleID"
colnames(hybrid)[1]<-"SampleID"

asci.scores<-merge(d,sba, by="SampleID", all.x=T, all.y=T)
asci.scores<-merge(asci.scores,hybrid, by="SampleID", all.x=T, all.y=T)

asci.scores$MMI.hybrid<-asci.scores$MMI_scaled.hybrid
asci.scores$MMI.d<-asci.scores$MMI_scaled.d
asci.scores$MMI.sba<-asci.scores$MMI_scaled.sba

asci.scores$MMI.null.hybrid<-asci.scores$MMI_scaled.null.hybrid
asci.scores$MMI.null.d<-asci.scores$MMI_scaled.null.d
asci.scores$MMI.null.sba<-asci.scores$MMI_scaled.null.sba



write.csv(asci.scores, "asci.scores.csv")


#################

stations<-read.csv("~/Documents/R/ASCI/DATA/algae.site.data.04052019.csv")

asci.scores$SampleID_old<-asci.scores$SampleID

asci.scores.forRafi<-merge(asci.scores, stations, by="SampleID_old", all.x=T)

write.csv(asci.scores.forRafi, "asci.scores.forRafi.csv")
names(asci.scores.forRafi) <- gsub(x = names(asci.scores.forRafi), pattern = "\\Means", replacement = "MMI")  
write.csv(asci.scores.forRafi, "asci.scores.forRafi2.csv")

