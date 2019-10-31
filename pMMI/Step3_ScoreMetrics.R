# After you have attempted to model all your metrics, 
# decide which to keep RAW and which to SCALE 

setwd("~/Documents/R/ASCI/pMMI/")

diatoms = F # pick one
sba = F
hybrid = T

if (diatoms == T) {assemblage <- "diatoms"}
if (sba == T) {assemblage <- "sba"}
if (hybrid == T) {assemblage <- "hybrid"}

if (diatoms == T) {setwd("~/Documents/R/ASCI/pMMI/diatoms")}
if (sba == T) {setwd("~/Documents/R/ASCI/pMMI/sba")}
if (hybrid == T) {setwd("~/Documents/R/ASCI/pMMI/hybrid")}

file1 <- paste0(assemblage, ".all.metrics.csv")

all.obs<-read.csv(file1, row.names=1, stringsAsFactors = F) 
all.obs<-all.obs[, colSums(is.na(all.obs)) < nrow(all.obs) * 0.5] #remove cols with > 50% NAs
all.obs<-all.obs[order(row.names(all.obs)),]

stations<-read.csv("~/Documents/R/ASCI/DATA/algae.site.data.07202019.csv", row.names=1)
stations<-subset(stations, row.names(stations) %in% row.names(all.obs))
stations<-stations[order(row.names(stations)),]

row.names(all.obs) == row.names(stations)

rc<-subset(stations, stations$RefCalVal=="Cal")
rv<-subset(stations, stations$RefCalVal=="Val")
str.c<-subset(stations, stations$StrCalVal=="Cal")
str.v<-subset(stations, stations$StrCalVal=="Val")
stations.rc<-subset(stations, row.names(stations) %in% row.names(rc))
stations.rv<-subset(stations, row.names(stations) %in% row.names(rv))

obs.rc<-subset(all.obs, row.names(all.obs) %in% row.names(rc))
obs.int<-subset(all.obs, stations$SiteSetSample2=="Intermediate")
obs.str.c<-subset(all.obs, row.names(all.obs) %in% row.names(str.c))

metricslist<-colnames(all.obs)

all.stations.list<-row.names(all.obs)
rc.list<-row.names(rc)
rv.list<-row.names(rv)
str.c.list<-row.names(str.c)
str.v.list<-row.names(str.v)

#######################################################
#### Determine if metric is increaser or decreaser 
#######################################################
increasers<-list()
decreasers<-list()

for (i in metricslist) { # i<-"cnt.spp.BCG45"
  mean.ref<-mean(obs.rc[[i]], na.rm=T)
  mean.int<-mean(obs.int[[i]], na.rm=T)
  mean.str<-mean(obs.str.c[[i]], na.rm=T)
  if ((mean.ref < mean.str)) {increasers[i]<-i}
  if ((mean.ref > mean.str)) {decreasers[i]<-i}
}

write.csv(increasers, file=paste0(assemblage, ".increasers.txt"))
write.csv(decreasers, file=paste0(assemblage, ".decreasers.txt"))


##########################################################################
#### Let's score all RAW metrics 
##########################################################################

scored.metrics<-data.frame(all.stations.list)
increaser.raw.quants<-list()

for (i in increasers){      # i<-"prop.Cyclotella"
  min<-as.numeric(quantile(obs.rc[[i]], 0.05, na.rm=T)) 
  max<-as.numeric(quantile(obs.str.c[[i]], 0.95, na.rm=T)) # should be str cal ** one sample per site 
  increaser.raw.quants[[i]]<-paste(min, max) 
  observed<-all.obs[[i]]
  
  a <-(observed-max)
  b <-(min-max)
  c <- as.data.frame(a/b)
  names(c)[1]<-i
  scored.metrics<-cbind(scored.metrics,c)
}
write.csv(increaser.raw.quants, paste0(assemblage, ".increaser.raw.quants.csv"))

decreaser.raw.quants<-list()
for (i in decreasers){
  min<-as.numeric(quantile(obs.str.c[[i]], 0.05, na.rm=T))
  max<-as.numeric(quantile(obs.rc[[i]], 0.95, na.rm=T))
  decreaser.raw.quants[[i]]<-paste(min, max)
  observed<-all.obs[[i]]
  
  a <-(observed-min)
  b <-(max-min)
  c <- as.data.frame(a/b)
  names(c)[1]<-i
  scored.metrics<-cbind(scored.metrics,c)
  
}
write.csv(decreaser.raw.quants, paste0(assemblage, ".decreaser.raw.quants.csv"))


write.csv(scored.metrics, paste0(assemblage,".scored.raw.metrics.csv"), row.names=F)

##########################################################################
#### Now let's score the MODELED metrics 
##########################################################################

# make dataframe with all the predicted values for metrics
if(diatoms==T) {all.predz<-read.csv("diatoms.predicted.metrics.csv", row.names=1)}
if(sba==T) {all.predz<-read.csv("sba.predicted.metrics.csv", row.names=1)}
if(hybrid==T) {all.predz<-read.csv("hybrid.predicted.metrics.csv", row.names=1)}

# align pred and obs
foo<-which(row.names(all.obs) %in% row.names(all.predz))
all.obs<-all.obs[foo,]
foo<-which(row.names(all.predz) %in% row.names(all.obs))
all.predz<-all.predz[foo,]
foo<-which(colnames(all.predz) %in% colnames(all.obs))
all.predz<-all.predz[,foo]

metrics<-colnames(all.predz)


##### find the difference between pred and obs 

all.predz<-as.data.frame(all.predz[,metrics]) # dataframe of observed values for modeled metrics 
colnames(all.predz) <- metrics

all.diff<-as.data.frame(row.names(all.obs)) # make an empty dataframe of all diffs between pred/obs

for (i in metrics) { 
diff<- (as.numeric(all.obs[[i]]) - as.numeric(all.predz[[i]])) # changed this 6/21/17
all.diff<-cbind(all.diff, diff)
}

row.names(all.diff)<-all.diff[,1]
all.diff<-data.frame(all.diff[,-1])
colnames(all.diff)<-metrics

write.csv(all.diff, paste0(assemblage,".modeled.residuals.csv"))

##### score the modeled metrics  

metrics 
all.diff 

diff.ref <- subset(all.diff, row.names(all.diff) %in% row.names(obs.rc))
diff.int <- subset(all.diff, row.names(all.diff) %in% row.names(obs.int))
diff.str <- subset(all.diff, row.names(all.diff) %in% row.names(obs.str.c))

scored.mod.metrics<-data.frame(row.names(all.obs)) #empty data frame


increaser.mod.quants<-list()
for (i in increasers){    # i<-"cnt.spp.BCG5"
  if(i %in% colnames(all.diff)) { 
  min<-as.numeric(quantile(diff.ref[[i]], 0.05, na.rm=T)) # try min as zero 
  max<-as.numeric(quantile(diff.str[[i]], 0.95, na.rm=T))
  increaser.mod.quants[[i]]<-paste(min,max)
  observed<-all.diff[[i]] 
  
  a <-(observed-max)
  b <-(min-max)
  c <- as.data.frame(a/b)
  names(c)[1]<-i
  scored.mod.metrics<-cbind(scored.mod.metrics,c)
}}
write.csv(increaser.mod.quants, paste0(assemblage, ".increaser.mod.quants.csv"))


decreaser.mod.quants<-list()
for (i in decreasers){ 
  if(i %in% colnames(all.diff)) { 
  min<-as.numeric(quantile(diff.str[[i]], 0.05, na.rm=T))
  max<-as.numeric(quantile(diff.ref[[i]], 0.95, na.rm=T)) # try max as zero 
  decreaser.mod.quants[[i]]<-paste(min, max)
  observed<-all.diff[[i]]  
  
  a <-(observed-min)
  b <-(max-min)
  c <- as.data.frame(a/b)
  names(c)[1]<-i
  scored.mod.metrics<-cbind(scored.mod.metrics,c)
  }}

write.csv(decreaser.mod.quants, paste0(assemblage, ".decreaser.mod.quants.csv"))

write.csv(scored.mod.metrics, paste0(assemblage,".scored.modeled.metrics.csv"), row.names = F)


 
############################################################################################################
# Now let's try renumbering ------------
############################################################################################################

# scored metrics 
raw<-read.csv(paste0(assemblage,".scored.raw.metrics.csv"), row.names=1)
modeled<-read.csv(paste0(assemblage,".scored.modeled.metrics.csv"), row.names=1)

modeled.renum<-modeled
modeled.renum[modeled.renum < 0] <- 0 
modeled.renum[modeled.renum > 1] <- 1 
modeled.renum$Type<-modeled$Type

raw.renum<-raw
raw.renum[raw.renum < 0] <- 0 
raw.renum[raw.renum > 1] <- 1 
raw.renum$Type<-raw$Type


write.csv(modeled.renum, paste0(assemblage,".scored.modeled.metrics.renum.csv"))
write.csv(raw.renum, paste0(assemblage,".scored.raw.metrics.renum.csv"))

####################################
# SCALE metrics -------------------
####################################

# Combine the raw and modeled combined renumbered metrics 
raw<-read.csv(paste0(assemblage,".scored.raw.metrics.renum.csv"), row.names=1)
modeled<-read.csv(paste0(assemblage,".scored.modeled.metrics.renum.csv"), row.names=1)
row.names(raw) == row.names(modeled)

colnames(raw) <- paste(colnames(raw), "raw", sep = "_")
del<-c("Type", "Type_raw")
foo<-which(colnames(raw)%in%del)
if (length(foo) > 0) {raw<-raw[,-foo]}
foo<-which(colnames(modeled)%in%del)
if (length(foo) > 0) {modeled<-modeled[,-foo]}

combined<-cbind(modeled, raw) #this is the modeled metrics + raw metrics 
combined<-combined[order(row.names(combined)),]

write.csv(combined, "combined.not.scaled.csv")

combined.rc<-subset(combined, row.names(combined) %in% row.names(rc))
combined.str<-subset(combined, row.names(combined) %in% row.names(str.c))

combined.sc <- sweep(combined,MARGIN=2,FUN="/",STATS=colMeans(combined.rc, na.rm = T)) # this should be refcal 7/12/18

refcalmean<-colMeans(combined.rc, na.rm=T)
write.csv(refcalmean, paste0(assemblage, "_refcalmean.csv"))

saveas<-paste0(assemblage, ".combined.scaled.csv")
write.csv(combined.sc, saveas)

foo<-which(row.names(combined.sc) %in% row.names(combined.rc)) 
combined.sc.rc<-combined.sc[foo,]
row.names(combined.sc.rc) == row.names(stations.rc)

foo<-which(row.names(combined.sc) %in% row.names(combined.str))
combined.sc.str<-combined.sc[foo,]

combined.sc.ref<-combined.sc.rc
stations.ref<-stations.rc







