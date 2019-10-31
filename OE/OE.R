# Build an O/E model

#Set a working directory
setwd("~/Documents/R/ASCI/OE/")

diatoms = F  # pick 1
sba =    F #
combo =   T # 

#### KNOBS TO TURN 
prezabs = T
genuslev = T
plotz = T
exportfiles = T
remlowtaxa = F # default F
trySOM = F # default F

cutoff = T  # deletes diatom samples with <200 valve counts , default T 
removeNA = T # deletes taxa with the name 'NA'
removeQual = T # remove qual sample
removeplanktonic = F # remove planktonic species


removeoutliersites=F # default F

if (diatoms==T) { 
#predictornumz<-c(10) ; clusternumz<-c(4:15) ; Pcapz <- c(0.4,0.5) } 
predictornumz<-c(10) ; clusternumz<-c(6) ; Pcapz <- c(0.5) } 

if (sba==T) {
  #predictornumz<-c(10) ; clusternumz<-c(4:15) ; Pcapz <- c(0.4,0.5) }
  predictornumz<-c(10) ; clusternumz<-c(8) ; Pcapz <- c(0.4) }
 
if (combo==T) { 
  #predictornumz<-c(10) ; clusternumz<-c(4:15) ; Pcapz <- c(0.4,0.5) }
  predictornumz<-c(10) ; clusternumz<-c(15) ; Pcapz <- c(0.4) }
  
#### Load packages 
library("sampling")
library("rms")
library("plyr")
library("ggplot2")
library("gridExtra")
library(cluster)
library(fpc)
source("~/Documents/R/bin/OE.load.and.source.R")
source("~/Documents/R/bin/OE.cand.vars.R")
source("~/Documents/R/bin/OE.caret.load.and.source.R")
source("~/Documents/R/bin/model.predict.RanFor.4.2_ed.r") #overwrites earlier model.predict

#export order 
#Scores.cal.rfe
#Scores.rv
#Scores.new.int
#Scores.new.str

#Import your bug data
bugs<-read.csv("~/Documents/R/ASCI/DATA/algae.bug.data.10172019.csv")
bugs<-subset(bugs,  ProtocolCode!="SNARL_2008_WS")
bugs<-droplevels(bugs)
bugs.orig<-bugs



#Import stations data, make sure #N/A are NA
stations<-read.csv("~/Documents/R/ASCI/DATA/algae.site.data.07202019.csv", stringsAsFactors = F)
row.names(stations)<-stations$SampleID_old
stations.orig<-stations

stations<-subset(stations, PSA6c!="NA", drop=T)
stations<-subset(stations, PSA6c!="NV", drop=T)
stations<-droplevels(stations)
stations<-subset(stations, MostRecentSample=="y", drop=T)

if (removeoutliersites==T) {source("remoutliers.R") 
  foo<-which(row.names(stations) %in% rem)
  stations<-stations[-foo,]
  }

###################################################################################
# ASSIGN CAL/VAL -------------------------------------------------
###################################################################################

stations.rc<-subset(stations, RefCalVal=="Cal")
stations.rv<-subset(stations, RefCalVal=="Val")
stations.r<-subset(stations, RefCalVal!="NA")

# 
# stations.r<-subset(stations, SiteSetSample2=="Reference", drop=T)
# df<-stations.r
# 
# set.seed(1)   # for reproducible example
# for (i in unique(df$PSA6c)) {
#   print(i)
#   dfi<-subset(df,PSA6c==i)
#   len <- length(dfi$PSA6c)
#   dfi$CalVal<-sample(x=(c("Cal","Val")),size=len,replace=T,p=c(0.8,0.2)) 
#   dfiname <- paste("out.dfi", i, sep="")
#   assign(dfiname,dfi) }
# 
# dfi.cat<-rbind(out.dfiNC,out.dfiCH,out.dfiSN,out.dfiSC,
#                out.dfiCV,out.dfiDM)
# 
# sum(dfi.cat$CalVal=="Cal")
# sum(dfi.cat$CalVal=="Val")
# 
# stations.r<-dfi.cat
# stations.rc<-subset(stations.r, CalVal=="Cal")
# stations.rv<-subset(stations.r, CalVal=="Val")
# 
# write.csv(stations.r, "Stations.r.calval.csv")
# 
# rm(df, out.dfiNC,out.dfiCH,out.dfiSN,out.dfiSC,
#    out.dfiCV,out.dfiDM, dfi, dfi.cat)
# 

# for stressed 
#stations.str<-subset(stations, SiteSetSample2=="Stressed", drop=T)
#df<-stations.str

#set.seed(1)   # for reproducible example
#for (i in unique(df$PSA6c)) {
#  print(i)
#  dfi<-subset(df,PSA6c==i)
#  len <- length(dfi$PSA6c)
#  dfi$CalVal<-sample(x=(c("Cal","Val")),size=len,replace=T,p=c(0.8,0.2)) 
#  dfiname <- paste("out.dfi", i, sep="")
#  assign(dfiname,dfi) }

#dfi.cat<-rbind(out.dfiNC,out.dfiCH,out.dfiSN,out.dfiSC,
#               out.dfiCV,out.dfiDM)

#sum(dfi.cat$CalVal=="Cal")
#sum(dfi.cat$CalVal=="Val")

#stations.str<-dfi.cat
#stations.str.c<-subset(stations.str, CalVal=="Cal")
#stations.str.v<-subset(stations.str, CalVal=="Val")

#write.csv(stations.str, "Stations.str.calval.csv")

#rm(df, out.dfiNC,out.dfiCH,out.dfiSN,out.dfiSC,
#   out.dfiCV,out.dfiDM, dfi, dfi.cat)



###################################################################################
# Step 1 DATA PREP -------------------------------------------------
###################################################################################

# delete poorly identified taxa -- how does this impact total vols? 
bugs<-subset(bugs, Phylum!="NA")

if (removeQual == T) { 
  bugs<-subset(bugs, SampleTypeCode!="Qualitative", drop=T) } 
if (removeQual == F ) { 
 bugs<- within(bugs, Result[SampleTypeCode == 'Qualitative'] <- 1) }

# delete low count diatom samples
if (cutoff==T) { 
  if (diatoms==T) {  
    bugs<-subset(bugs, Phylum=="Bacillariophyta") 
    bugs<-subset(bugs, BAResult!="NULL")
    bugs<-subset(bugs, BAResult!=0) 
    bugs.m0<-acast(bugs,SampleID_old~FinalIDassigned, value.var="BAResult", fun.aggregate=sum )
    foo<-which(rowSums(bugs.m0)<200)
    rownames.to.delete<-row.names(bugs.m0[foo,])
    foo<-which(bugs$SampleID_old %in% rownames.to.delete)
    bugs<-bugs[-foo,]
  }}

if (diatoms==T) { 
  bugs<-subset(bugs, Phylum=="Bacillariophyta") 
  bugs<-subset(bugs, BAResult!="NULL")
  bugs<-subset(bugs, BAResult!=0)
  
  #remove low number taxa in each sample
  if (remlowtaxa ==T) { bugs<-subset(bugs, bugs$BAResult>5) }
  
  #set.seed(10)
  #bugs.sub<-rarify(inbug=bugs, sample.ID="SampleID_old", abund="BAResult", subsiz=500) 
  bugs.sub<-bugs
  if (genuslev==T) {bugs.m<-acast(bugs.sub, SampleID_old~Genus, value.var="BAResult", fun.aggregate=sum)}
  if (genuslev==F) {bugs.m<-acast(bugs.sub, SampleID_old~FinalIDassigned, value.var="BAResult", fun.aggregate=sum)}
  bugs.m.sp<-acast(bugs.sub, SampleID_old~FinalIDassigned, value.var="BAResult", fun.aggregate=sum)
  write.csv(bugs.m.sp, "diatoms.m.species.csv")
  
  }

if (sba==T) { 
  bugs<-subset(bugs, Phylum!="Bacillariophyta") 
  bugs<-subset(bugs, Result!="NULL") # aug 2019 commented out 
  bugs<-subset(bugs, Result!=0) # aug 2019 commented out 
  bugs<-subset(bugs, Result!="NA") # aug 2019 commented out 
  #bugs$ComboResult<-as.numeric(pmax(bugs$BAResult,bugs$Result, na.rm=T))
  #bugs<-subset(bugs, ComboResult!="NULL", drop = T)
  #bugs<-subset(bugs, ComboResult!=0, drop=T)
  if (genuslev==T) {bugs.m<-acast(bugs, SampleID_old~Genus, value.var="Result", fun.aggregate=sum)} 
  if (genuslev==F) {bugs.m<-acast(bugs, SampleID_old~FinalIDassigned, value.var="Result", fun.aggregate=sum)}
  bugs.m.sp<-acast(bugs, SampleID_old~FinalIDassigned, value.var="Result", fun.aggregate=sum)
  write.csv(bugs.m.sp, "sba.m.species.csv")
  
  }

if (combo==T) { 
  bugs$ComboResult<-as.numeric(pmax(bugs$BAResult,bugs$Result, na.rm=T))
  bugs<-subset(bugs, ComboResult!="NULL", drop = T)
  bugs<-subset(bugs, ComboResult!=0, drop=T)
  if (genuslev==T) {bugs.m<-acast(bugs, SampleID_old~Genus, value.var="ComboResult", fun.aggregate=sum) }
  if (genuslev==F) {bugs.m<-acast(bugs, SampleID_old~FinalIDassigned, value.var="ComboResult", fun.aggregate=sum) }
  bugs.m.sp<-acast(bugs, SampleID_old~FinalIDassigned, value.var="ComboResult", fun.aggregate=sum)
  write.csv(bugs.m.sp, "hybrid.m.species.csv")
  
  }

if (removeNA == T) { 
del<-"NA"
bugs.m<-as.data.frame(bugs.m)
foo<-which(colnames(bugs.m) %in% del)
if (length(foo) >0) { bugs.m<-bugs.m[,-foo] } 
 }

if (removeplanktonic == T) {
  traits<-read.csv("~/Documents/R/ASCI/DATA/algae.traits.data.02102018_3.csv")
  benthic<-subset(traits, traits$Habitat_genus=="BN")
  bugs.m<-subset(bugs.m, colnames(bugs.m) %in% benthic$FinalIDassigned)
}


#subset bug matrices to just the refcal sites
bugs.rc<-as.data.frame(subset(bugs.m, row.names(bugs.m) %in% row.names(stations.rc), drop = T))
foo<-which(colSums(bugs.rc)==0) # no zero sum columns
bugs.rc<-bugs.rc[,-foo]

#subset stations.rc matrices to just the bug sites
stations.rc<-subset(stations.rc, row.names(stations.rc) %in% rownames(bugs.rc))

#Check to make sure rows are aligned correctly
row.names(bugs.rc)==row.names(stations.rc)

#If they aren't, you can do this:
bugs.rc<-bugs.rc[row.names(stations.rc),]
row.names(bugs.rc)==row.names(stations.rc)

#Make a presence-absence version
if (prezabs==T) {bugs.rc<-ifelse(bugs.rc>0,1,0)}
if (prezabs==F) {bugs.rc<-bugs.rc}
bugs.rc<-as.data.frame(bugs.rc)

# export matrices 
if (diatoms==T) {write.csv(bugs.m, "diatoms.bugs.m.csv")}
if (sba==T) {write.csv(bugs.m, "sba.bugs.m.csv")}
if (combo==T) {write.csv(bugs.m, "hybrid.bugs.m.csv")}

bugs.rc.norare<-bugs.rc
bugs.rc.norare<-ifelse(bugs.rc.norare>0,1,0) #always cluster on a PA 

# lower threshold 
thresh<-(nrow(bugs.rc)*0.025)
foo<-which(colSums(bugs.rc.norare) < thresh)
bugs.rc.norare<-as.data.frame(bugs.rc.norare[,-foo])

# upper threshold 
thresh<-(nrow(bugs.rc)*0.95)
foo<-which(colSums(bugs.rc.norare) > thresh)
if (length(foo) > 0) { bugs.rc.norare<-as.data.frame(bugs.rc.norare[,-foo]) } 

foo<-which(rowSums(bugs.rc.norare)==0)
if (length(foo)>0) {bugs.rc.norare<-bugs.rc.norare[-foo,]}

stations.rc<-stations.rc[which(row.names(stations.rc) %in% row.names(bugs.rc.norare)),]
bugs.rc<-bugs.rc[which(row.names(bugs.rc) %in% row.names(stations.rc)),]

foo<-paste("prezabs", prezabs, "genuslev", genuslev, "remlowtaxa", remlowtaxa, "cutoff", cutoff, 
           "removeNA", removeNA, "removeQual", removeQual, "removeplanktonic", removeplanktonic)

write.csv(foo, "settings.txt")

#################################################
####### Start Loop ----- ***********
#################################################

for (x in predictornumz) { 
  for (y in clusternumz) { 
    for (p in Pcapz) {
    predictornum<-x
    clusternum<-y
    Pcap<-p
    
##################################################################;
# CLUSTERING -----------------------------------
##################################################################;

#Recommend doing clustering on data sets with rare taxa removed. 
#But build models using the data set that includes rare taxa
#https://www.r-statistics.com/2013/08/k-means-clustering-from-r-in-action/
#fit.km <- kmeans(bugs.rc,clusternum, nstart=25)     
#stations.rc$BG<-as.factor(fit.km$cluster)
        
set.seed(10)
fit.km <- kmeans(bugs.rc.norare,clusternum, nstart=25)     
stations.rc$BG<-as.factor(fit.km$cluster)
saveas<-paste0("Prednums_",x,"_Clusternums_",y,"_Pcap_", p, ".BG.csv")
write.csv(data.frame(row.names(stations.rc), stations.rc$BG), saveas)

saveas<-paste0("Prednums_",x,"_Clusternums_",y,"_Pcap_", p, "clusterplot.pdf")
pdf(saveas,width=6,height=4) 
plotcluster(bugs.rc.norare, fit.km$cluster)
dev.off()


#clusplot(bugs.rc.norare, fit.km$cluster, color=TRUE, shade=TRUE, 
#         labels=2, lines=0)

saveas<-paste0("Prednums_",x,"_Clusternums_",y,"_Pcap_", p, "_anosim")
set.seed(10)
anosim1<-anosim(bugs.rc.norare, stations.rc$BG, permutations = 99, distance="bray")
anosim.out<-paste(saveas, "R", anosim1$statistic, "Sig", anosim1$signif) 
capture.output(anosim.out, file="anosim.out.txt", append=T)


if (trySOM==T) { 
  clustass<-read.csv("clustass.csv", row.names=1)
  clustass<-as.data.frame(cbind(rownames(stations.rc), clustass$x))
  clustass<-droplevels(clustass)
  stations.rc$BG<-as.factor(clustass$V2)
  }


###################################################################################
# EVALUATE CANDIDATE VARS ----------
###################################################################################

#import from file above
source("~/Documents/R/ASCI/DATA/cand.var.bu.R")
cand.var
cand.var %in% names(stations.rc) #do you have them all?

#############################
# CARET analysis to select env vars -------
#############################

# prep CARET dataframe 
dat<-data.frame(stations.rc[,cand.var])
dat.BG<-stations.rc$BG

set.seed(10)
rfe.out<-rfe(y=dat.BG, x=dat, sizes=2:predictornum,  rfeControl=ctrl, maximize=T)

#record optvars
optvars.rfe <- rfe.out$optVariables
optvars.rfe

#view results
#ggplot(data=rfe.out$results, aes(x=Accuracy, y=Kappa))+
 # geom_text(aes(label=Variables))

##########################################################
# RF with env vars from CARET -----------------
##########################################################

#from RFE
set.seed(10)
rfe.out2<-randomForest(x = dat[,optvars.rfe], y = dat.BG, 
                      ntree=5000, importance=TRUE, norm.votes=TRUE, keep.forest=TRUE)
print(rfe.out2)


if (diatoms == T) { 
diatom.rf.oe<-rfe.out2
save(diatom.rf.oe, file="diatom.RF.OE.Rdata") 
write.csv(row.names(diatom.rf.oe$importance), "diatoms.predvars.txt")
}

if (sba == T)  { 
sba.rf.oe<-rfe.out2
save(sba.rf.oe, file="sba.RF.OE.Rdata")
write.csv(row.names(sba.rf.oe$importance), "sba.predvars.txt")}

if (combo == T) { 
hybrid.rf.oe<-rfe.out2
save(hybrid.rf.oe, file="hybrid.RF.OE.Rdata") 
write.csv(row.names(hybrid.rf.oe$importance), "hybrid.predvars.txt")}

############################################
# Step 9 Calculate scores -----------------
############################################

# need dataframes to be numeric
if (prezabs==F) {bugs.rc<-ifelse(bugs.rc>0,1,0)} 

# save opt vars
saveas<-paste0("PredNums_",x,"_Clustnum_",y,"_Pcap_", p, "_RFEoptvars", ".txt")
foo<-paste(saveas,unlist(optvars.rfe))
capture.output(foo, file="optvars.rfe.all.csv")

# save RFE scores
set.seed(10)
Scores.cal<-model.predict.RanFor.4.2(
  bugcal.pa = bugs.rc,
  grps.final = dat.BG,
  preds.final = optvars.rfe, # for caret
  ranfor.mod = rfe.out2, # for caret
  prednew = stations.rc,
  bugnew = (bugs.rc),
  Pc=Pcap,
  Cal.OOB=TRUE);


# Accuracy: Intercept should be zero
# Precision: Points should be tightly clustered around the regression line
saveas2<-paste0("Prednums_",x,"_Clusternums_",y,"_Pcap_", p,".pdf")
ggplot(data=Scores.cal$OE.scores, aes(x=E, y=O), cex=2)+
  geom_point() + theme_bw() + xlab("Expected") + ylab("Observed") + 
  geom_smooth(method=lm, color="red") #Observed prediction
if (plotz == T) { ggsave(file=saveas2, width = 6, height = 6) }

#plot the O/E for ref sites by region 
scores.out<-as.data.frame(cbind(Scores.cal$OE.scores, stations.rc$PSA6c))
saveas2<-paste0("Prednums_",x,"_Clusternums_",y,"_Pcap_", p,"PSAregion.pdf")
ggplot(scores.out, aes(y=as.numeric(scores.out$OoverE))) + 
  geom_boxplot(aes(x=factor(scores.out$`stations.rc$PSA6c`))) +
  geom_jitter(aes(x=factor(scores.out$`stations.rc$PSA6c`)), cex=0.5) +
  geom_hline(yintercept = 1) + theme_bw() + 
  xlab("") + ylab("O/E")
if (plotz == T) { ggsave(file=saveas2, width = 6, height = 6) }

# Calculate replicate sampling SD of O/E
#Then execute the function, using either in-bag or OOB predicted occurrence probs ('Capture probs') for the calibration data;
saveas2<-paste0("Prednums_",x,"_Clusternums_",y,"_Pcap_",p,"_repsamsd")
repsamsd1<-rep.sam.sd(occprb.cal=Scores.cal$Capture.Probs,Pc=Pcap);
repsamsd1.out<-paste(saveas2, repsamsd1)
capture.output(print(repsamsd1.out), file="rep.sam.sd.txt", append=T)

######################################################
# Apply model to reference validation  ---------
######################################################
bugs.rv<-as.data.frame(subset(bugs.m, row.names(bugs.m) %in% row.names(stations.rv), drop = T))
stations.rv<-as.data.frame(subset(stations.rv, row.names(stations.rv) %in% row.names(bugs.rv), drop=T))

set.seed(10)
Scores.rv<-model.predict.RanFor.4.2(
       bugcal=bugs.rc,
       grps.final=dat.BG,
       preds.final=optvars.rfe, 
       ranfor.mod=rfe.out2,
       prednew=stations.rv,
       bugnew=(bugs.rv),
       Pc=Pcap,
       Cal.OOB=F)

saveas2<-paste0("Prednums_",x,"_Clusternums_",y,"_Pcap_",p,"_repsamsd_rv")
repsamsd2<-rep.sam.sd(occprb.cal=Scores.rv$Capture.Probs,Pc=Pcap);
repsamsd2.out<-paste(saveas2, repsamsd2)
capture.output(print(repsamsd2.out), file="rep.sam.sd.rv.txt", append=T)

####################################
# Apply model to intermediate sites -------
####################################
stations.int<-subset(stations, stations$SiteSetSample2=="Intermediate", drop = T)
bugs.int<-as.data.frame(subset(bugs.m, row.names(bugs.m) %in% row.names(stations.int), drop = T))
stations.int<-subset(stations.int, row.names(stations.int) %in% row.names(bugs.int))

set.seed(10)
Scores.int<-model.predict.RanFor.4.2(
  bugcal=bugs.rc,
  grps.final=dat.BG,
  preds.final=optvars.rfe, 
  ranfor.mod=rfe.out2,
  prednew=stations.int,
  bugnew=(bugs.int),
  Pc=Pcap,
  Cal.OOB=F)

####################################
# Apply model to stressed sites --------
####################################

stations.str<-subset(stations, stations$SiteSetSample2=="Stressed", drop = T)
stations.str<-stations.str[,optvars.rfe]
foo<-which(is.na(rowSums(stations.str)))
if (length(foo) > 0) { stations.str<-stations.str[-foo,] }

bugs.str<-as.data.frame(subset(bugs.m, row.names(bugs.m) %in% row.names(stations.str), drop = T))
foo<-which(colSums(bugs.str)==0)
if (length(foo) > 0) {bugs.str<-bugs.str[,-foo] }

stations.str<-subset(stations.str, row.names(stations.str) %in% row.names(bugs.str), drop = T)

set.seed(10)
Scores.str<-model.predict.RanFor.4.2(
  bugcal=bugs.rc,
  grps.final=dat.BG,
  preds.final=optvars.rfe, 
  ranfor.mod=rfe.out2,
  prednew=stations.str,
  bugnew=(bugs.str),
  Pc=Pcap,
  Cal.OOB=F)

####################################
# Apply model to ALL OTHERS sites --------
####################################

foo<-which(stations.orig$MostRecentSample=="n")
foo1<-stations.orig[foo,]

foo<-which(is.na(stations.orig$SiteSetSample2))
foo2<-stations.orig[foo,]

stations.all<-rbind(foo1, foo2)
stations.all<-stations.all[,optvars.rfe]
foo<-which(is.na(rowSums(stations.all)))
stations.all<-stations.all[-foo,]

bugs.all<-bugs.m
bugs.all<-ifelse(bugs.all>0,1,0)
foo<-which(colSums(bugs.all)==0)
if (length(foo) > 0) {bugs.all<-bugs.all[,-foo] }

foo<-which(row.names(stations.all) %in% row.names(bugs.all))
if (length(foo)>0) {stations.all <-stations.all[foo,]}
foo<-which(row.names(bugs.all) %in% row.names(stations.all))
if (length(foo)>0) {bugs.all <-bugs.all[foo,]}

foo<-which(colSums(bugs.all)==0)
bugs.all<-bugs.all[,-foo]
foo<-which(rowSums(bugs.all)==0)

bugs.all<-as.data.frame(bugs.all)
bugs.all<-bugs.all[order(row.names(bugs.all)),]
stations.all<-stations.all[order(row.names(stations.all)),]

row.names(bugs.all) == row.names(stations.all)

set.seed(10)
Scores.all<-model.predict.RanFor.4.2(
  bugcal=bugs.rc,
  grps.final=dat.BG,
  preds.final=optvars.rfe, 
  ranfor.mod=rfe.out2,
  prednew=stations.all,
  bugnew=(bugs.all),
  Pc=Pcap,
  Cal.OOB=F)

####################################
# Combined Scores --------
####################################

Scores.cal<-as.data.frame(Scores.cal); Scores.cal$Type<-"Reference"
Scores.rv<-as.data.frame(Scores.rv); Scores.rv$Type<-"Reference"
Scores.int<-as.data.frame(Scores.int); Scores.int$Type<-"Intermediate"
Scores.str<-as.data.frame(Scores.str); Scores.str$Type<-"Stressed"
Scores.all<-as.data.frame(Scores.all); Scores.all$Type<-"Other"

Scores.cat<-rbind(Scores.cal, Scores.rv, Scores.int, Scores.str, Scores.all)

type<-data.frame(row.names(stations.orig),stations.orig$SiteSetSample2, stations.orig$RefCalVal, stations.orig$StrCalVal )
row.names(type)<-type$row.names.stations.orig.

#foo<-which(colnames(Scores.cat)=="Type")
#Scores.cat<-Scores.cat[,-foo]

source("~/Documents/R/bin/cbind.na.R")
Scores.cat<-cbind.na(Scores.cat, type)

####################################
# Plot ref/int/str --------
####################################

# plot O/E by Ref/Int/Stressed 
saveas<-paste0("Prednums_",x,"_Clusternums_",y,"_Pcap_", p,"_OEbyType.pdf")
ggplot(Scores.cat, aes(y=Scores.cat$OE.scores.OoverE, x=Scores.cat$Type)) +
  geom_boxplot(aes(fill=Scores.cat$Type)) + 
  geom_jitter(cex=0.02) +
  ylim(0,1.5) + scale_fill_manual(values=c("#FFCC66", "#339966", "#CC0000")) + 
  scale_x_discrete(limits=c("Reference", "Intermediate", "Stressed")) +
  geom_hline(yintercept = 1) + theme_bw() + ylab("O/E") + xlab("")
if (plotz == T) { ggsave(file=saveas, width = 6, height = 6) }


########################################
# Stats -----------
########################################

cat1<-Scores.cat
cat1<-subset(cat1, Type!="Other")

#calculating stats by ref/int/str
saveas<-paste0("Prednums_",x,"_Clusternums_",y,"_Pcap_", p,"_anova")
set.seed(10)
fit1<-aov(cat1$OE.scores.OoverE~cat1$Type); summary(fit1)
anova.P<-summary(fit1)[[1]][["Pr(>F)"]][[1]]
anova.F<-summary(fit1)[[1]][["F value"]][[1]]
anova.out<-paste(saveas, "Pv", anova.P, "Fv", anova.F) 
capture.output(anova.out, file="anova.out.txt", append=T)

saveas<-paste0("Prednums_",x,"_Clusternums_",y,"_Pcap_", p,"_ttest")
set.seed(10)
ttest1<-pairwise.t.test(cat1$OE.scores.OoverE, cat1$Type)
ttest.RefInt<-ttest1$p.value[1]
ttest.RefStr<-ttest1$p.value[2]
ttest.IntStr<-ttest1$p.value[4]
ttest.RefStr.t<-ttest1$p.value[2]
ttest.out<-paste(saveas, "RefInt", ttest.RefInt, "RefStr", ttest.RefStr, "IntStr", ttest.IntStr) 
capture.output(ttest.out, file="pw.ttest.out.txt", append=T)

cat1.sub<-subset(cat1, Type!="Intermediate")
set.seed(10)
ttest1<-t.test(cat1.sub$OE.scores.OoverE~cat1.sub$Type)
ttest.RefStr<-ttest1$statistic[1]
ttest.out<-paste(saveas,"RefStr", ttest.RefStr) 
capture.output(ttest.out, file="ttest.out.txt", append=T)
  

  }}}



########## Export files -----------

if (exportfiles==T) { 
  
  if(diatoms==T) {write.csv(stations.rc, "diatoms.stations.rc.csv")}
  if(diatoms==T) {write.csv(bugs.rc, "diatoms.bugs.rc.csv")}
  if(diatoms==T) {write.csv(stations, "diatoms.stations.csv")}
  if(diatoms==T) {write.csv(bugs.m, "diatoms.bugs.m.csv")}
  if (diatoms==T) { if (genuslev==F) { write.csv(bugs.m, "diatoms.bugs.m.species.csv")}}
  
  if(sba==T) {write.csv(stations.rc, "sba.stations.rc.csv")}
  if(sba==T) {write.csv(bugs.rc, "sba.bugs.rc.csv")}
  if(sba==T) {write.csv(stations, "sba.stations.csv")}
  if(sba==T) {write.csv(bugs.m, "sba.bugs.m.csv")}
  if (sba==T) { if (genuslev==F) { write.csv(bugs.m, "sba.bugs.m.species.csv")}}
  
  if(combo==T) {write.csv(stations.rc, "combo.stations.rc.csv")}
  if(combo==T) {write.csv(bugs.rc, "combo.bugs.rc.csv")}
  if(combo==T) {write.csv(stations, "combo.stations.csv")}
  if(combo==T) {write.csv(bugs.m, "combo.bugs.m.csv")}
  if(combo==T) { if (genuslev==F) { write.csv(bugs.m, "combo.bugs.m.species.csv")}}
  
  
  
  write.csv(Scores.cat, "Scores.all.csv")

  
}


# import the RFE scores into df
rfe.final<-read.table("RFEscores.txt", sep = " ", stringsAsFactors = F)
rfe.final <- data.frame(do.call('rbind', strsplit(as.character(rfe.final$V2),' ',fixed=TRUE)))
colnames(rfe.final)<-c("combo", "model.mean", "model.stdev", "null.mean", "null.stdev", "dunno")
combo<-data.frame(do.call('rbind', strsplit(as.character(rfe.final$combo),'_',fixed=TRUE)))                               
colnames(combo)<-c("a", "prednums", "b", "clusternums", "c", "pcap", "d")
del<-c("a", "b", "c", "d")
foo<-which(colnames(combo) %in% del)
combo<-combo[,-foo]
rfe.final<-cbind(combo, rfe.final)                               
write.csv(rfe.final, "rfe.all.csv")

# import the ANOVA scores into df
anova.final<-read.table("anova.out.txt", sep = " ", stringsAsFactors = F)
anova.final <- data.frame(do.call('rbind', strsplit(as.character(anova.final$V2),' ',fixed=TRUE)))
colnames(anova.final)<-c("combo", "Pv", "pvalue", "Fv", "Fvalue")
combo<-data.frame(do.call('rbind', strsplit(as.character(anova.final$combo),'_',fixed=TRUE)))                               
colnames(combo)<-c("a", "prednums", "b", "clusternums", "c", "pcap", "d")
del<-c("a", "b", "c", "d")
foo<-which(colnames(combo) %in% del)
combo<-combo[,-foo]
anova.final<-cbind(combo, anova.final)                               
write.csv(anova.final, "anova.all.csv")

##### 
# make compare file 

foo1<-read.csv("rfe.all.csv")
foo2<-read.csv("anova.all.csv")

foo1<-subset(foo1, foo1$null.mean==1)
foo3<-cbind(foo1, foo2$Fvalue)
foo3$diff<-(foo3$model.stdev - foo3$null.stdev)
foo3<-subset(foo3, foo3$diff<0)

write.csv(foo3, "compare.csv")

