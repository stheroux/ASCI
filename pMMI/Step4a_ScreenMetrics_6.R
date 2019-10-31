# screen metrics

setwd("~/Documents/R/ASCI/pMMI/")
source("~/Documents/R/bin/cbind.na.R")
library(ggplot2)
source("~/Documents/R/bin/cbind.na.R")


diatoms = F # pick one
sba = F
hybrid = T

if (diatoms == T) {assemblage <- "diatoms"}
if (sba == T) {assemblage <- "sba"}
if (hybrid == T) {assemblage <- "hybrid"}

if (diatoms == T) {setwd("~/Documents/R/ASCI/pMMI/diatoms/")}
if (sba == T) {setwd("~/Documents/R/ASCI/pMMI/sba/")}
if (hybrid == T) {setwd("~/Documents/R/ASCI/pMMI/hybrid/")}


#################################
# Import data ------------------
#################################

combined<-read.csv(paste0(assemblage, ".combined.scaled.csv"), row.names=1, header=T)

# Import stations data 
stations<-read.csv('~/Documents/R/ASCI/DATA/algae.site.data.07202019.csv')
row.names(stations)<-stations$SampleID_old
stations<-stations[order(row.names(stations)),]
stations<-subset(stations, row.names(stations) %in% row.names(combined), drop=T)
foo<-sapply(stations, is.numeric)
stations.num<-stations[,foo] # only numeric columns 

row.names(combined) == row.names(stations)
row.names(combined) == row.names(stations.num)

z<-length(row.names(combined))
foo<-which (colSums(!is.na(combined) == 0) > (z-400)) # find columns with 2000+ NAs 
if (length(foo) > 0 ) {combined<-combined[,-foo]}

metrics.list<-colnames(combined)

cal<-subset(stations, stations$RefCalVal=="Cal")
val<-subset(stations, stations$RefCalVal=="Val")

stations.str<-subset(stations, stations$SiteSetSample2=="Stressed")
stations.rc<-subset(stations, row.names(stations) %in% row.names(cal))
stations.rc<-droplevels(stations.rc)

combined.str<-subset(combined, row.names(combined) %in% row.names(stations.str))
combined.rc<-subset(combined, row.names(combined) %in% row.names(stations.rc))

##########################################################
# Performance of individual metrics -----------------------
##########################################################
# ANOVA across PSA6c regions for ref stations ----------
# < 2 (or 3? )
Pv.list<-list()
Fv.list<-list()
metrics1<-list()

  for (i in metrics.list) {
    fit<-aov(combined.rc[[i]]~stations.rc$PSA6ced)
    Pv<-summary(fit)[[1]][["Pr(>F)"]][[1]]
    Fv<-summary(fit)[[1]][["F value"]][[1]]
    Pv.list[[i]]<-Pv
    Fv.list[[i]]<-Fv
  } 

anova.psa<-cbind(Pv.list, Fv.list)
#write.csv(anova.psa, paste0(assemblage,".anova.PSA6c.ref.csv"))

i<-"prop.ind.most.tol"
ggplot(combined.rc) + 
  geom_boxplot(aes(x=stations.rc$PSA6ced, y=combined.rc[[i]])) +
  geom_jitter(aes(x=stations.rc$PSA6ced, y=combined.rc[[i]]), cex=0.4)


# ANOVA across ref/int/str ------------------

Pv.list<-list()
Fv.list<-list()

for (i in metrics.list) {
  foo<-as.matrix(combined[,i])
  fit<-aov(foo~stations$SiteSetSample2)
  Pv<-summary(fit)[[1]][["Pr(>F)"]][[1]]
  Fv<-summary(fit)[[1]][["F value"]][[1]]
  Pv.list[[i]]<-Pv
  Fv.list[[i]]<-Fv
}

anova.refintstr<-cbind(Pv.list, Fv.list)
#write.csv(anova.refintstr, file=paste0(assemblage,".anova.ref.stress.csv"))

# ttest across ref/str ------------------------------

tt<-list()
ttp<-list()

for (i in metrics.list) {
tt.out<-t.test(combined.rc[,i], combined.str[,i])
tt[i]<-tt.out$statistic
ttp[i]<-tt.out$p.value
}

tt.out2<-cbind(tt,ttp)
#write.csv(tt.out2, paste0(assemblage,".ttest.ref.str.csv"))

# Frequence of 1's and 0's in the ref dataset -------------------------
# MaxFreqZero<-.33
# MaxCount<-5

library(plyr)

FreqZeroRefCal<-list()
for (i in metrics.list) { # i<-"prop.spp.BCG5"
  length1<-length(combined.rc[[i]])
  foo<-sum(combined.rc[i]==0)
  if (is.na(foo)) { FreqZeroRefCal[i]<- 0 }
  if (is.numeric(foo)) {FreqZeroRefCal[i]<-(foo/length1)*100 }
}

FreqOneRefCal<-list()
for (i in metrics.list) { # i<-"prop.spp.BCG5"
  length1<-length(combined.rc[[i]])
  foo<-sum(combined.rc[i]==1)
  if (is.na(foo)) { FreqOneRefCal[i]<- 0 }
  if (is.numeric(foo)) {FreqOneRefCal[i]<-(foo/length1)*100 }
}

Freqs<-cbind(FreqZeroRefCal, FreqOneRefCal)
#write.csv(Freqs, paste0(assemblage,".Freqs.csv")) 

# Rafi SN --------------------------------------------------------------
# F value less than 3
# variance across all sites, variance at repeat samplings

DF<-cbind(stations$StationCode, combined)
StationCode.mean<-ddply(DF, .(`stations$StationCode`), colwise(mean))
StationCode.sd<-ddply(DF, .(`stations$StationCode`), colwise(sd))
StationCode.var<-ddply(DF, .(`stations$StationCode`), colwise(var))
SN.rafi<-cbind(StationCode.mean, StationCode.sd, StationCode.var)
#write.csv(SN.rafi, paste0(assemblage,".Rafi.SN.csv"))

# Rafi ANOVA, variance within station code, want F < 3 -------------------------------
Pv.list<-list()
Fv.list<-list()

refz<-droplevels(subset(stations, stations$SiteSetSample2=="Reference"))
refz.stations.codes<-unique(refz$StationCode)
combined.rcz<-cbind(combined, stations$StationCode)
foo<-which(combined.rcz$`stations$StationCode` %in% refz.stations.codes)
combined.rcz<-combined.rcz[foo,]

for (i in metrics.list) { 
  fit<-aov(combined.rcz[[i]]~combined.rcz$`stations$StationCode`)
  Pv<-summary(fit)[[1]][["Pr(>F)"]][[1]]
  Fv<-summary(fit)[[1]][["F value"]][[1]]
  Pv.list[[i]]<-Pv
  Fv.list[[i]]<-Fv
    }
rafi.anova.out<-cbind(Pv.list, Fv.list)
#write.csv(rafi.anova.out, file=paste0(assemblage,".anova.within.stationcode.csv"))


# Range test (Stevenson 2013) --------------------------------------------
# the median value of a candidate metric was > 0 in either reference or disturbed sites
Range.ref<-list()
Range.str<-list()
for (i in metrics.list) { 
  a <- median(combined.rc[[i]], na.rm=T) 
  b <- median(combined.str[[i]], na.rm=T) 
  Range.ref[i] <- a
  Range.str[i] <- b 
  }
Range<-cbind(Range.ref, Range.str)
#write.csv(Range, file = paste0(assemblage,".Range.Stevenson.csv"))

freq.max<-list()
for (i in metrics.list) { # i<-"prop.spp.IndicatorClass_Ref_raw"
  max1<-max(combined[[i]])
  freq1<-length(which(combined[[i]]==max1))
  freq2<-freq1/length(combined[[i]])*100
  freq.max[i] <- freq2
}


Mode = function(x){ 
  ta = table(x)
  tam = max(ta)
  if (all(ta == tam))
    mod = NA
  else
    if(is.numeric(x))
      mod = as.numeric(names(ta)[ta == tam])
  else
    mod = names(ta)[ta == tam]
  return(mod)
}

freq.mode<-list()
for (i in metrics.list) { # i<-"prop.spp.IndicatorClass_Ref_raw"
  max1<-Mode(combined[[i]])
  freq1<-length(which(combined[[i]]==max1))
  freq2<-freq1/length(combined[[i]])*100
  freq.mode[i] <- freq2
}


# signal to noise (Stoddard) ----------------------------------------
# variance across all sites /  variance at repeat samplings 
# > 2 is best , periphyton 1 or 1.5 okay 
SN<-list()
for (i in metrics.list) { 
SN[i] <- var(combined[[i]], na.rm=T) / mean(StationCode.var[[i]], na.rm=T) }
#write.csv(SN, file=paste0(assemblage,".SN.stoddard.csv"))

# Now calculate lm to disturbance gradients ----------------------
# from Cao et al 2007
#list3<- c("DayOfYear", "Elevation", "XSLOPE", "NHDSLOPE", "New_Lat", "New_Long", "LST32AVE", "MEANP_WS", "PPT_00_09", "CondQR50", "TMAX_WS", "XWD_WS", 
#  "KFCT_AVE", "BDH_AVE", "PRMH_AVE" )
list3<-c("MaxOfW1_HALL", "Ag_2000_5K", "URBAN_2000_5K", "CODE_21_2000_5K","RoadDens_5K",
               "PAVED_INT_1K","PerManMade_WS")

r2.list<-list()
slope.list<-list()

for (x in metrics.list) { 
  for (y in list3) {  #or list2 for all vars
    NAsum<-sum(is.na(stations.num[[y]]))
    z<-length(row.names(stations.num))
    if((NAsum/z) < .20) {  
      fit<-lm(stations.num[[y]]~combined[[x]], na.action = "na.exclude")
      r2 <- format(summary(fit)$adj.r.squared, digits=3) 
      slope <-format(coef(fit)[2], digits=3)
      saveas<-paste(x,y)
      r2.list[[saveas]]<-r2
      slope.list[[saveas]]<-slope
    }}}

lm.out<-as.matrix(cbind(r2.list, slope.list))
#write.csv(lm.out, "lm.out.csv")

if (diatoms==T) {write.csv(lm.out, "diatoms.lm.out.gradients.csv")}
if (sba==T) {write.csv(lm.out, "sba.lm.out.gradients.csv")}
if (hybrid==T) {write.csv(lm.out, "hybrid.lm.out.gradients.csv")}

## try a distrubance PCA like Fetscher 2016 --------------
# disturb.vars<-c("MaxOfW1_HALL", "Ag_2000_5K", "URBAN_2000_5K", "CODE_21_2000_5K","RoadDens_5K",
#                "PAVED_INT_1K","PerManMade_WS")
#stations.disturb.vars<-stations[,disturb.vars]
#foo<-which(is.na(rowSums(stations.disturb.vars)))
#rda1<-rda(stations.disturb.vars[-foo,], cor=T)
#summary(rda1)

# SD
SD<-list()
#cal<-droplevels(subset(calval, calval$CalVal=="Cal"))
combined.rc<-droplevels(subset(combined, row.names(combined) %in% row.names(cal)))

for (i in metrics.list) { 
  SD[i] <- sd(combined.rc[[i]], na.rm=T) }




# combine results ------------------------------------------------------------
results<-cbind(anova.psa,
               anova.refintstr,
               tt.out2,
               FreqZeroRefCal,
               FreqOneRefCal,
               freq.max,
               freq.mode,
               #SN.rafi,
               rafi.anova.out,
               Range,
               SN, 
               SD)
results<-as.data.frame(results)

colnames(results)<-c("anova.psa.pv", "anova.psa.fv",
                     "anova.refintstr.pv","anova.refintstr.fv",
                     "tt.out2.t", "tt.out2.p",
                     "FreqZero", "FreqOne", "Freq.Max", "Freq.Mode",
                     "anova.rafi.pv", "anova.rafi.fv",
                     "Range.ref", "Range.str",
                     "SN", "SD")
results<-as.matrix(results)
write.csv(results, file=paste0(assemblage,".results.csv"))
results<-as.data.frame(results)

results.PF<-data.frame(row.names(results))
results.PF$anova.psa <-ifelse (results$anova.psa.fv > 10, "fail", "pass" )
results.PF$anova.refintstr <- ifelse (results$anova.refintstr.fv > 10, "pass", "fail") 
foo<-as.numeric(results$tt.out2.t)
results.PF$tt.out2 <- ifelse ( sqrt(foo^2) < 10, "fail", "pass")
results.PF$FreqZero <- ifelse (results$FreqZero > 33, "fail", "pass")
results.PF$FreqOne <- ifelse (results$FreqOne > 33, "fail", "pass")                         
#results.PF$rafi.anova.out <- ifelse (results$anova.rafi.fv > 1, "pass", "fail")
results.PF$Range.ref <- ifelse (results$Range.ref > 0, "pass", "fail")
results.PF$Range.str <- ifelse (results$Range.str > 0, "pass", "fail")
results.PF$SN <- ifelse (results$SN > 1, "pass", "fail")
results.PF$SD<-ifelse(results$SD<0.25, "pass", "fail")
#write.csv(results.PF, file=paste0(assemblage,".results.PF.csv"))

results.PF.sum<-data.frame(row.names(results))
results.PF.sum$Pass<-apply(results.PF, 1, function(x) length(which(x=="pass")))
results.PF.sum$Fail<-apply(results.PF, 1, function(x) length(which(x=="fail")))

#write.csv(results.PF.sum, file=paste0(assemblage,".results.PF.sum.csv") )

results2<-as.matrix(cbind(results, results.PF, results.PF.sum))
results2<-as.data.frame(results2)



##################################
# add rsq info --------------
##################################

rsq<-read.csv(file = paste0(assemblage, ".refcal.rsq.csv"), header=T, row.names=2)
#rsq<-rsq[,-1]
#rsq<-as.data.frame(rsq)
#foo<-which(rsq$rsq>0.001)
#rsq<-rsq[foo,]
results3<-cbind.na(results2, rsq)
#results3<-merge(results2,rsq,by="row.names",all.x=TRUE)
del1<-which(results3$rsq<0.2)
results3<-results3[-del1,]
write.csv(as.matrix(results3), file=paste0(assemblage, ".results3.csv") )

##################################
# add rsq info --------------
##################################

rsq<-read.csv(file = paste0(assemblage, ".refcal.rsq.csv"), header=T, row.names=2)
#rsq<-rsq[,-1]
#rsq<-as.data.frame(rsq)
#foo<-which(rsq$rsq>0.001)
#rsq<-rsq[foo,]
results3<-cbind.na(results2, rsq)
#results3<-merge(results2,rsq,by="row.names",all.x=TRUE)
del1<-which(results3$rsq<0.2)
results3<-results3[-del1,]
write.csv(as.matrix(results3), file=paste0(assemblage, ".results3.csv") )






