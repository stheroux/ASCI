# screen metrics

setwd(paste0("/Users/susannat/Documents/R/ASCI/pMMI_newnew/",assemblage))


#################################
# Import data ------------------
#################################

combined<-read.csv("combined.not.scaled.csv", row.names=1, header=T)


# Import stations data 
stations<-read.csv('~/Documents/R/ASCI/DATA/algae.site.data.04052019.csv')
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

i<-"prop.spp.BCG45"
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

########################################################################################
# Frequence of one's and zero's in the ref dataset -------------------------
########################################################################################
# MaxFreqZero<-.33
# MaxCount<-5

library(plyr)

FreqZeroRefCal<-list()
for (i in metrics.list) { # i<-"prop.spp.Salinity.BF"
  length1<-length(combined.rc[[i]])
  foo<-sum(combined.rc[i]==0, na.rm = T)
  if (is.na(foo)) { FreqZeroRefCal[i]<- 0 }
  if (is.numeric(foo)) {FreqZeroRefCal[i]<-(foo/length1)*100 }
}

FreqOneRefCal<-list()
for (i in metrics.list) { # i<-"Salinity.BF.richness"
  length1<-length(combined.rc[[i]])
  foo<-sum(combined.rc[i]==1, na.rm = T)
  if (is.na(foo)) { FreqOneRefCal[i]<- 0 }
  if (is.numeric(foo)) {FreqOneRefCal[i]<-(foo/length1)*100 }
}

Freqs<-cbind(FreqZeroRefCal, FreqOneRefCal)

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
               #freq.max,
               #freq.mode,
               #SN.rafi,
               rafi.anova.out,
               Range,
               SN, 
               SD)
results<-as.data.frame(results)

colnames(results)<-c("anova.psa.pv", "anova.psa.fv",
                     "anova.refintstr.pv","anova.refintstr.fv",
                     "tt.out2.t", "tt.out2.p",
                     "FreqZero", "FreqOne",
                     "anova.rafi.pv", "anova.rafi.fv",
                     "Range.ref", "Range.str",
                     "SN", "SD")
results<-as.matrix(results)
write.csv(results, file=paste0(assemblage,".results.csv"))
results<-as.data.frame(results)

results2<-results
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






