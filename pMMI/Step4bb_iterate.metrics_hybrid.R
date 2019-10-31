#find best custom combination of metrics for pMMI 
#step 1:  import your data from Step4b so it is in your global environment 

setwd("~/Documents/R/ASCI/pMMI/hybrid/")

source("~/Documents/R/bin/cbind.na.R")
metrics.cats<-read.csv("~/Documents/R/ASCI/DATA/metric.categories.csv", row.names=1)

# screen metrics
results3<-read.csv("hybrid.results3.csv", row.names=1, stringsAsFactors = F)
results3<-as.data.frame(results3)
foo<-subset(results3, results3$anova.psa.fv<10)
foo<-subset(foo, anova.refintstr.fv>14)
foo<-subset(foo, Range.ref.1=="pass")
foo<-subset(foo, Range.str.1=="pass")
foo<-subset(foo, tt.out2.t>15)
foo<-subset(foo, SN>2)
foo$Freq.Max<-as.numeric(as.character(foo$Freq.Max))
foo$Freq.Mode<-as.numeric(as.character(foo$Freq.Mode))
del<-which(foo$Freq.Max>10)
if (length(foo)>0) {foo<-foo[-del,]}
del<-which(foo$Freq.Mode>20)
if (length(foo)>0) {foo<-foo[-del,]}
foo$FreqZero<-as.numeric(as.character(foo$FreqZero))
foo$FreqOne<-as.numeric(as.character(foo$FreqOne))
foo<-subset(foo, FreqZero<10)
foo<-subset(foo, FreqOne<10)

mypick<-c(NULL)
mypick<-row.names(foo)

win.metrics<-mypick
##########################################
# iterate winning metrics 
##########################################

MMI.min<-4
MMI.max<-6

MMI.all.subsets<-unlist(
  lapply(MMI.min:MMI.max, function(i)
    combn(x=win.metrics, m=i, simplify=F)),recursive=F)

#### run parallel 

source("~/Documents/R/ASCI/pMMI/1screen.metrics.R")
library(parallel)
rc<-subset(calval, calval$RefCalVal=="Cal")

#foo<-MMI.all.subsets[1:3]
#sapply(foo, function(x) fun=screen.metrics(x)) # make sure code is working

Pv.list<-list()
Fv.list<-list()
SN.list<- list()
sd1.list<-list()
psa.anova.list<-list()

# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores)
base<-screen.metrics
clusterExport(cl, c("base", "combined.sc", 
                    "stations", "calval", "rc",
                    "Pv.list",
                    "Fv.list",
                    "SN.list",
                    "sd1.list",
                    "psa.anova.list"
                    ))

test.out<-parSapply(cl, MMI.all.subsets,
          function(x)
            fun=base(x))

stopCluster(cl)

out1<-do.call(rbind, strsplit(as.character(test.out), ","))
out1<-as.data.frame(out1)
out1$CAT<-as.vector(MMI.all.subsets)
out1$V1<-as.numeric(as.character(out1$V1))
out1$V2<-as.numeric(as.character(out1$V2))
out1$V3<-as.numeric(as.character(out1$V3))
out1<-subset(out1, out1$V1>100)
out1<-subset(out1, out1$V2<0.18)
out1<-subset(out1, out1$V3<3)
write.csv(as.matrix(out1), "out1.csv")

list1<-list()

for (i in win.metrics){
  foo1<-length(grep(i, out1$CAT))
  print(paste(i,foo1))
  list1[i]<-foo1
  }

list1<-as.data.frame(t(t(list1)))

write.csv(as.matrix(list1), "iterate.out.csv")


##########################################
# looking at results
##########################################

Pv.list<-list()
Fv.list<-list()
SN.list<- list()
sd1.list<-list()
psa.anova.list<-list()


#### now the winners 

win.metrics<-c(  "prop.spp.BCG4",
                 "Salinity.BF.richness",
                 "cnt.spp.IndicatorClass_TP_low_raw",
                 "prop.spp.IndicatorClass_DOC_high_raw",
                 "OxyRed.DO_30.richness",
                 "cnt.spp.IndicatorClass_TN_low_raw",
                NULL) # most often in tops 



MMI.min<-4
MMI.max<-6

MMI.all.subsets<-unlist(
  lapply(MMI.min:MMI.max, function(i)
    combn(x=win.metrics, m=i, simplify=F)),recursive=F)

cl <- makeCluster(no_cores)
base<-screen.metrics
clusterExport(cl, c("base", "combined.sc", 
                    "stations", "calval", "rc",
                    "Pv.list",
                    "Fv.list",
                    "SN.list",
                    "sd1.list",
                    "psa.anova.list"
))

test.out<-parSapply(cl, MMI.all.subsets,
          function(x)
            fun=base(x))

stopCluster(cl)

out1<-do.call(rbind, strsplit(as.character(test.out), ","))
out1<-as.data.frame(out1)
out1$CAT<-as.vector(MMI.all.subsets)
out1$V1<-as.numeric(as.character(out1$V1))
out1$V2<-as.numeric(as.character(out1$V2))
out1$V3<-as.numeric(as.character(out1$V3))
out1<-subset(out1, out1$V1>100)
out1<-subset(out1, out1$V2<0.18)
out1<-subset(out1, out1$V3<3)
write.csv(as.matrix(out1), "out1.best.csv")

list1<-list()

for (i in win.metrics){
  foo1<-length(grep(i, out1$CAT))
  print(paste(i,foo1))
  list1[i]<-foo1
}

list1<-as.data.frame(t(t(list1)))

write.csv(as.matrix(list1), "iterate.out.best.csv")

### 

win.metrics<-c("prop.spp.BCG4",
                "Salinity.BF.richness",
               "prop.spp.IndicatorClass_DOC_high_raw",
               "OxyRed.DO_30.richness",
               NULL)

write.csv(win.metrics, "win.metrics.csv")


win.metrics.raw<-c("prop.spp.BCG4_raw",
                   "Salinity.BF.richness_raw",
                   "prop.spp.IndicatorClass_DOC_high_raw",
                   "OxyRed.DO_30.richness_raw",
                   NULL)
write.csv(win.metrics.raw, "win.metrics.raw.csv")

# export your random forest models (update aug 2019)




