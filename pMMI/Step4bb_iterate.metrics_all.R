#find best custom combination of metrics for pMMI 
#step 1:  import your data from Step4b so it is in your global environment 

setwd(paste0("~/Documents/R/ASCI/pMMI_newnew/", assemblage, "/"))

# screen metrics
results3<-read.csv(paste0(assemblage,".results3.csv"), row.names=1, stringsAsFactors = F)
results3<-as.data.frame(results3)
foo<-subset(results3, results3$anova.psa.fv<3)
foo<-subset(foo, tt.out2.t>10)
foo<-subset(foo, FreqZero<33)
foo<-subset(foo, FreqOne<33)
foo<-subset(foo, Range.ref >0)
foo<-subset(foo, Range.str >0)
foo<-subset(foo, SN>2) # or 1.5? 
foo<-subset(foo, foo$anova.rafi.fv<3) # or 1.5? <- change to anova on repeat samplings
mypick<-c(NULL)
mypick<-row.names(foo)

win.metrics<-mypick
##########################################
# iterate winning metrics 
##########################################

MMI.min<-4
MMI.max<-ifelse(length(win.metrics)>8,8,max(length(win.metrics)))

MMI.all.subsets<-unlist(
  lapply(MMI.min:MMI.max, function(i)
    combn(x=win.metrics, m=i, simplify=F)),recursive=F)

#### run parallel 

source("~/Documents/R/ASCI/pMMI_newnew/1screen.metrics.R")
library(parallel)
rc<-subset(calval, calval$RefCalVal=="Cal")

Pv.list<-list()
Fv.list<-list()
SN.list<- list()
sd1.list<-list()
psa.anova.list<-list()

foo<-MMI.all.subsets[1:3]
sapply(foo, function(x) fun=screen.metrics(x)) # make sure code is working

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

foo<-read.csv("iterate.out.csv", row.names=1)
foo<-subset(foo, foo$x %in% "x")

win.metrics<-row.names(foo)

MMI.min<-4
MMI.max<-ifelse(length(win.metrics)>8,8,max(length(win.metrics)))

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
write.csv(as.matrix(out1), "out1.best.csv")

list1<-list()

for (i in win.metrics){
  foo1<-length(grep(i, out1$CAT))
  print(paste(i,foo1))
  list1[i]<-foo1
}

list1<-as.data.frame(t(t(list1)))

write.csv(as.matrix(list1), "iterate.out.best.csv")


